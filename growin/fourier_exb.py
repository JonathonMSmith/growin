#time for to be obtaining vertical drifts from satellite data using pysat

import warnings
import numpy as np
from scipy.optimize import curve_fit
import xarray as xr
from datetime import datetime
from pysat import Instrument
from pysat.ssnl.avg import median2D as med2D

def make_fourier(na, nb):
    """ The function for the curve fit

    Parameters
    ----------
    na: (int)
        number of cosine terms/coefficients
    nb: (int)
        number of sin terms/coefficients
    """
    def fourier(x, *a):
        ret = a[0]
        for deg in range(0, na):
            ret += a[deg + 1] * np.cos((deg + 1) * np.pi * x / 12)
        for deg in range(na, na+nb):
            ret += a[deg + 1] * np.sin((deg - na + 1) * np.pi * x / 12)
        return ret
    return fourier

def fourier_fit(local_times, median_drifts, num_co):
    """ Here the terms in the fourier fit are actually determined

    Parameters
    ----------
    local_times : (array-like)
        xdim for fit; local time values
    median_drifts : (array-like)
        ydim for fit; median drift values from data
    num_co : (int)
        'number of coefficients) how many sin/cosine pairs for the fit
    """
    exb_drifts = np.zeros((num_co, 2))
    ind, = np.where(~np.isnan(median_drifts))
    if ind.size < num_co*2+1:
        warnings.warn('not enough viable drift data, '
                      'returning zero value \"flat fit\"', Warning)
        return 0, exb_drifts
    #popt contains the coeficients. First ten are cosines, second ten are sins
    popt, pcov = curve_fit(make_fourier(num_co, num_co), local_times[ind], 
                           median_drifts[ind], [0.0]*(num_co*2+1))
    #format the coefficients for input ito the SAMI2 model
    #the shape is np.zeroes((10,2))
    ve01 = popt[0]
    for n in range(1, num_co*2):
        i = (n - 1) % num_co
        j = int((n - 1) / num_co)
        exb_drifts[i, j] = popt[n]

    return ve01, exb_drifts

class DriftInstrument(Instrument):
    """Class that inherits from a pysat instrument that will have an attribute
       corresponding to the median drifts, their deviation, and the fourier
       curve fit for use in the SAMI2 model
    """
    def fit_drifts(self, start, stop, coords):
        """This gets the median drifts and the coefficients and puts them in
           an xarray Dataset. The drifts are first obtained using the pysat
           function pysat.ssnl.avg.median2D and then the fits are performed
           on the output and all of it is placed in xarray.DataArrays that
           are then combined into a data set for this Year/Season.
           The longitude coordinate is just the longitude values but will/can
           be later specified as names if desired.

        Parameters
        ----------
        start : (datetime)
            start date for median2d
        stop : (datetime)
            stop date for median2d
        """
        self.bounds = (start, stop)
        #compute medians for each longitude sector in hourly LT bins
        drift_medians = med2D(self, self.slt_bins, 'slt',
                              self.lon_bins, 'glon', [self.drift_key],
                              auto_bin=False)
        #create dataArrays for drifts, deviations, and fourier fits
        tmp_drift = xr.DataArray(drift_medians[self.drift_key]['median'],
                                 [coords[0], coords[1]])
        tmp_dev = xr.DataArray(drift_medians[self.drift_key]['avg_abs_dev'],
                               [coords[0], coords[1]])
        fit_out = [fourier_fit(self.slt_bins[:-1], x, self.num_co)
                   for x in tmp_drift.values]
        ve01 = xr.DataArray([i[0] for i in fit_out], [coords[0]])
        coef = xr.DataArray([i[1] for i in fit_out], [coords[0], coords[2],
                                                      coords[3]])
        #return a Dataset containing all of the dataArrays
        return xr.Dataset({'drifts': tmp_drift, 'dev': tmp_dev, 've01': ve01,
                           'coefficients': coef})

    def compile_drifts(self):
        """Builds larger Dataset containg all of the drifts and fits over the
           entire user-specified annual range.
        """
        year_drift = []
        #create coordinate tuples for xarray
        coords = [('longitude', self.lon_bins[:-1]),
                  ('slt', self.slt_bins[:-1]),
                  ('term', np.linspace(1, self.num_co, self.num_co)),
                  ('func', ['cos', 'sin'])]
        for year in range(self.start_year, self.stop_year+1):
            #iterate through the four seasons
            tmpyear = year
            sea_drift = []
            for season in range(len(self.season_bins)-1):
                #constant to subtract from year in start date to maintain
                #contiguous December solstices
                start = datetime(tmpyear, self.season_bins[season], 1)
                if self.season_bins[season] > self.season_bins[season+1]:
                    tmpyear += 1
                stop = datetime(tmpyear, self.season_bins[season+1], 1)

                sea_drift.append(self.fit_drifts(start, stop, coords))
                #store all of these bins for all years
                #change this to be setting xarray coord names
            year_drift.append(xr.concat(sea_drift, 'season'))
        drift_arr = xr.concat(year_drift, 'year')
        if self.season_names is not None:
            drift_arr = drift_arr.assign_coords(season=(self.season_names))
        if self.zone_labels is not None:
            drift_arr = drift_arr.assign_coords(longitude=(self.zone_labels))

        year_coord = (list(range(self.start_year, self.stop_year+1)))
        drift_arr = drift_arr.assign_coords(year=year_coord)
        return drift_arr

    def get_drifts(self, drift_key=None, num_co=10,
                   lon_bins=np.linspace(0, 10, 2, dtype=int),
                   slt_bins=np.linspace(0, 24, 49),
                   season_bins=np.linspace(1, 2, 2, dtype=int),
                   season_names=None, zone_labels=None,
                   start_year=None, stop_year=None):
        """The Big Function. This is the main function that generates all of
           the data and organizes it for the drift attribute in the object.

        Parameters
        ----------
        drift_key : (string)
            instrument key for vertical drift
        num_co : (int)
            how many sin/cosine pairs for the fitdownload
        slt_bins : (array-like)
            bin edges in local time
        lon_bins : (array-like)
            bin edges in longitude
        season_bins : (array-like)
            months used as edges of the season "bins"
        season_names : (array-like)
            string names for season xarray coords
        zone_labels : (array-like)
            string names for longitude xarray coords
        start_year : (int)
            the year to start the seasonal averaging
        stop_year : (int)
            the year to stop the seasonal averaging
        """
        if drift_key is None:
            raise AttributeError('drift_key must be specified')
        #assign variables to the object
        self.drift_key = drift_key
        self.num_co = num_co
        self.lon_bins = lon_bins
        self.slt_bins = slt_bins
        self.season_bins = season_bins
        self.season_names = season_names
        self.zone_labels = zone_labels
        self.start_year = start_year
        self.stop_year = stop_year
        #sets the available limits as the start and stop if nothing specified
        if start_year is None:
            self.start_year = self.bounds[0][0].year
        if stop_year is None:
            self.stop_year = self.bounds[1][0].year
        self.drifts = self.compile_drifts()
