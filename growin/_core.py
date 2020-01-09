"""functions and variables to run the requisite models and functions to
calculate the drifts and and from them run the model and then compute the rtgr
"""
import os
import pickle
import datetime
import numpy as np
import sami2py
import growin

# custom functions for pysat instrument to modify the data for use here
def shift_local_time(inst):
    """shift local times so that they are centered on midnight"""
    idx, = np.where(inst['slt'] < 12.)
    inst[idx, 'slt'] += 24.


def shift_longitude(inst, offset=None):
    """shift longitude values some degrees for longitude sector binning
        there has got to be a better way"""
    if offset is None:
        offset = 15
    inst['glon'] += offset
    idx, = np.where(inst['glon'] > 360)
    inst[idx, 'glon'] -= 360
    idx_neg, = np.where(inst['glon'] < 0)
    inst[idx_neg, 'glon'] += 360


def drift_fix(inst):
    """in the cnofs data positive drifts are toward earth, so sign change
       not needed for meridional drifts"""
    inst['ionVelocityZ'] *= -1


def get_drifts(start=2008, stop=2014, clean_level='none', drift_inst=None,
               drift_key='ionVelocityZ', season_names: list = None,
               season_bounds: list = None, zone_names: list = None,
               zone_bounds: list = None, offset: int = None):

    """create/load the instrument and obtain the drifts then save the drifts

    Parameters
    ----------
    start : (int)
         start year of the survey
    stop : (int)
         stop year of the survey
    clean_level : (string)
         specify cleaning routine for pysat
    drift_inst: (growin.fourier_exb.DriftInstrument)
         drift instrument if different custom modifier functions are desired
    drift_key : (int)
         dictionary key for the pysat instrument drift values to use
    season_names : (array-like of strings)
         array-like containing the names of the specified seasons
    season_bounds : (array-like of int or float)
         array-like ocntaining the days that delineate the season bounds
    zone_names : (array-like of strings)
         array-like of stings specifying the names of longitude zones used
    zone_bounds : (array-like of int or float)
         array-like of longitudes in degrees that delineate the zone bounds
    """
    path = growin.utils.generate_path('drift', year=start, end_year=stop)
    #change this to be more specific... maybe add a tag or something
    drift_f_name = os.path.join(path, clean_level+'_'+drift_key+'.p')

    if os.path.isfile(drift_f_name):
        drift_inst = pickle.load(open(drift_f_name, 'rb'))
        return drift_inst

    if isinstance(drift_inst, growin.fourier_exb.DriftInstrument):
        drift_inst = drift_inst
    else:
        drift_inst = growin.DriftInstrument(platform='cnofs', name='ivm',
                                            clean_level=clean_level)
        drift_inst.custom.add(drift_fix, 'modify')
        drift_inst.custom.add(shift_longitude, 'modify', offset)
    drift_inst.get_drifts(drift_key=drift_key,
                          lon_bins=zone_bounds,
                          season_bins=season_bounds,
                          season_names=season_names,
                          zone_labels=zone_names,
                          start_year=start, stop_year=stop)
    if not os.path.isdir(path):
        os.makedirs(path)
    pickle.dump(drift_inst, open(drift_f_name, 'wb'))
    return drift_inst


def get_growth(tag, day, year, lon, exb_drifts, ve01=0, f10=120.0):
    '''get the sami instrument with growth rates calculated
       checks if there is an existing sami instrument with the appropriate tag
       and loads it. Otherwise it runs the growth rate calculation.

    Parameters
    ----------
    tag : (string)
        name of run where growth is/will be archived
    day : (int)
        day of year for SAMI run
    year : (int)
        year for SAMI run
    lon : (int)
        geo longitude in degrees for SAMI run
    exb_drifts : (10x2 ndarray of floats)
        Matrix of Fourier series coefficients dependent on solar local time
        (SLT) in hours where
        exb_total = exv_drifts[i,0]*cos((i+1)*pi*SLT/12)
                  + exb_drifts[i,1]*sin((i+1)*pi*SLT/12)
    ve01 : (float)
        offset for Fourier exb drifts, not used by default therefore
        we are assuming net zero vertical drift
    '''
    path = growin.utils.generate_path('growth', year=year,
                                      lon=lon, day=day)
    sami_filename = os.path.join(path, 'sami'+tag+'.p')

    if os.path.isfile(sami_filename):
        sami = pickle.load(open(sami_filename, 'rb'))
        return sami
    if exb_drifts is not None:
        sami2py.run_model(tag=tag, day=day, year=year, lon=lon, fejer=False,
                          ExB_drifts=exb_drifts, ve01=0, outn=True,
                          f107=f10, f107a=f10)
    else:
        sami2py.run_model(tag=tag, day=day, year=year, lon=lon, fejer=True,
                          outn=True, f107=f10, f107a=f10)

    sami = sami2py.Model(tag=tag, day=day,
                         year=year, lon=lon, outn=True)
    sami.gamma = growin.growth_rate.run_growth_calc(sami, exb_drifts)

    if not os.path.isdir(path):
        os.makedirs(path)
    pickle.dump(sami, open(sami_filename, 'wb'))
    return sami


def fit_fejer(year, day, lon):
    '''Compute the fourier coefficients for the Fejer-Scherliess model

    Parameters
    ----------
    year : (int)
        year to use for Fejer-Scherliess drifts
    day : int
        julian day
    lon : int or double
        longitude in degrees
    '''
    import pyglow
    delta_t = datetime.timedelta(day-1)
    slt_step = np.linspace(0, 23.5, 48)
    # convert from LT to UT because IRI uses UT exclusively
    ut_step = slt_step - lon/15
    ut_step = [t+24 if t < 0 else t for t in ut_step]
    drifts = []
    # extract the Fejer-Scherliess drifts from IRI via pyglow
    for t in ut_step:
        hour = int(t)
        minute = int(60 * ((t) % 1))
        day = datetime.datetime(year, 1, 1, hour, minute) + delta_t
        point = pyglow.Point(day, 0, lon, 250)
        point.run_iri()
        drifts.append(point.exb)
    drifts = np.array(drifts)
    # compute the coefficients
    ve01, exb_drifts = growin.fourier_exb.fourier_fit(slt_step, drifts, 10)
    return exb_drifts


def get_growth_rates_survey(start: int, stop: int, clean_level: str,
                            drift_key: str, season_names: list,
                            season_bounds: list, season_days: dict,
                            zone_names: list, zone_bounds: list,
                            zone_lons: dict):
    """calculate the growth rate from the sami model using computed drifts
       run the model for each year and season
       compute the growth rate and plot

    Parameters
    ----------
    start : (int)
         start year of the survey
    stop : (int)
         stop year of the survey
    clean_level : (string)
         specify cleaning routine for pysat
    drift_key : (int)
         dictionary key for the pysat instrument drift values to use
         a good default for cnofs is 'IonVelmeridional'
    season_names : (array-like of strings)
         array-like containing the names of the specified seasons
    season_bounds : (array-like of int or float)
         array-like ocntaining the days that delineate the season bounds
    season_days : (dict)
         dictionary with season names as keys, and as values the day to be
         used by SAMI
    zone_names : (array-like of strings)
         array-like of stings specifying the names of longitude zones used
    zone_bounds : (array-like of int or float)
         array-like of longitudes in degrees that delineate the zone bounds
    zone_lons : (dict)
         dictionary with zone names as keys, and as values the longitude to
         be used by SAMI
    """

    if drift_key != 'Fejer':
        drift_inst = get_drifts(start=start, stop=stop,
                                clean_level=clean_level, drift_key=drift_key,
                                season_names=season_names,
                                season_bounds=season_bounds,
                                zone_bounds=zone_bounds,
                                zone_names=zone_names)
        drift = drift_inst.drifts
    else:
        drift = None

    for year in range(start, stop):
        for season in season_names:
            day = season_days[season]
            for zone in zone_names:
                lon = zone_lons[zone]
                if drift:
                    exb_drifts = drift.coefficients.sel(year=year,
                                                        season=season,
                                                        longitude=zone).values
                else:
                    exb_drifts = fit_fejer(year, day, lon)

                tag = clean_level + '_' + drift_key
                sami = get_growth(tag=tag, day=day, year=year, lon=lon,
                                  exb_drifts=exb_drifts)
    return sami
