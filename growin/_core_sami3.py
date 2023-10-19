"""functions and variables compute the rtgr using the sami3 model
"""
import os
import datetime
import numpy as np
import warnings
import xarray as xr
import growth_rate as gr
import utils

def merge_sami3_files(date, lon):
    """
    Notes
    Sami3_f4.nc file contains the magnetic coords for the grid.
    Sami3_f1*.nc contains electron density, O+ motion along the field line, ExB drifts
    Sami3_f2*.nc contains neutral densities and neutral wind in the phi-direction
    Sami3_f3*.nc contains electron temperature and additional ion densities.

    the apex longitude is near the midpoint of the fluxtube every time.
    So we use that as a proxy for determining which flux tube is at the
    desired glon
    """
    year = date.year
    day = date.timetuple().tm_yday
    date_str = '{:d}{:03d}'.format(year, day)

    os.chdir('/Users/jmsmit37/data/sami3')
    sami4 = xr.load_dataset('sami3_f4_2014.nc', decode_times=False)
    sub_dat = sami4.isel(nz=80)
    ind = abs(sub_dat.glon[:, 0] - lon).argmin()
    sami4 = sami4.isel(nl=ind)

    sami1 = xr.load_dataset(''.join(['sami3_f1_', date_str, '.nc']), decode_times=False)
    sami1 = sami1.isel(nl=ind)
    sami2 = xr.load_dataset(''.join(['sami3_f2_', date_str, '.nc']), decode_times=False)
    sami2 = sami2.isel(nl=ind)
    sami3 = xr.load_dataset(''.join(['sami3_f3_', date_str, '.nc']), decode_times=False)
    sami3 = sami3.isel(nl=ind)
    sami_out = xr.merge([sami1, sami2, sami3, sami4])

    # xarray needs to be transposed so all the growth calculations don't change
    sami_out = sami_out.transpose()

    sami_out = sami_out.rename({'hrut':'ut', 'vnphi':'u4'})
    sami_out.attrs['lon0'] = lon
    sami_out.to_netcdf(''.join(['sami3_merged_', date_str, '_', str(lon), '.nc']))

def get_growth(sami_filename):
    '''get the sami instrument with growth rates calculated
       checks if there is an existing sami instrument with the appropriate tag
       and loads it. Otherwise it runs the growth rate calculation.

    Parameters
    ----------
    sami_filename : (str)
        full path and name of sami3 file
    lon : (int)
        geo longitude in degrees for SAMI run
    '''

    if os.path.isfile(sami_filename):
        sami = xr.load_dataset(sami_filename, decode_times=False)
    else:
        raise NameError("Invalid Filname, the file:" + sami_filename +
                        " does not exist")

    sami_out = gr.run_growth_calc(sami)
    sami_out.to_netcdf(''.join([sami_filename[:-3], 'growth.nc']))

def get_growth_wedge(sami_filename):
    """Get growth rates for a wedge of longitudes"""
    if os.path.isfile(sami_filename):
        sami = xr.load_dataset(sami_filename, decode_times=False)
        sami = sami.rename({'hrut':'ut', 'vnphi':'u4'})
        sami = sami.transpose()
    else:
        raise NameError("Invalid Filname, the file:" + sami_filename +
                        " does not exist")
    slice_list = []
    for lon in sami.nl:
        slice_growth = gr.run_growth_calc(sami.sel(nl=lon))
        slice_list.append(slice_growth)
        break
    grow_out = xr.concat(slice_list, dim='lon')
    grow_out.to_netcdf(''.join([sami_filename[:-3], 'growth.nc']))
    return grow_out
