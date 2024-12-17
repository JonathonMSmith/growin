import datetime as dt
import json
import numpy as np
import os
import pandas as pd
import psutil
import sys
import urllib.request
import xarray as xr

storm = 'march2023'
if storm == 'feb2022':
    RUN_NAME = 'SAMI3-TIEGCM-01_2022-02-TP-01_082823_IT_1'
    NUM_DAYS = 3
    start_date = dt.datetime(2022, 2, 2)
elif storm == 'march2023':
    RUN_NAME = 'SAMI3-TIEGCM-01_2023-03-TP-01_081823_IT_1'
    NUM_DAYS = 6
    start_date = dt.datetime(2023, 3, 20)
elif storm == 'april2023':
    RUN_NAME = 'SAMI3-TIEGCM-01_2023-04-TP-01_011724_IT_1'
    NUM_DAYS = 5
    start_date = dt.datetime(2023, 4, 21)
# RUN_NAME = sys.argv[1]
#if RUN_NAME[0:3] != 'Jon':
#    print(' '.join(['invalid run name:', RUN_NAME]))


laptop = True
if laptop:
    prefix = '/Users/jklenzin'
else:
    prefix = '/Volumes/Expansion'
SAMI3PATH = ''.join([prefix, '/data/sami3/', RUN_NAME, '/'])

# NUM_DAYS=6
# this is the size for the iconTIEGCM runs sz = [304, 124, 96, 25]
SZ = [304, 124, 96, 1+96*NUM_DAYS]
RG_SZ = [100,100,96,1+96*NUM_DAYS]
MOD_VARS = {'deneu.dat': 'electron density cm^-3',
            'u1pu.dat': 'meridional ExB velocity cm/s',
            'denn1u.dat': 'H density cm^-3',
            'denn2u.dat': 'O density cm^-3',
            'denn3u.dat': 'NO density cm^-3',
            'denn4u.dat': 'O2 density cm^-3',
            'denn5u.dat': 'He density cm^-3',
            'denn6u.dat': 'N2 density cm^-3',
            'denn7u.dat': 'N density cm^-3',
            'u1u.dat': 'zonal neutral wind velocity cm/s',
            'teu.dat': 'electron temperature K',
            'deni1u.dat': 'H+ density cm^-3',
            'deni2u.dat': 'O+ density cm^-3',
            'deni3u.dat': 'NO+ density cm^-3',
            'deni4u.dat': 'O2+ density cm^-3',
            'deni5u.dat': 'He+ density cm^-3',
            'deni6u.dat': 'N2+ density cm^-3',
            'deni7u.dat': 'N+ density cm^-3'}

MET_VARS = {'time.dat': 'time seconds',
            'glatu.dat': 'geo latitude degrees',
            'glonu.dat': 'geo longitude degrees',
            'blatu.dat': 'geo latitude degrees',
            'blonu.dat': 'geo longitude degrees',
            'zaltu.dat': 'altitude km',
            'SAMI3_List': 'namelist variant',
            'sami3-3.22.namelist': 'namelist for fortran code',
            'namelist_input.dat': 'namelist variant'}

REG_VARS = {'glon0.dat': 'geo longitude degrees',
            'glat0.dat': 'geo latitude degrees',
            'hmf2u.dat': 'height of the f layer peak km',
            'tecu.dat': 'total electron content',
            'zalt0.dat': 'altitude km'}

def download_run(data_dir, run_name, mod_vars, dim_vars):
    """download ror_file from ccmc runs on request output"""
    run_url = ''.join(['https://ccmc.gsfc.nasa.gov/results/',
                       'run_files.php?runnumber=', run_name])
    json_url = ''.join([run_url, '&format=json'])
    if not os.path.exists(data_dir):
        os.mkdir(data_dir)
    # TODO use a session object with requests instead of urllib
    # get json
    urllib.request.urlretrieve(json_url, run_name)
    # iterate through files and download
    with open(run_name) as run:
        son = json.load(run)
        file_list = son['files']
        for file in file_list:
            full_name = os.path.join(data_dir, file['file'])
            if os.path.exists(full_name):
                continue
            print(full_name)
            if file['file'] in mod_vars.keys() or file['file'] in dim_vars.keys():
                print('downloading')
                urllib.request.urlretrieve(file['link'], full_name)

def sami3data(sami3file, sz):
    "reshape and format unformatted binary .dat file"
    fid = open(sami3file, 'rb')
    data = fid.read()
    tempdata = np.frombuffer(data, np.float32)
    fid.close()
    if len(sz) == 4:
        tempdata = np.reshape(tempdata,
                               (sz[0] * sz[1] * sz[2] + 2, sz[3]), order='F')
        xx = tempdata[1:-1, :]
    else:
        xx = tempdata[1:-1]
    samidata = np.reshape(xx, sz, order='F')
    process = psutil.Process(os.getpid())
    print(f"inside the function, it is using {process.memory_info().rss / 1024 / 1024} MB")  # in Megabytes
    return samidata


def sami3data_grid(sami3file, sz):
    "reshape and format regridded binary .dat file"
    fid = open(sami3file, 'rb')
    data = fid.read()
    tempdata = np.frombuffer(data, np.float32)
    fid.close()
    xx = tempdata[1:-1]
    samidata = np.reshape(xx, sz, order='F')
    return samidata

def combine_regridded_netcdf(sami3path, sz, reg_vars, zone = (270, 310)):
    """
    Parameters
    ----------
    sami3path: str
        path to sami3 data output
    sz: list
        dimensions of sami3 data to reshape .dat files
    reg_vars: dict
        dictionary containing regridded variables of interest and their units
    zone: tuple
        longitude bounds for downslection
    """

    time = np.loadtxt(sami3path + 'time.dat')
    ut = time[:, 1] + time[:, 2] / 60 + time[:, 3] / 3600
    glat = os.path.join(sami3path, 'glat0.dat')
    glon = os.path.join(sami3path, 'glon0.dat')
    zalt = os.path.join(sami3path, 'zalt0.dat')
    hmf2 = os.path.join(sami3path, 'hmf2u.dat')
    tec = os.path.join(sami3path, 'tecu.dat')
    lat_coord = sami3data_grid(glat, sz[0:3])
    lon_coord = sami3data_grid(glon, sz[0:3])
    zalt_coord = sami3data_grid(zalt, sz[0:3])
    hmf2_var = sami3data_grid(hmf2, sz[1:4])
    tec_var = sami3data_grid(tec, sz[1:4])
    sami_out = xr.Dataset(coords=dict(ut=(['nt'], ut),
                                      glat=(['nx', 'ny', 'nl'], lat_coord),
                                      glon=(['nx', 'ny', 'nl'], lon_coord),
                                      zalt=(['nx', 'ny', 'nl'], zalt_coord)),
                          data_vars=dict(hmf2=(['nx', 'nl', 'nt'], hmf2_var),
                                         tec=(['nx', 'nl', 'nt'], tec_var)))
    apex_ind = sami_out.zalt[:, 0, 0].argmax()
    lon_ind, = np.where((sami_out.glon[apex_ind, 0, :] > zone[0]) &
                       (sami_out.glon[apex_ind, 0, :] < zone[1]))
    sami_out = sami_out.isel(nl=lon_ind)
    lon = int(sami_out.glon.mean().values)
    sami_out.to_netcdf(''.join([sami3path, 'sami3_reg_merged_',
                                '_', str(lon), '.nc']))


def combine_global_regridded_netcdf(sami3path, sz, reg_vars):
    """
    Parameters
    ----------
    sami3path: str
        path to sami3 data output
    sz: list
        dimensions of sami3 data to reshape .dat files
    reg_vars: dict
        dictionary containing regridded variables of interest and their units
    zone: tuple
        longitude bounds for downslection
    """

    time = np.loadtxt(sami3path + 'time.dat')
    ut = time[:, 1] + time[:, 2] / 60 + time[:, 3] / 3600
    delta_day = np.floor(time[:,4]/24)
    num_vals = len(delta_day)
    time = pd.to_datetime({'year': np.ones(num_vals) * start_date.year,
                           'month': np.ones(num_vals) * start_date.month,
                           'day': np.ones(num_vals) * start_date.day + delta_day,
                           'hour': time[:, 1],
                           'minute': time[:, 2],
                           'second': time[:, 3]})

    glat = os.path.join(sami3path, 'glat0.dat')
    glon = os.path.join(sami3path, 'glon0.dat')
    zalt = os.path.join(sami3path, 'zalt0.dat')

    hmf2 = os.path.join(sami3path, 'hmf2u.dat')
    tec = os.path.join(sami3path, 'tecu.dat')

    # Reshape and drop superfluous dimension
    lat_coord = sami3data_grid(glat, sz[0:3])[:, 0, :]
    lon_coord = sami3data_grid(glon, sz[0:3])[:, 0, :]
    zalt_coord = sami3data_grid(zalt, sz[0:3])[:, 0, :]
    hmf2_var = sami3data_grid(hmf2, sz[1:4])
    tec_var = sami3data_grid(tec, sz[1:4])

    sami_out = xr.Dataset(coords=dict(time=(['nt'], time),
                                      glat=(['nx', 'nl'], lat_coord),
                                      glon=(['nx', 'nl'], lon_coord),
                                      zalt=(['nx', 'nl'], zalt_coord)),
                          data_vars=dict(hmf2=(['nx', 'nl', 'nt'], hmf2_var),
                                         tec=(['nx', 'nl', 'nt'], tec_var)))
    sami_out.to_netcdf(os.path.join(sami3path, 'sami3_reg_merged.nc'))

    return sami_out



def combine_in_netcdf(sami3path, sz, mod_vars, zone= (270, 310)):
    """
    Parameters
    ----------
    sami3path: str
        path to sami3 data output
    sz: list
        dimensions of sami3 data to reshape .dat files
    mod_vars: dict
       dictionary containing variables of interest and their units
    zone: tuple
        longitude bounds for downslection
    """
    time = np.loadtxt(sami3path + 'time.dat')
    ut = time[:, 1] + time[:, 2] / 60 + time[:, 3] / 3600
    glat = os.path.join(sami3path, 'glatu.dat')
    glon = os.path.join(sami3path, 'glonu.dat')
    mlat = os.path.join(sami3path, 'blatu.dat')
    mlon = os.path.join(sami3path, 'blonu.dat')
    zalt = os.path.join(sami3path, 'zaltu.dat')
    lat_coord = sami3data(glat, sz[0:3])
    lon_coord = sami3data(glon, sz[0:3])
    mlat_coord = sami3data(mlat, sz[0:3])
    mlon_coord = sami3data(mlon, sz[0:3])
    zalt_coord = sami3data(zalt, sz[0:3])
    sami_out = xr.Dataset(coords=dict(ut=(['nt'], ut),
                                      glat=(['nz', 'nf', 'nlt'], lat_coord),
                                      glon=(['nz', 'nf', 'nlt'], lon_coord),
                                      mlat=(['nz', 'nf', 'nlt'], mlat_coord),
                                      mlon=(['nz', 'nf', 'nlt'], mlon_coord),
                                      zalt=(['nz', 'nf', 'nlt'], zalt_coord)))
    apex_ind = sami_out.zalt[:, 0, 0].argmax()
    lon_ind, = np.where((sami_out.glon[apex_ind, 0, :] > zone[0]) &
                       (sami_out.glon[apex_ind, 0, :] < zone[1]))
    print(sami_out.glon[apex_ind, 0, lon_ind])
    sami_out = sami_out.isel(nlt=lon_ind)
    nc_flist = []

    for var_file in mod_vars:
        buff = os.path.join(sami3path, var_file)
        varname = var_file[:-5]
        print(varname)
        var_fname = ''.join([sami3path, 'trimmed_', varname, '.nc'])
        nc_flist.append(var_fname)
        if os.path.isfile(var_fname):
            continue
        tmp_var = sami3data(buff, sz)
        process = psutil.Process(os.getpid())
        print(f"outside the function, it is using {process.memory_info().rss / 1024 / 1024} MB")
        var_arr = xr.DataArray(dims=['nz', 'nf', 'nlt', 'nt'], data=tmp_var,
                               attrs={'desc': mod_vars[var_file]})
        print(lon_ind)

        var_arr = var_arr.isel(nlt=lon_ind)
        print(var_arr)
        var_arr.to_netcdf(var_fname)
       # For reference, for some data there is reduced dimension
       # TEC = sami3data_grid(sami3_tec,sz[1:4])
    print(lon_ind)
    for f in nc_flist:
        tmp_var = xr.load_dataset(f)
        tmp_var = tmp_var['__xarray_dataarray_variable__']
        fname = os.path.split(f)
        print(fname[-1][8:-3])
        sami_out[fname[-1][8:-3]] = tmp_var

    sami_out = sami_out.rename({'u1': 'u4'})

    try:
        with open(os.path.join(sami3path, 'SAMI3_list')) as f:
            lines = f.readlines()
            date = date = lines[1][11:21]
            day = dt.datetime.strptime(date, '%Y/%m/%d').timetuple().tm_yday
            year = dt.datetime.strptime(date, '%Y/%m/%d').year
    except FileNotFoundError:
        try:
            with open(os.path.join(sami3path, 'namelist_input.dat')) as f:
                lines = f.readlines()
                year = int(lines[0])
                date = ','.join([lines[0][:-1], lines[1][:-1]])
                day = dt.datetime.strptime(date, '%Y,%m,%d').timetuple().tm_yday
        except FileNotFoundError:
            with open(os.path.join(sami3path, 'sami3-3.22.namelist')) as f:
                lines = f.readlines()
                year = int(lines[11][10:14])
                day = int(lines[12][8:].split(',')[0])
    date_str = '{:d}{:03d}'.format(year, day)
    sami_out['day'] = day
    sami_out['year'] = year
    lon = int(sami_out.glon.mean().values)
    sami_out.attrs['lon0'] = lon
    sami_out.to_netcdf(''.join([sami3path, 'sami3_merged_', date_str,
                                '_', str(lon), '.nc']))

# download_run(SAMI3PATH, RUN_NAME, MOD_VARS, MET_VARS)
download_run(SAMI3PATH, RUN_NAME, REG_VARS, MET_VARS)
# combine_in_netcdf(SAMI3PATH, SZ, MOD_VARS)
# combine_regridded_netcdf(SAMI3PATH, RG_SZ, REG_VARS)
sami_out = combine_global_regridded_netcdf(SAMI3PATH, RG_SZ, {})
# Sort out arrays by glon
sami_out = sami_out.sortby(sami_out['glon'][0,:])