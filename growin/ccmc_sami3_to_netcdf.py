import datetime as dt
import numpy as np
import sys
import xarray as xr

run_name = sys.argv[1]
if run_name[0:3] != 'Jon':
    print(' '.join(['invalid run name:', run_name]))
sami3path = ''.join(['/data/sami3/', run_name, '/'])

sz = [304, 124, 96, 25]
mod_vars = {'deneu.dat': 'electron density cm^-3',
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


def sami3data(sami3file, sz):
    fid = open(sami3file, 'rb')
    data = fid.read()
    tempdata = np.frombuffer(data, np.float32)
    fid.close()
    if len(sz) == 4:
        tempdata1 = np.reshape(tempdata,
                               (sz[0] * sz[1] * sz[2] + 2, sz[3]), order='F')
        xx = tempdata1[1:-1, :]
    else:
        xx = tempdata[1:-1]
    samidata = np.reshape(xx, sz, order='F')
    return samidata


def sami3data_grid(sami3file, sz):
    fid = open(sami3file, 'rb')
    data = fid.read()
    tempdata = np.frombuffer(data, np.float32)
    fid.close()
    xx = tempdata[1:-1]
    samidata = np.reshape(xx, sz, order='F')
    return samidata


time = np.loadtxt(sami3path + 'time.dat')
ut = time[:, 1] + time[:, 2] / 60 + time[:, 3] / 3600
glat = '{0}{1}'.format(sami3path, 'glatu.dat')
glon = '{0}{1}'.format(sami3path, 'glonu.dat')
zalt = '{0}{1}'.format(sami3path, 'zaltu.dat')
lat_coord = sami3data(glat, sz[0:3])
lon_coord = sami3data(glon, sz[0:3])
zalt_coord = sami3data(zalt, sz[0:3])

sami_out = xr.Dataset(coords=dict(ut=(['nt'], ut),
                                  glat=(['nz', 'nf', 'nlt'], lat_coord),
                                  glon=(['nz', 'nf', 'nlt'], lon_coord),
                                  zalt=(['nz', 'nf', 'nlt'], zalt_coord)))

for var_file in mod_vars:
    print(var_file)
    buff = '{0}{1}'.format(sami3path, var_file)
    tmp_var = sami3data(buff, sz)
    # For reference, for some data there is reduced dimension
    # TEC = sami3data_grid(sami3_tec,sz[1:4])
    sami_out[var_file[:-5]] = (['nz', 'nf', 'nlt', 'nt'], tmp_var,
                               {'desc': mod_vars[var_file]})

sami_out = sami_out.rename({'u1': 'u4'})

apex_ind = sami_out.zalt[:, 0, 0].argmax()
lon_ind = abs(sami_out.glon[apex_ind, 0, :] - 284).argmin()
sami_out = sami_out.isel(nlt=lon_ind)

try:
    with open('{0}{1}'.format(sami3path, 'SAMI3_list')) as f:
        lines = f.readlines()
        date = date = lines[1][11:21]
        day = dt.datetime.strptime(date, '%Y/%m/%d').timetuple().tm_yday
        year = dt.datetime.strptime(date, '%Y/%m/%d').year
except FileNotFoundError:
    with open('{0}{1}'.format(sami3path, 'namelist_input.dat')) as f:
        lines = f.readlines()
        year = int(lines[0])
        date = ','.join([lines[0][:-1], lines[1][:-1]])
        day = dt.datetime.strptime(date, '%Y,%m,%d').timetuple().tm_yday

date_str = '{:d}{:03d}'.format(year, day)
sami_out['day'] = day
sami_out['year'] = year
lon = int(sami_out.glon.mean().values)
sami_out.attrs['lon0'] = lon
sami_out.to_netcdf(''.join([sami3path, 'sami3_merged_', date_str,
                            '_', str(lon), '.nc']))
