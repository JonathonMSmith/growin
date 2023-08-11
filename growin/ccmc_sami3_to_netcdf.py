import datetime as dt
import numpy as np
import xarray as xr

sami3path = '../models/ccmc_run/'
sz = [304, 124, 96, 193]
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

def sami3data(sami3file,sz):
    fid= open(sami3file,'rb')
    data = fid.read()
    tempdata = np.frombuffer(data, np.float32)
    fid.close()
    if len(sz)==4:
       tempdata1 = np.reshape(tempdata,(sz[0]*sz[1]*sz[2]+2,sz[3]),order='F')
       xx = tempdata1[1:-1,:]
    else:
       xx = tempdata[1:-1]
    samidata= np.reshape(xx,sz,order='F')
    return samidata

def sami3data_grid(sami3file,sz):
    fid= open(sami3file,'rb')
    data = fid.read()
    tempdata = np.frombuffer(data, np.float32)
    fid.close()
    xx = tempdata[1:-1]
    samidata= np.reshape(xx,sz,order='F')
    return samidata

time = np.loadtxt(sami3path+'time.dat')
ut = time[:,1] + time[:,2] / 60 + time[:,3] / 3600
glat = '{0}{1}'.format(sami3path, 'glatu.dat')
glon = '{0}{1}'.format(sami3path, 'glonu.dat')
zalt = '{0}{1}'.format(sami3path, 'zaltu.dat')
lat_coord = sami3data(glat, sz[0:3])
lon_coord = sami3data(glon, sz[0:3])
zalt_coord = sami3data(zalt, sz[0:3])

sami_dat = xr.Dataset(coords = dict(ut = (['nt'], ut), 
                                    glat=(['nz', 'nf', 'nlt'], lat_coord),
                                    glon=(['nz', 'nf', 'nlt'], lon_coord),
                                    zalt=(['nz', 'nf', 'nlt'], zalt_coord)))

for var_file in mod_vars:
    print(var_file)
    buff = '{0}{1}'.format(sami3path, var_file)
    # why is this throwing a damned value error dude (everything with a u in th efilename?
    tmp_var = sami3data(buff, sz)
    # For reference, for some data there is reduced dimension
    # TEC = sami3data_grid(sami3_tec,sz[1:4])
    sami_dat[var_file[:-5]] = (['nz', 'nf', 'nlt', 'nt'], tmp_var,
                               {'desc': mod_vars[var_file]})

sami_dat = sami_dat.rename({'u1':'u4'})

apex_ind = sami_dat.zalt[:, 0, 0].argmax()
lon_ind = abs(sami_dat.glon[apex_ind, 0, :] - 284).argmin()
sami_dat = sami_dat.isel(nlt=lon_ind)
with open('{0}{1}'.format(sami3path, 'SAMI3_list')) as f:
    lines  = f.readlines()
    date = date = lines[1][11:21]
    day = dt.datetime.strptime(date, '%Y/%m/%d').timetuple().tm_yday
    year = dt.datetime.strptime(date, '%Y/%m/%d').year
sami_dat['day'] = day
sami_dat['year'] = year
sami_dat.attrs['lon0'] = sami_dat.glon.mean().values
sami_dat.to_netcdf('ccmc_sami.nc')
