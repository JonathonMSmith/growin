"""functions and variables to run the requisite models and functions to
calculate the drifts and and from them run the model and then compute the rtgr
year-1,11,1     Dec Sol
year,1,31

year,2,1        Mar Eq
year,4,30

year,5,1        Jun Sol
year,7,31

year,8,1        Sep Eq
year,10,31
"""
import os
import pickle
import numpy as np
import xarray as xr
import sami2py
import irfl

# default values to organize seasons and zones
# months that bound each season
SEASON_BOUNDS = np.array([11, 2, 5, 8, 11])
SEASON_NAMES = ['DecSol', 'MarEq', 'JunSol', 'SepEq']
# midpoint days for each season for sami
SEASON_DAYS = {'DecSol':355, 'MarEq':80, 'JunSol':155, 'SepEq':266}

# longitude regions shifted by +15 degrees for ease of binning
LONGITUDE_BOUNDS = np.array([0, 75, 145, 300, 360])
LON_SECTOR_NAMES = ['African', 'Indian', 'Pacific', 'South American']
# midpoints from longitude_bounds, subtract 15 again for sami
MODEL_SECTORS = {'African':22, 'Indian':95, 'Pacific':207,
                 'South American':315}

# custom functions for pysat instrument to modify the data for use here
def shift_local_time(inst):
    """shift local times so that they are centered on midnight"""
    idx, = np.where(inst['slt'] < 12.)
    inst[idx, 'slt'] += 24.

def shift_longitude(inst):
    """shift longitude values 15 degrees for longitude sector binning"""
    inst['glon'] += 15.

def drift_fix(inst):
    """in the cnofs data positive drifts are toward earth, so sign change"""
    inst['ionVelocityZ'] *= -1

def filter_ivm(inst):
    """some simple critereon to filter the data"""
    idx, = np.where((inst['Ni'] >= 3000) & (inst['apex_altitude'] <= 550))
    inst.data = inst.data.iloc[idx]

def get_drifts(start=2008, stop=2014, clean_level='none',
               drift_key='ionVelocityZ'):
    """create/load the instrument and obtain the drifts then save the drifts"""
    path = irfl.utils.generate_path('drift', year=start, end_year=stop)
    drift_f_name = os.path.join(path, clean_level+'_'+drift_key+'.p')
    if os.path.isfile(drift_f_name):
        drift_inst = pickle.load(open(drift_f_name, 'rb'))
    else:
        drift_inst = irfl.DriftInstrument(platform='cnofs', name='ivm',
                                          clean_level=clean_level)
        drift_inst.custom.add(drift_fix, 'modify')
        drift_inst.custom.add(shift_longitude, 'modify')
#        drift_inst.custom.add(filter_ivm, 'modify')
        drift_inst.exb_fourier_fit(drift_key=drift_key,
                                   lon_bins=LONGITUDE_BOUNDS,
                                   season_bins=SEASON_BOUNDS,
                                   season_names=SEASON_NAMES,
                                   zone_labels=LON_SECTOR_NAMES,
                                   start_year=start, stop_year=stop)
        if not os.path.isdir(path):
            os.makedirs(path)
        pickle.dump(drift_inst, open(drift_f_name, 'wb'))
    return drift_inst

def get_growth(tag, day, year, lon, exb_drifts, ve01=0):

    path = irfl.utils.generate_path('growth', year=year,
                                    lon=lon, day=day)
    sami_filename = os.path.join(path, 'sami'+tag+'.p')

    if os.path.isfile(sami_filename):
        sami = pickle.load(open(sami_filename, 'rb'))
        return sami
    sami2py.run_model(tag=tag, day=day, year=year, lon=lon, fejer=False,
                      ExB_drifts=exb_drifts, ve01=0, outn=True)
    sami = sami2py.Model(tag=tag, day=day,
                         year=year, lon=lon, outn=True)
    sami.gamma = irfl.growth_rate.run_growth_calc(sami,
                                                  exb_drifts)
    if not os.path.isdir(path):
        os.makedirs(path)
    pickle.dump(sami, open(sami_filename, 'wb'))
    return sami



# run the model for each year and season and compute the growth rate and plot
def get_growth_rates_survey(start=2008, stop=2014, clean_level='none',
                            drift_key='ionVelocityZ'):
    """calculate the growth rate from the sami model using computed drifts"""
    drift_inst = get_drifts(start=start, stop=stop, clean_level=clean_level,
                            drift_key=drift_key)

    drift = drift_inst.drifts
    for year in range(start, stop):
        for season in SEASON_NAMES:
            day = SEASON_DAYS[season]
            for zone in LON_SECTOR_NAMES:
                lon = MODEL_SECTORS[zone]
                exb_drifts = drift.coefficients.sel(year=year, season=season,
                                                    longitude=zone).values
                ve01 = drift.ve01.sel(year=year, season=season,
                                      longitude=zone).values
                sami = get_growth(tag=clean_level, day=day, year=year, lon=lon,
                                  exb_drifts=exb_drifts)
                irfl.generate_plots.plot_drifts(drift_inst, year=year,
                                                season=season, zone=zone,
                                                day=day, lon=lon)
                irfl.generate_plots.plot_growth_term(sami, season, drift=True)
                irfl.generate_plots.plot_growth_term(sami, season, 'N',
                                                     vmax=None)
                irfl.generate_plots.plot_growth_term(sami, season, 'K',
                                                     vmax=8*10**(-5))
                irfl.generate_plots.plot_growth_term(sami, season, 'g_nu',
                                                     vmax=20)
