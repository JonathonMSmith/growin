"""A place for python script to run the growin calculation under various
   conditions to explore the sensitivity of the results
   
"""

import pickle
import numpy as np
import growin
import f10_tool
import fake_drifts
import generate_plots


# default values to organize seasons and zones
# months that bound each season
SEASON_BOUNDS = np.array([11, 2, 5, 8, 11])
SEASON_NAMES = ['December', 'March', 'June', 'September']
SEASONS = {'December': 355, 'March': 80, 'June': 155, 'September': 266}
# longitude regions shifted by +15 degrees for ease of binning
ZONE_BOUNDS = np.array([0, 75, 145, 300, 360])
ZONE_NAMES = ['African', 'Asian', 'Pacific', 'American']
ZONES = {'African': 22, 'Asian': 95, 'Pacific': 207, 'American': 315}

def sin_growth(year, season, zone, f10):
    '''growth driven by sinusoidal drifts'''
    print('sinusoidal')
    lon = ZONES[zone]
    day = SEASONS[season]
    drift = fake_drifts.sinusoidal_drifts()
    time = np.linspace(0, 24, 49)
    offset, fit_coef = growin.fourier_exb.fourier_fit(time, drift, 10)
    sami = growin.get_growth(tag='sinusoidal', day=day, year=year, lon=lon,
                             exb_drifts=fit_coef, f10=f10)
    generate_plots.plot_growth_term(sami=sami, season=season, drift=True,
                                    peak=False)

def flat_sin_growth(year, season, zone, f10):
    '''growth driven by flat sinusoidal drifts'''
    print('flat')
    lon = ZONES[zone]
    day = SEASONS[season]
    drift = fake_drifts.sinusoidal_drifts(amplitude=1)
    time = np.linspace(0, 24, 49)
    offset, fit_coef = growin.fourier_exb.fourier_fit(time, drift, 10)
    sami = growin.get_growth(tag='sinusoidal_flat', day=day, year=year, lon=lon,
                             exb_drifts=fit_coef, f10=f10)
    generate_plots.plot_growth_term(sami=sami, season=season, drift=True,
                                    peak=False)

def pre_sin_growth(year, season, zone, f10):
    '''growth drien by sinusoidal drifts with pre'''
    print('PRE')
    lon = ZONES[zone]
    day = SEASONS[season]
    drift = fake_drifts.sinusoidal_pre_drifts()
    time = np.linspace(0, 24, 49)
    offset, fit_coef = growin.fourier_exb.fourier_fit(time, drift, 10)
    sami = growin.get_growth(tag='sinusoidal_pre', day=day, year=year, lon=lon,
                             exb_drifts=fit_coef, f10=f10)
    generate_plots.plot_growth_term(sami=sami, season=season, drift=True,
                                    peak=False)

def cnofs_growth(year, season, zone, bubbs, f10):
    '''growth driven by cnofs drifts'''
    #get cnofs drifts
    print('CNOFS')
    lon = ZONES[zone]
    day = SEASONS[season]
    exb = growin._core.get_drifts(start=2008, stop=year, clean_level='dusty',
                                  drift_key='ionVelmeridional',
                                  season_names=SEASON_NAMES,
                                  season_bounds=SEASON_BOUNDS,
                                  zone_names=ZONE_NAMES,
                                  zone_bounds=ZONE_BOUNDS)
    exb_coef = exb.drifts.coefficients.sel(year=year-1, season=season,
                                           longitude=zone).values
    sami = growin.get_growth(tag='cnofs', day=day, year=year, lon=lon,
                             exb_drifts=exb_coef, f10=f10)
    if bubbs:
        sea_idx = SEASON_NAMES.index(season)
        months = [SEASON_BOUNDS[sea_idx], SEASON_BOUNDS[sea_idx+1]]
        lon_idx = ZONE_NAMES.index(zone)
        lons = [ZONE_BOUNDS[lon_idx], ZONE_BOUNDS[lon_idx+1]]
        generate_plots.plot_growth_term(sami=sami, season=season, drift=True,
                                        bubbs=True, peak=False,
                                        months=months, lons=lons)
    else:
        generate_plots.plot_growth_term(sami=sami, season=season, drift=True,
                                        peak=False)


def fejer_growth(year, season, zone, bubbs, f10):
    '''growth driven by fejer drifts'''
    print('FEJER')
    lon = ZONES[zone]
    day = SEASONS[season]
    fejer_coef = growin._core.fit_fejer(year, day, lon)
    sami = growin.get_growth(tag='fejer', day=day, year=year, lon=lon,
                             exb_drifts=fejer_coef, f10=f10)
    if bubbs:
        sea_idx = SEASON_NAMES.index(season)
        months = [SEASON_BOUNDS[sea_idx], SEASON_BOUNDS[sea_idx+1]]
        lon_idx = ZONE_NAMES.index(zone)
        lons = [ZONE_BOUNDS[lon_idx], ZONE_BOUNDS[lon_idx+1]]
        generate_plots.plot_growth_term(sami=sami, season=season, drift=True,
                                        bubbs=True, peak=False,
                                        months=months, lons=lons)
    else:
        generate_plots.plot_growth_term(sami=sami, season=season, drift=True,
                                        peak=False)
def run_survey(year):
    year = year
    for z in ZONES:
        for s in SEASONS:
            day = SEASONS[s]
            if day == 355:
                f10 = f10_tool.get_f10(year-1, day)
            else:
                f10 = f10_tool.get_f10(year, day)
            sin_growth(year=year, season=s, zone=z, f10=f10)
            flat_sin_growth(year=year, season=s, zone=z, f10=f10)
            pre_sin_growth(year=year, season=s, zone=z, f10=f10)
            cnofs_growth(year=year, season=s, zone=z, bubbs=True, f10=f10)
            fejer_growth(year=year, season=s, zone=z, bubbs=True, f10=f10)
