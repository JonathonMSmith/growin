'''
year-1,11,1     Dec Sol
year,1,31

year,2,1        Mar Eq
year,4,30

year,5,1        Jun Sol
year,7,31

year,8,1        Sep Eq
year,10,31
'''
# default values to organize seasons and zones
# months that bound each season
SEASON_BOUNDS = np.array([11, 2, 5, 8, 11])
SEASON_NAMES = ['DecSol', 'MarEq', 'JunSol', 'SepEq']
# midpoint days for each season for sami
SEASON_DAYS = {'DecSol': 355, 'MarEq': 80, 'JunSol': 155, 'SepEq': 266}

# longitude regions shifted by +15 degrees for ease of binning
LONGITUDE_BOUNDS = np.array([0, 75, 145, 300, 360])
LON_SECTOR_NAMES = ['African', 'Indian', 'Pacific', 'South American']
# midpoints from longitude_bounds, subtract 15 again for sami
MODEL_SECTORS = {'African': 22, 'Indian': 95, 'Pacific': 207,
                 'South American': 315}


