import tirtil 
import numpy as np
import pickle
import xarray as xr
import matplotlib.pyplot as plt
import seaborn as sns
#year-1,11,1     Dec Sol
#year,1,31
#
#year,2,1        Mar Eq
#year,4,30
#
#year,5,1        Jun Sol
#year,7,31
#
#year,8,1        Sep Eq
#year,10,31
#static touple for seasonal iteration (month, day, month, day)
ivm = sami2py.DriftInstrument(platform='cnofs', name='ivm', clean_level='none')

season_bounds = np.array([11, 2, 5, 8, 11])
#season_bounds = np.array([2, 5, 8])
season_names = ['DecSol', 'MarEq', 'JunSol', 'SepEq']
season_days = {'DecSol':355, 'MarEq':80, 'JunSol':155, 'SepEq':266}
#longitude regions shifted by +15 degrees for ease of binning
longitude_bounds = np.array([0, 75, 145, 300, 360])
#longitude_bounds = np.array([300, 360])
lon_sector_names = ['African', 'Indian', 'Pacific', 'South American']
#midpoints from longitude_bounds, subtract 15 again for sami
model_sectors = {'African':22, 'Indian':95, 'Pacific':207, 
                 'South American':315}

def shift_local_time(inst):
    idx, = np.where(inst['slt'] < 12.)
    inst[idx, 'slt'] += 24.
    return
#shift longitude values 15 degrees for ease of longitude sector binning
def shift_longitude(inst):
    inst['glon'] += 15.
    return
def drift_fix(inst):
    inst['ionVelocityZ'] *= -1
    return
def filter_ivm(inst):
    idx, = np.where((inst['Ni'] >= 3000) & (inst['apex_altitude'] <= 550))
    inst.data = inst.data.iloc[idx]
    return

def plot_drifts(ivm, year, season, lon):
    fit_drift = []
    for lt in ivm.drifts.slt:
        vals = ivm.drifts.sel(year=year, season=season, longitude=lon, 
                              slt=float(lt))
        coeffs = vals.coefficients
        ve01 = vals.ve01
        fit_drift.append(sami2py.growth_rate.exb_calc(coeffs.values[0], 
                         float(ve01), float(lt)))
    lt = ivm.drifts.slt.values
    plt.errorbar(lt, 
                 ivm.drifts.drifts.sel(year=year, 
                                       season=season, longitude=lon).values[0], 
                 ivm.drifts.dev.sel(year=year, 
                                    season=season, longitude=lon).values[0], 
                 color='k')
    plt.plot(lt, fit_drift, color='r', zorder = 10)
    plt.plot([0,12,24],[0,0,0])
    plt.plot([0,12,24],[ve01,ve01,ve01], linestyle='-.', color='violet')
    plt.xlim(0,23)
    plt.ylim(-75, 75)
    plt.xticks(np.arange(0, 30, 6),fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.xlabel(r'solar local time [h]', fontsize=40)
    plt.ylabel(r'vertical drift [ms$^{-1}]$', fontsize=40)
    plt.title(season + '-' + str(year) + lon, fontsize=40)
    fig = plt.gcf()
    fig.set_size_inches(18, 6)
    plt.tight_layout()
    filename = season + '-' + str(year) + '-' + lon + '.png'
    fig.savefig(filename, dpi=100)

def plot_growth_term(S, season, term='gamma', vmin=0, vmax=8*10**(-4)):
    #plot F peak
    nz = np.shape(S.glat)[0]
    nf = np.shape(S.glat)[1]
    f_peak = []
    for i,t in enumerate(S.ut):
        apex_alt = []
        apex_deni = []
        for ft in range(nf):
            ind = np.argmax(S.zalt[:, ft])
            apex_alt.append(S.zalt[ind,ft])
            apex_deni.append(np.sum(S.deni[ind,ft,:,i]))
        ind = np.argmax(apex_deni)
        f_peak.append([t, apex_alt[ind], apex_deni[ind]])    
    peak_times = np.array([(peak[0] + S.lon0/15)%24 for peak in f_peak])
    peak_alts = np.array([peak[1] for peak in f_peak])
    i, = np.where(peak_times == np.min(peak_times))
    peak_times = np.roll(peak_times, -i)
    peak_alts = np.roll(peak_alts, -i)
    plt.plot(peak_times, peak_alts,  color='k', linestyle='--')
    #plot bottomside
    s = S.gamma.assign_coords(ut=((S.gamma.ut + S.gamma.lon.values / 15)%24))
    i, = np.where(s.ut == np.min(s.ut))
    s = s.roll(ut=-i[0])
    ind, = np.argmax(s.K, axis=2)
    alts = s.K.alt[ind]
    plt.plot(s.ut, alts, color='r', linestyle='--')
    #plot growth rates
    growth_term = s[term][0]
    fig = growth_term.plot(x = 'ut', y='alt', vmin=vmin, vmax=vmax, 
                   cmap='cubehelix_r')
    fig.colorbar.set_label(label=r'$\gamma$ $[s^{-1}]$', fontsize=40)
    cbar = fig.colorbar
    ax = plt.gca()
    cbar.ax.tick_params(labelsize=30)
    plt.ylim(210, 600)
    plt.xlabel('solar local time [h]', fontsize=40)
    plt.xticks(np.arange(0, 30, 6),fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.ylabel('altitude [km]', fontsize=40)
    plt.title(season+str(S.year)+', lon = '+str(S.lon0), fontsize=40)
    fig = plt.gcf()
    fig.set_size_inches(24, 6)
    plt.tight_layout()
    fig.savefig('gamma_'+season+str(S.year)+', lon = '+str(S.lon0)+'.png', 
                dpi=100)


ivm.custom.add(drift_fix, 'modify')
ivm.custom.add(shift_longitude, 'modify')
ivm.custom.add(filter_ivm, 'modify')

ivm.exb_fourier_fit(drift_key = 'ionVelocityZ', lon_bins = longitude_bounds, 
                    season_bins = season_bounds, season_names = season_names, 
                    zone_labels = lon_sector_names, 
                    start_year = 2008, stop_year = 2008)
drift_filename = 'drifts_09-12.p'
pickle.dump(ivm, open(drift_filename, 'wb'))
for year in range(2008, 2009):    
    for season in season_names:
        day = season_days[season]
        for lon in lon_sector_names:
            lon_deg = model_sectors[lon]
            ExB_drifts = ivm.drifts.coefficients.sel(year=year, 
                                                     season=season, 
                                                     longitude=lon).values
            ve01 = ivm.drifts.ve01.sel(year=year, 
                                       season=season,longitude=lon).values
            sami2py.run_model(day=day, year=year, lon=lon_deg, fejer=False, 
                              ExB_drifts=ExB_drifts, ve01=ve01, outn=True)
            S = sami2py.model(tag='test', day=day, 
                              year=year, lon=lon_deg, outn=True)
            S.growth_rate(ExB_drifts, ve01)
            sami_filename = 'S_'+year+'_'+season+'_'+str(lon)+'.p'
            pickle.dump(S, open(sami_filename, 'wb'))
            plot_drifts(ivm, year, season, lon)
            plot_growth(S, season)
