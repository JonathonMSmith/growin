"""collection of premade plot routines for the growth rate
"""
import numpy as np
import matplotlib.pyplot as plt
import irfl


def plot_drifts(ivm, year, season, lon):
    """plot the vertical drift as a function of local time and the results of
       the ExB drift fourier fit
    """
    fit_drift = []
    for local_time in ivm.drifts.slt:
        vals = ivm.drifts.sel(year=year, season=season, longitude=lon,
                              slt=float(local_time))
        coeffs = vals.coefficients
        ve01 = vals.ve01
        fit_drift.append(irfl.growth_rate.exb_calc(coeffs,
                                                   float(ve01),
                                                   float(local_time)))
    local_time = ivm.drifts.slt.values
    plt.errorbar(local_time,
                 ivm.drifts.drifts.sel(year=year,
                                       season=season, longitude=lon).values[0],
                 ivm.drifts.dev.sel(year=year,
                                    season=season, longitude=lon).values[0],
                 color='k')
    plt.plot(local_time, fit_drift, color='r', zorder=10)
    plt.plot([0, 12, 24], [0, 0, 0])
    plt.plot([0, 12, 24], [ve01, ve01, ve01], linestyle='-.', color='violet')
    plt.xlim(0, 23)
    plt.ylim(-75, 75)
    plt.xticks(np.arange(0, 30, 6), fontsize=30)
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


def plot_bottomside(sami):
    """attains and plots the F layer bottomside altitude
    """
    nf = np.shape(sami.glat)[1]
    f_peak = []
    for i, t in enumerate(sami.ut):
        apex_alt = []
        apex_deni = []
        for ft in range(nf):
            ind = np.argmax(sami.zalt[:, ft])
            apex_alt.append(sami.zalt[ind, ft])
            apex_deni.append(np.sum(sami.deni[ind, ft, :, i]))
        ind = np.argmax(apex_deni)
        f_peak.append([t, apex_alt[ind], apex_deni[ind]])
    peak_times = np.array([(peak[0] + sami.lon0/15) % 24 for peak in f_peak])
    peak_alts = np.array([peak[1] for peak in f_peak])
    i, = np.where(peak_times == np.min(peak_times))
    peak_times = np.roll(peak_times, -i)
    peak_alts = np.roll(peak_alts, -i)
    plt.plot(peak_times, peak_alts, color='k', linestyle='--')


def plot_f_peak(rt_eq):
    """attains and plots the F layer peak density altitude
    """
    ind, = np.argmax(rt_eq.K, axis=2)
    alts = rt_eq.K.alt[ind]
    plt.plot(rt_eq.ut, alts, color='r', linestyle='--')


def plot_growth_term(sami, season, term='gamma', vmin=0, vmax=8*10**(-4),
                     bottomside=True, peak=True):
    """plot one of the terms or the result of the RT growth rate equation
    """
    if bottomside:
        plot_bottomside(sami)
    rt_eq = sami.gamma.assign_coords(ut=((sami.gamma.ut +
                                          sami.gamma.lon.values / 15) % 24))
    i, = np.where(rt_eq.ut == np.min(rt_eq.ut))
    rt_eq = rt_eq.roll(ut=-i[0])
    if peak:
        plot_f_peak(rt_eq)
    growth_term = rt_eq[term][0]
    fig = growth_term.plot(x='ut', y='alt', vmin=vmin, vmax=vmax,
                           cmap='cubehelix_r')
    fig.colorbar.set_label(label=r'$\gamma$ $[s^{-1}]$', fontsize=40)
    cbar = fig.colorbar
    cbar.ax.tick_params(labelsize=30)
    plt.ylim(210, 600)
    plt.xlabel('solar local time [h]', fontsize=40)
    plt.xticks(np.arange(0, 30, 6), fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.ylabel('altitude [km]', fontsize=40)
    plt.title(season+str(sami.year)+', lon = '+str(sami.lon0), fontsize=40)
    fig = plt.gcf()
    fig.set_size_inches(24, 6)
    plt.tight_layout()
    fig.savefig('gamma_'+season+str(sami.year)+', lon = ' +
                str(sami.lon0)+'.png', dpi=100)
