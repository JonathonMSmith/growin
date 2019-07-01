"""collection of premade plot routines for the growth rate
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import irfl
from matplotlib.lines import Line2D


def plot_drifts(ivm=None, year=None, season=None, zone=None, day=None,
                lon=None):
    """plot the vertical drift as a function of local time and the results of
       the ExB drift fourier fit
    """
    fit_drift = []
    for local_time in ivm.drifts.slt:
        vals = ivm.drifts.sel(year=year, season=season, longitude=zone,
                              slt=float(local_time))
        coeffs = vals.coefficients
        ve01 = vals.ve01
        fit_drift.append(irfl.growth_rate.exb_calc(coeffs,
                                                   float(ve01),
                                                   float(local_time)))
    local_time = ivm.drifts.slt.values
    plt.errorbar(local_time,
                 ivm.drifts.drifts.sel(year=year, season=season,
                                       longitude=zone).values,
                 ivm.drifts.dev.sel(year=year, season=season,
                                    longitude=zone).values, color='k')
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
    plt.title(season + '-' + str(year) + zone, fontsize=40)
    fig = plt.gcf()
    fig.set_size_inches(18, 6)
    plt.tight_layout()
    filename = season + '.png'
    file_path = irfl.utils.generate_path('drift', lon=lon, year=year, day=day)
    if not os.path.isdir(file_path):
        os.makedirs(file_path)
    full_name = os.path.join(file_path, filename)
    fig.savefig(full_name, dpi=100)
    plt.close()

def plot_f_peak(sami):
    """attains and plots the F layer peak density altitude
    """
    nz = np.shape(sami.glat)[0]
    nf = np.shape(sami.glat)[1]
    f_peak = []
    for i,t in enumerate(sami.slt):
        apex_alt = []
        apex_deni = []
        for ft in range(nf):
            ind = np.argmax(sami.zalt[:, ft])
            apex_alt.append(sami.zalt[ind, ft])
            apex_deni.append(np.sum(sami.deni[ind, ft, :, i]))
        ind = np.argmax(apex_deni)
        f_peak.append([t, apex_alt[ind], apex_deni[ind]])    
    peak_times = np.array([peak[0] for peak in f_peak])
    peak_alts = np.array([peak[1] for peak in f_peak])
    i, = np.where(peak_times == np.min(peak_times))
    peak_times = np.roll(peak_times, -i)
    peak_alts = np.roll(peak_alts, -i)
    plt.plot(peak_times, peak_alts,  color='k', linestyle='--')

def plot_apex_density(sami):
    nz = np.shape(sami.glat)[0]
    nf = np.shape(sami.glat)[1]
    f_peak = []
    apex_denis = []
    for i,t in enumerate(sami.slt):
        apex_alt = []
        apex_deni = []
        for ft in range(nf):
            ind = np.argmax(sami.zalt[:, ft])
            apex_alt.append(sami.zalt[ind, ft])
            apex_deni.append(np.sum(sami.deni[ind, ft, :, i]))
        apex_denis.append(apex_deni)
    slt = sami.slt
    i, = np.where(slt == np.min(slt))
    slt = np.roll(np.array(slt), -i)
#    alt = np.roll(np.array(apex_alt), -i)
    alt = np.array(apex_alt)
    xx, yy = np.meshgrid(slt, alt)
    apex_denis = np.array(apex_denis).transpose()
    apex_denis = np.roll(apex_denis, -i)
    plt.pcolormesh(xx, yy, apex_denis)
#    cs = plt.contourf(slt, alt, apex_denis, levels=15)
    plt.colorbar()
    plot_f_peak(sami)
#    plt.contour(cs, colors='k')
    plt.tight_layout()
    plt.show()

def plot_bottomside(rt_eq):
    """attains and plots the F layer bottomside altitude
    """
    ind, = np.argmax(rt_eq.K, axis=2)
    alts = rt_eq.K.alt[ind]
    plt.plot(rt_eq.ut, alts, color='r', linestyle='--')

def plot_growth_term(sami, season, term='gamma', vmin=0, vmax=8*10**(-4),
                     bottomside=True, peak=True, drift=False):
    """plot one of the terms or the result of the RT growth rate equation
    """
    rt_eq = sami.gamma.assign_coords(ut=((sami.gamma.ut +
                                          sami.gamma.lon.values / 15) % 24))
    i, = np.where(rt_eq.ut == np.min(rt_eq.ut))
    rt_eq = rt_eq.roll(ut=-i[0])
    if peak:
        plot_f_peak(sami)
    if bottomside:
        plot_bottomside(rt_eq)

    growth_term = rt_eq[term][0]
    g_ax = growth_term.plot(x='ut', y='alt', vmin=vmin, vmax=vmax,
                            cmap='cubehelix_r')
    fig = plt.gcf()
    cbar = g_ax.colorbar
    cbar.set_label(label=r'$\gamma$ $[s^{-1}]$', fontsize=40)
    cbar.ax.tick_params(labelsize=30)
    plt.ylim(210, 600)
    plt.xlabel('solar local time [h]', fontsize=40)
    plt.xticks(np.arange(0, 30, 6), fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.ylabel('altitude [km]', fontsize=40)
    plt.title(season+str(sami.year)+', lon = '+str(sami.lon0), fontsize=40)
    fig.set_size_inches(24, 6)
    if drift:
        host = plt.gca()
        new = host.twinx()
        v_drift = rt_eq['V'][0, :, 0]
        v_drift.plot(x='ut', ax=new)
        new.set_ylim(-75, 75)
        new.yaxis.set_label(None)
        new.yaxis.label.set_color('tab:blue')
        new.tick_params(axis='y', which='major', labelsize=30,
                        color='tab:blue', labelcolor='tab:blue')
        new.set_title(None)
        legend_elements = [Line2D([0], [0], color='tab:blue', label='V [m/s]'),
                           Line2D([0], [0], linestyle='--', color='k',
                                  label='F peak'),
                           Line2D([0], [0], linestyle='--', color='r',
                                  label='bottomside')]
        new.legend(handles=legend_elements, prop={'size': 20})


    plt.tight_layout()
    file_path = irfl.utils.generate_path('growth', lon=sami.lon0,
                                         year=sami.year,
                                         day=sami.day)
    filename = season + '_' + term + '_' + sami.tag + '.png'
    full_name = os.path.join(file_path, filename)
    fig.savefig(full_name, dpi=100)
    plt.close()
