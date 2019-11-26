"""collection of premade plot routines for the growth rate
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import growin
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
        fit_drift.append(growin.growth_rate.exb_calc(coeffs,
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
    file_path = growin.utils.generate_path('drift', lon=lon, year=year, day=day)
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
    plt.plot(rt_eq.ut, alts, color='k', linestyle='--')

def plot_bubble_histogram(host, sami, months, lons):
    '''
        Parameters:
        host : matplotlib.pyplot.axes
        sami : sami2py.Model
        months : list
            monthly extent of bubbles to be used for histogram
        lons : list
            zonal extent of bubbles to be used for histograms
    '''
    import bubbles
    out_bub, out_freq = bubbles.make_bubble_histogram(sami.year,
                                                      months, lons,
                                                      frequency_type='orbit')
    bins = np.linspace(0, 24, 25)
    center = (bins[:-1] + bins[1:])/2
    #if you want to scale the y axis to show relative occurrence in season
    new = host.twinx()
    #factor of 4 keeps the hist in the bottom quarter of the figure
    new.set_ylim(0, 4)
    new.yaxis.set_label(None)
    new.set_yticks([0, .9])
    new.tick_params('y', labelsize=30, labelcolor='purple')
#    new.tick_params(axis='y', which='both', left=False, right=False,
#                    labelright=False, labelleft=False)
    new.set_title(None)
    new.bar(center, out_bub/out_freq, align='center', width=1, alpha=1,
            color='purple')


def plot_growth_term(sami, season, term='gamma', vmin=-8*10**(-4),
                     vmax=8*10**(-4), bottomside=True, peak=True, drift=False,
                     bubbs=False, save=True, bubble_scale=None, 
                     cmap='cividis', months=None, lons=None):
    """plot one of the terms or the result of the RT growth rate equation
    """
    rt_eq = sami.gamma.assign_coords(ut=((sami.gamma.ut +
                                          sami.gamma.lon.values / 15) % 24))
    i, = np.where(rt_eq.ut == np.min(rt_eq.ut))
    rt_eq = rt_eq.roll(ut=-i[0], roll_coords=True)
    fig = plt.figure()

    if peak:
        plot_f_peak(sami)
    if bottomside:
        plot_bottomside(rt_eq)
    
    host = plt.gca()

    growth_term = rt_eq[term].values[0]
    growth_term = growth_term.transpose()
    ut = rt_eq.ut.values
    alt = rt_eq.alt.values
    x, y = np.meshgrid(ut, alt)
    if (np.min(growth_term) < 0) & (np.max(growth_term) > 0):
        if (vmin is None) & (vmax is None):
            vmax = np.percentile(growth_term, 98)
            vmin = -1*vmax
        cmap = 'RdBu_r' 
    g_ax = host.pcolormesh(x, y, growth_term, vmin=vmin, vmax=vmax, cmap=cmap)
    plt.ylim(np.min(alt), 600)
    plt.xlim(np.min(ut), np.max(ut))
    plt.xlabel('solar local time [h]', fontsize=40)
    plt.xticks([1, 3, 6, 12, 18, 21, 23], fontsize=30)
    plt.yticks([200, 300, 400, 500, 600], fontsize=30)
    plt.ylabel('altitude [km]', fontsize=40)
    plt.title(season+'-'+str(sami.year)+', lon:'+str(sami.lon0)+', F10.7:'+
              str(sami.MetaData['F10.7']),
              fontsize=40)
    fig.set_size_inches(16, 9)

    cbar = fig.colorbar(g_ax, ax=host, pad=.1)
    if term == 'gamma':
        cbar.set_label(label=r'$\gamma$ $[s^{-1}]$', fontsize=40)
    else:
        cbar.set_label(label=term, fontsize=40)
    if vmin == -8*10**(-4):
        cbar.set_ticks([-.0008, -.0004, 0, .0004, .0008])
        cbar.ax.set_yticklabels(['-8e-4', '-4e-4', '0', '4e-4', '8e-4'])
    cbar.formatter.set_powerlimits((0, 0))
    cbar.ax.tick_params(labelsize=30)
    cbar.ax.yaxis.get_offset_text().set_fontsize(30)
    if bubbs:
        '''
            this plots bubble histograms at the bottom of the growth figure
        '''
        #host = plt.gca()
        plot_bubble_histogram(host, sami, months, lons)
        host.set_ylim(40, 600)
        print('allegedly plotted a histogram')
    if drift:
        '''
            this plots drifts over the growth rate
        '''
        new = host.twinx()
        v_drift = rt_eq['V'][0, :, 0].values
        ut = rt_eq.ut.values
        new.plot(ut, v_drift)
        new.plot([0, 12, 24], [0, 0, 0], color='tab:blue', linestyle='-.')
        new.set_ylim(-140, 75)
        new.set_yticks([-60, -30, 0, 30, 60])
        new.yaxis.set_label(None)
        new.yaxis.label.set_color('tab:blue')
        new.tick_params(axis='y', which='major', labelsize=30,
                        color='tab:blue', labelcolor='tab:blue')
        new.set_title(None)
        #legend_elements = [Line2D([0], [0], color='tab:blue', label='V [m/s]'),
        #                   Line2D([0], [0], linestyle='--', color='k',
        #                          label='bottomside')]
        #new.legend(handles=legend_elements, prop={'size': 20})
    plt.tight_layout()
    if save:
        file_path = growin.utils.generate_path('growth', lon=sami.lon0,
                                               year=sami.year,
                                               day=sami.day)
        filename = season + '_' + term + '_' + sami.tag + '.png'
        full_name = os.path.join(file_path, filename)
        fig.savefig(full_name, dpi=100)
        plt.close()
    else:
        plt.show()
