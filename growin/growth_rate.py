"""
Calculation of altitude profiles of flux tube integrated growth rates of the 
Rayleigh Taylor instability using SAMI2 output for ion densities, and electron
temperatures

Author
-------------------------------------------------------------------------------
Jonathon Smith (JS), 20 Sep 2018, Goddard Space Flight Center (GSFC)
-------------------------------------------------------------------------------
"""
import math
from datetime import datetime, timedelta
import numpy as np
import igrf12
import xarray as xr

"""
CONSTANTS
m_e : electron mass in kg
q_e : electron charge coulombs
m_i : ion/neut mass (amu)
g_0 : gravity at earth's surface m/s**2
R_e : (km) Earth's radius
"""
m_e = 9.1 * 10**(-31)
amu = 1.7 * 10**(-27)
q_e = 1.6 * 10**(-19)
m_i = {'H':1, 'He':4, 'N':14, 'O':16, 'N2':28, 'NO':30, 'O2':32, 'Ar':40}
sami_enumerate = {0:'H', 4:'He', 6:'N', 1:'O', 2:'NO', 5:'N2', 3:'O2'}
g_0 = 9.8
R_e = 6.38 * 10**3

class FluxTube():
    """
    flux tube integrated (fti) quantities
    gr_arr: growth rate
    sig_ratio: ratio of F region conductivity to total conductivity
    V: vertical exb drift
    U: neutral wind in the L direction
    g: gravity at apex altitude
    nu_ef: effective f region collision frequency
    g_nu: g over nu
    K: altitude gradient in flux tube integrated electron density
    alt: apex altitude
    N: electron density
    sig_F: f region conductivity
    sig_total: total flux tube integrated conductivity
    """
    def __init__(self, sami_out, ft, max_alt, exb):
        self.U = 0
        self.N = 0
        self.sig_total = 0
        self.sig_F = 0
        self.nu_ef = 0
        self.alt = max_alt
        self.V = exb
        self.g = g_e_L(sami_out, ft)
        self.gamma = 0
        self.sig_ratio = 0
        self.g_nu = 0
        self.K = 0

class FluxTubeCell():
    """one cell or bin of a flux tube from the sami model and all of its
       attributes are stored in an object for easy reference and use
       this includes all of the density values and the location coordinates
       to specify this cell
    """
    def __init__(self, sami_out, ftl, ft, iyd, d_time, d_str, t_step):
        lat, lat_2, lon, alt, alt_2 = ft_bin_loc(sami_out, ftl, ft)
        mag, atmos, hwm = run_models(sami_out, lat, lon, alt, ftl, ft,
                                     d_str, t_step)
        denis = sami_out.deni[ftl, ft, :, t_step]

        self.alt = alt
        self.len = ft_length(alt, alt_2, lat, lat_2)
        self.n_n, species = get_n_n(atmos)
        self.n_e = np.sum(denis)
        self.t_e = sami_out.te[ftl, ft, t_step]
        self.A = np.mean([m_i[i] for i in species])
        self.B = float(mag.total) * 10**(-9) #convert nT output to T
        self.phi = (90 - float(mag.incl)) * math.pi / 180 #dip angle radians
        self.sig = sigma_tot(denis=denis, n_n=self.n_n, n_e=self.n_e, B=self.B,
                             A=self.A, T_e=self.t_e)
        self.wind = hwm / 10**2 #convert to m/s
        self.nu = 0
        for n_i in denis:
            self.nu += nu_i(n_i, self.n_n, self.A)

def ft_bin_loc(sami_out, ftl, ft):
    """
    returns the location and spatial extent of current bin
    """
    lat = sami_out.glat[ftl, ft]
    lat_2 = sami_out.glat[ftl + 1, ft]
    lon = sami_out.glon[ftl, ft]
    alt = sami_out.zalt[ftl, ft]
    alt_2 = sami_out.zalt[ftl+1, ft]
    return lat, lat_2, lon, alt, alt_2

def format_dates(sami_out, t_step):
    day = sami_out.day
    year = sami_out.year
    ut = sami_out.ut[t_step]
    iyd = int((year - (2000 if year > 1999 else 1900)) * 1000) + day
    d_time = datetime(year, 1, 1) + timedelta(days=day-1, seconds=ut*3600)
    d_str = d_time.strftime('%Y-%m-%d')
    return iyd, d_time, d_str

def ft_length(alt_1, alt_2, lat_1, lat_2):
    """
    law of cosines for determining linear extent of flux tube (ft) bin
    alt_1 & 2 in km
    lat_1 & 2 in radians?
    """
    delta_lat = (lat_2-lat_1) * math.pi / 180
    alt_1 += R_e
    alt_2 += R_e
    lenkm = alt_1**2 + alt_2**2 - 2 * alt_1 * alt_2 * math.cos(delta_lat)
    return math.sqrt(lenkm) * 10**3

def nu_i(n_i, n_n, A):
    """
    approximate calculation of ion collision frequency from Kelley 89
    """
    return 2.6 * 10**(-9) * (n_i + n_n) * A**(-1/2)

def nu_e(n_n, n_e, T_e):
    """
    approximate calculation of electron collision frequency from Kelly 89
    """
    nu_e_n = 5.4 * 10**(-10) * n_n * T_e**(1/2)
    nu_e_i = (34 + 4.18 * math.log(T_e**3 / n_e)) * n_e * T_e**(-3/2)
    return nu_e_n + nu_e_i

def g_e_L(sami_out, ft):
    """
    gravity at the bin altitude
    L is geocentric distance in earth radii
    the L shell at the local altitude not apex.
    """
    apex_alt = np.amax(sami_out.zalt[:, ft])
    L = (apex_alt + R_e) / R_e
    return g_0 / L**2


def omega(B, particle):
    if particle == 'electron':
        return q_e * B / m_e
    else:
        return q_e * B / (m_i[particle] * amu)

def sigma_tot(denis, n_n, n_e, B, A, T_e):
    """
    inputs :
    ion densities from sami cm-3
    total neutral density cm-3
    total electron density cm-3
    magnetic field in Teslas
    A average neutral density in amus
    T_e electron temperature in Kelvin
    output :
    total Pedersen conductivity at location in mho/m
    """
    sig_tot = 0
    k_e = omega(B, 'electron') / nu_e(n_n, n_e, T_e)
    for i, n_i in enumerate(denis):
        if n_i > 0:
            k_i = omega(B, sami_enumerate[i]) / nu_i(n_i, n_n, A)
            sig_tot += n_i * k_i / (1 + k_i**2)
    sig_tot += n_e * k_e / (1 + k_e**2)
    return 10**6 * sig_tot * q_e / B

def get_n_n(nn):
    """
    nn is sami attribute denn (atmos)
    neutral number densities are in m-3 so the must be converted to cm-3
    """
    species = []
    n_n = 0
    for n in range(len(nn)):
        if nn[n] > 0:
            n_n += float(nn[n])
            species.append(sami_enumerate[n])
    return n_n, species

def calc_growth_rate(tube):
    """
    the growth rate equation from Sultan 96
    inputs
    sig_F_P : flux tube integrated pedersen cond. F region in mho
    sig_total : '                                ' total in mho
    V_P : flux tube integrated vertical drift or drift at apex altitude m/s
    U_L : flux tube integrated neutral wind perp. B in L direction in m/s
    g_e : gravtiy at apex altitude in m/s^2
    nu_eff : collision frequency in s-1
    K_F : altitude gradient in density in m-1
    """
    sig_F_P = tube.sig_F
    sig_total = tube.sig_total
    V = tube.V
    U_L = tube.U
    g_e = tube.g
    nu_eff = tube.nu_ef
    K_F = tube.K

    gamma = sig_F_P / sig_total * (V - U_L + g_e/nu_eff) * K_F
    return gamma

def run_models(sami, lat, lon, alt, cell, flux_tube, d_str, t_step):
    mag = igrf12.igrf(d_str, glat=lat, glon=lon, alt_km=alt)
    atmos = sami.denn[cell, flux_tube, :, t_step]
    #only the meridional component of wind is used as per Sultan1996
    hwm = sami.u[cell, flux_tube, t_step]
    return mag, atmos, hwm

def eval_tubes(sami_out, exb, t_step=0):
    """
    calculate the flux tube integrated quantities for each flux tube needed for
    the growth rate calculation
    inputs:
    sami_out : sami2py object
    t_step : array index for sami2py object timestep variable
    """
    iyd, d_time, d_str = format_dates(sami_out, t_step)
    nz = np.shape(sami_out.glat)[0]
    nf = np.shape(sami_out.glat)[1]

    tube_list = []
    #print('progress %:')
    for ft in range(nf):
        #progress = int(100 * ft / range(nf)[-1])
        #print('%03d' % progress)
        max_alt = np.amax(sami_out.zalt[:, ft])
        if max_alt <= 200:
            continue
        tube = FluxTube(sami_out, ft, max_alt, exb)
        for ftl in range(nz-1):
            ftc = FluxTubeCell(sami_out, ftl, ft, iyd, d_time, d_str, t_step)
            #Reimann sum values for total flux tube
            tube.U += ftc.wind * math.cos(ftc.phi) * ftc.len * ftc.sig
            tube.N += ftc.n_e * ftc.len * 10**2
            tube.sig_total += ftc.sig * ftc.len
            #Reimann sum values for F region
            if ftc.alt > 200:
                tube.sig_F += ftc.sig * ftc.len
                tube.nu_ef += ftc.nu * ftc.n_e * ftc.len * 10**2

        tube.U = tube.U / tube.sig_total
        tube.nu_ef = tube.nu_ef / tube.N
        tube.sig_ratio = tube.sig_F / tube.sig_total
        tube.g_nu = tube.g / tube.nu_ef
        tube_list.append(tube)
    return tube_list, d_time

def rt_growth_rate(sami_out, exb, t_step=0):
    """
    calculate flux tube integrated electron density altitude gradient
    and flux tube integrated growth rate
    """
    tube_list, d_time = eval_tubes(sami_out, exb, t_step)
    nf = len(tube_list)
    for ft in range(nf):
        if ft + 1 == nf:
            h1 = tube_list[ft - 1].alt
            h2 = tube_list[ft].alt
            N_e_1 = tube_list[ft - 1].N
            N_e_2 = tube_list[ft].N
        else:
            h1 = tube_list[ft].alt
            h2 = tube_list[ft + 1].alt
            N_e_1 = tube_list[ft].N
            N_e_2 = tube_list[ft + 1].N

        dN_e = (N_e_2 - N_e_1)
        dh = (h2 - h1) * 10**3
        K = (1 / N_e_1) * (dN_e / dh)
        tube_list[ft].K = K

        gam = calc_growth_rate(tube_list[ft])
        tube_list[ft].gamma = gam
    return tube_list, d_time

def exb_calc(coefficients, ve01, t):
    exb = ve01
    for i, term in enumerate(coefficients):
        a = term[0]
        b = term[1]
        exb += ((a * np.cos((i+1) * t * np.pi / 12))
              + (b * np.sin((i+1) * t * np.pi / 12)))
    return exb

def run_growth_calc(sami, coefficients=None, ve01=0):
    ''' runs the growth rate calculation for a sami2 run
        Parameters:
        sami: SAMI2py model object
        coefficients: 10x2 matrix with fourier coefficients for drift fit
        ve01: offset or 0th term of fourier fit
    '''
    time_steps = len(sami.ut)
    rtgr_sets = []
    if coefficients is None:
        coefficients = np.zeros((10,2))
    lon0 = sami.lon0
    for i in range(time_steps):
        t = sami.ut[i]
        print(str(t)[:3])
        lt = t + lon0/15
        lt = lt % 24
        exb = exb_calc(coefficients, ve01, lt)
        tube_list, t = rt_growth_rate(sami_out=sami, exb=exb, t_step=i)
        tubes = []
        tube_dict = {}
        for tube in tube_list:
            tubes.append(tube.__dict__)
        for k in tubes[0]:
            tube_dict[k] = list(q[k] for q in tubes)
        coords = [('alt', tube_dict['alt'])]
        tmp_dict = {}
        for key in tube_dict:
            if key == 'alt':
                continue
            tmp_dict[key] = xr.DataArray(tube_dict[key], coords)
        rtgr_sets.append(xr.Dataset(tmp_dict))
    rtgr_arr = xr.concat(rtgr_sets, 'ut')
    rtgr_arr = rtgr_arr.assign_coords(ut=(sami.ut))
    rtgr_arr = rtgr_arr.expand_dims('lon')
    rtgr_arr = rtgr_arr.assign_coords(lon=([lon0]))
    return rtgr_arr
