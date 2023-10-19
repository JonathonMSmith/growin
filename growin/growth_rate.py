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
import igrf # switching from igrf12
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
    """flux tube integrated (fti) quantities
       all of these have the same shape that is the number of altitude bins by
       the number of time steps

    Object variables
    ----------------
    gr_arr:
        growth rate
    sig_ratio:
        ratio of F region conductivity to total conductivity
    V:
        vertical exb drift
    U:
        neutral wind in the L direction
    g:
        gravity at apex altitude
    nu_ef:
        effective f region collision frequency
    g_nu:
        g over nu
    K:
        altitude gradient in flux tube integrated electron density
    alt:
        apex altitude
    N:
        electron density
    sig_F:
        f region conductivity
    sig_total:
        total flux tube integrated conductivity
    """
    def __init__(self, sami_data, ft, max_alt, exb):
        '''
        Parameters
        ----------
        sami_data : (xarray.core.dataset.Dataset)
            the sami model being used to calculate growth rates
        ft : (int)
            'flux tube' flux tube index for SAMI
        max_alt : (float)
            maximum altitude in km of the flux tube AKA apex altitude
        exb : (float)
            vertical drift at the apex altitude
        '''

        self.U = 0
        self.N = 0
        self.sig_total = 0
        self.sig_F = 0
        self.nu_ef = 0
        self.alt = max_alt
        self.V = exb
        self.g = g_e_L(sami_data, ft)
        self.gamma = 0
        self.sig_ratio = 0
        self.g_nu = 0
        self.K = 0
        self.R = 0

class FluxTubeCell():
    """one cell or bin of a flux tube from the sami model and all of its
       attributes are stored in an object for easy reference and use
       this includes all of the density values and the location coordinates
       to specify this cell
    """
    def __init__(self, sami_data, ftl, ft, iyd, d_str, t_step):
        '''
        Parameters
        ----------
        sami_data : (xarray.core.dataset.Dataset)
            the sami model being used to calculate growth rates
        ftl : (float)
            'flux tube length' index along the length of a flux tube for SAMI
        ft : (int)
            'flux tube' flux tube index for SAMI
        iyd : (int)
            'integer year day' year and day integer like yyyyddd
        d_str : (str)
            'date string' fromat %Y-%m-%d
        t_step : (int)
            'time step'
        '''
        lat, lat_2, lon, alt, alt_2 = ft_bin_loc(sami_data, ftl, ft)
        mag, atmos, hwm, denis, nus = run_models(sami_data, lat, lon, alt, ftl,
                                                 ft, d_str, t_step)
        self.alt = alt
        self.len = ft_length(alt, alt_2, lat, lat_2)
        self.n_n, species = get_n_n(atmos)
        self.n_e = np.sum(denis)
        self.t_e = sami_data.te.values[ftl, ft, t_step]
        self.A = np.mean([m_i[i] for i in species])
        self.B = float(mag.total) * 10**(-9) #convert nT output to T
        self.phi = (90 - float(mag.incl)) * math.pi / 180 #dip angle radians
        self.wind = hwm / 10**2 #convert to m/s

        # TODO: test this collision frequency code with sami3 data that has nu
        self.sig = sigma_tot(denis=denis, n_n=self.n_n, n_e=self.n_e, B=self.B,
                             A=self.A, T_e=self.t_e, nus=nus)
        if nus is None:
            self.nu = 0
            for n_i in denis:
                self.nu += nu_i(n_i, self.n_n, self.A)
        else:
            self.nu = np.sum(nus)

        self.r_local = r_local(denis, alt)

def ft_bin_loc(sami_data, ftl, ft):
    """returns the location and spatial extent of current bin

    Parameters
    ----------
    sami_data : (xarray.core.dataset.Dataset)
        the sami model being used to calculate growth rates
    ftl : (float)
        'flux tube length' index along the length of a flux tube for SAMI
    ft : (int)
        'flux tube' flux tube index for SAMI
    """
    lat = sami_data.glat.values[ftl, ft]
    lat_2 = sami_data.glat.values[ftl + 1, ft]
    lon = sami_data.glon.values[ftl, ft]
    alt = sami_data.zalt.values[ftl, ft]
    alt_2 = sami_data.zalt.values[ftl+1, ft]
    return lat, lat_2, lon, alt, alt_2

def format_dates(sami, t_step):
    '''returns the date in all the required formats for different packages

    Parameters
    ----------
    sami : (sami2py.Model)
        the sami model being used to calculate growth rates
    t_step : (int)
        time step for the sami model object
    '''
    day = int(sami.day.values)
    year = sami.year.values
    ut = sami.ut[t_step].values
    iyd = int((year - (2000 if year > 1999 else 1900)) * 1000) + day
    d_time = datetime(year, 1, 1) + timedelta(days=day-1, seconds=ut*3600)
    d_str = d_time.strftime('%Y-%m-%d')
    return iyd, d_time, d_str

def ft_length(alt_1, alt_2, lat_1, lat_2):
    """law of cosines for determining linear extent of flux tube (ft) bin

    Parameters
    ----------
    alt_1 : (float)
        altitude of current flux tube cell in km
    alt_2 : (float)
        altitude of next flux tube cell on same flux tube in km
    lat_1 : (float)
        latitude of current flux tube cell in degrees
    lat_2 : (float)
        latitude of next flux tube cell on same flux tube in degrees
    """
    delta_lat = (lat_2-lat_1) * math.pi / 180
    alt_1 += R_e
    alt_2 += R_e
    lenkm = alt_1**2 + alt_2**2 - 2 * alt_1 * alt_2 * math.cos(delta_lat)
    return math.sqrt(lenkm) * 10**3

def nu_i(n_i, n_n, A):
    """approximate calculation of ion collision frequency from Kelley 89

    Parameters
    ----------
    n_i : (float)
        ion density cm-3
    n_n : (float)
        neutral density cm-3
    A : (int)
        mean neutral molecular mass in atomic mass units
    """
    return 2.6 * 10**(-9) * (n_i + n_n) * A**(-1/2)

def nu_e(n_n, n_e, T_e):
    """approximate calculation of electron collision frequency from Kelly 89

    Parameters
    ----------
    n_n : (float)
        neutral density cm-3
    n_e : (float)
        electron density cm-3
    T_e : (float)
        electron temperature K
    """
    nu_e_n = 5.4 * 10**(-10) * n_n * T_e**(1/2)
    nu_e_i = (34 + 4.18 * math.log(T_e**3 / n_e)) * n_e * T_e**(-3/2)
    return nu_e_n + nu_e_i

def g_e_L(sami_data, ft):
    """gravity at the bin altitude
       L is geocentric distance in earth radii
       the L shell at the local altitude not apex.

    Parameters
    ----------
    sami_data : (xarray.dataset.Dataset)
        the sami model being used to calculate growth rates       
    ft : (int)
        'flux tube' flux tube index for SAMI
    """
    apex_alt = np.amax(sami_data.zalt[:, ft])
    L = (apex_alt + R_e) / R_e
    return g_0 / L**2

def r_local(denis, alt):
    """Local recombination from Sultan eq 21 alpha*n_mol
       Risbeth & Garriott '69: Dissociative recombination is the principal
       E and F region loss mechansim.
       Huba '96: RTI not damped by recombination in F region
       n_mol is the concentration of molecular ions
       alpha = 2*10**(-7) according to Sultan '92
       alpha ~ 10**(-7) according to Risbeth & Garriott '69
    """
    n_mol = 0
    if alt < 200:
        return 0
    for i, n_i in enumerate(denis):
        if i == 2 | i == 3 | i == 5:
            n_mol += n_i
    return n_mol*2*10**(-7)


def omega(B, particle):
    '''
    Parameters
    ----------
    B : (float)
        total electron density cm-3
    particle : (str)
        particle name to get the correct gyrofrequency
    '''
    if particle == 'electron':
        return q_e * B / m_e
    else:
        return q_e * B / (m_i[particle] * amu)

def sigma_tot(denis, n_n, n_e, B, A, T_e):
    """calculate thetotal Pedersen conductivity at location in mho/m

    Parameters
    ----------
    denis : (list)
        ion densities from sami cm-3
    n_n : (float)
        total neutral density cm-3
    n_e : (float)
        total electron density cm-3
    B : (float)
        magnetic field in Teslas
    A : (float)
        average neutral density in amus
    T_e:
        electron temperature in Kelvin
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
    Parameters
    ----------
    nn : (array)
        neutral number densities in cm-3
        there was an old note here that they were in m-3, but this was likely
        before this quantity was build into sami2py
    """
    species = []
    n_n = 0
    for n in range(len(nn)):
        if nn[n] > 0:
            n_n += float(nn[n])
            species.append(sami_enumerate[n])
    return n_n, species

def calc_growth_rate(tube):
    """the growth rate equation from Sultan 96

    Parameters
    ----------
    tube : (FluxTube)
        flux tube object

    Variables used
    ----------
    sig_F_P : (float)
        flux tube integrated pedersen cond. F region in mho
    sig_total : (float)
        flux tube integrated pedersen cond. total in mho
    V_P : (float)
        flux tube integrated vertical drift or drift at apex altitude m/s
    U_L : (float)
        flux tube integrated neutral wind perp. B in L direction in m/s
    g_e : (float)
        gravtiy at apex altitude in m/s^2
    nu_eff : (float)
        collision frequency in s-1
    K_F : (float)
        altitude gradient in density in m-1
    """
    sig_F_P = tube.sig_F
    sig_total = tube.sig_total
    V = tube.V
    U_L = tube.U
    g_e = tube.g
    nu_eff = tube.nu_ef
    K_F = tube.K
    R_T = tube.R

    gamma = sig_F_P / sig_total * (V - U_L + g_e/nu_eff) * K_F - R_T
    return gamma

def run_models(sami, lat, lon, alt, cell, flux_tube, d_str, t_step):
    '''run all required models to get quantities not contained in SAMI2
       vestigial inclusion of neutral density here from before SAMI2 offered it

    Parameters
    ----------
    sami : (sami2py.Model)
        sami2 model output
    lat : (float)
        latitude where model is to be run
    lon : (float)
        longitude where model is to be run
    alt : (float)
        altitude where model is to be run
    cell : (int)
        index for the cell along the flux tube
    flux_tube : (int)
        index for the flux tube
    d_str : (string)
        string versionn of date
    t_step : (int)
        time step for sami2
    '''
    mag = igrf.igrf(d_str, glat=lat, glon=lon, alt_km=alt)
    hwm = sami.u4.values[cell, flux_tube, t_step]

    if 'denn' in sami.data_vars:
        atmos = sami.denn.values[cell, flux_tube, :, t_step]
        denis = sami.deni.values[cell, flux_tube, :, t_step]
    else:
        atmos = []
        denis = []
        for i1 in range(1, 8):
            n_var_name = ''.join(['denn', str(i1)])
            i_var_name = ''.join(['deni', str(i1)])
            atmos.append(sami[n_var_name][cell, flux_tube, t_step].values)
            denis.append(sami[i_var_name][cell, flux_tube, t_step].values)
    if 'nuin' in sami.data_vars:
        nus = []
        for i1 in range(1, 8):
            nu_var_name = ''.join(['nuin', str(i1)])
            nus.append(sami[nu_var_name][cell, flux_tube, t_step].values)
    else:
        nus = None

    #only the meridional component of wind is used as per Sultan1996
    return mag, atmos, hwm, denis, nus

def eval_tubes(sami, exb, t_step=0):
    """calculate the flux tube integrated quantities for each flux tube needed
       for the growth rate calculation

    Parameters
    ----------
    sami : (sami2py.Model)
        sami2py model object
    t_step : (int)
        array index for sami2py object timestep variable
    """
    if not isinstance(sami, xr.core.dataset.Dataset):
        sami_data = sami.data
    else:
        sami_data = sami
    iyd, d_time, d_str = format_dates(sami, t_step)
    nz = sami_data.nz.shape[0]
    nf = sami_data.nf.shape[0]
    tube_list = []
    for ft in range(nf):
        max_alt = np.amax(sami_data.zalt.values[:, ft])
        if max_alt <= 200:
            continue
        if max_alt > 650:
            continue
        tube = FluxTube(sami_data, ft, max_alt, exb)
        for ftl in range(nz-1):
            ftc = FluxTubeCell(sami_data, ftl, ft, iyd, d_str, t_step)
            #Reimann sum values for total flux tube
            tube.U += ftc.wind * math.cos(ftc.phi) * ftc.len * ftc.sig
            #factor of 100 to convert ftc.len from m to cm
            tube.N += ftc.n_e * ftc.len * 10**2
            tube.sig_total += ftc.sig * ftc.len
            tube.R += ftc.r_local * ftc.n_e * ftc.len * 10**2
            #Reimann sum values for F region
            if ftc.alt > 200:
                tube.sig_F += ftc.sig * ftc.len
                tube.nu_ef += ftc.nu * ftc.n_e * ftc.len * 10**2

        #U is weighted by the total flux tube integrated Pedersen conductivity
        tube.U = tube.U / tube.sig_total
        tube.sig_ratio = tube.sig_F / tube.sig_total
        #nu_ef and R are weighted by the flux tube integrated electron density
        tube.nu_ef = tube.nu_ef / tube.N
        tube.R = tube.R / tube.N
        tube.g_nu = tube.g / tube.nu_ef
        tube_list.append(tube)
    return tube_list, d_time

def rt_growth_rate(sami, exb, t_step=0):
    """calculate flux tube integrated electron density altitude gradient
       and flux tube integrated growth rate. These are done together to avoid
       iterating over all of the flux tubes twice to do each of these
       calculations individually. Further this is because K cannot be
       calculated using this method in the eval_tubes function.

    Parameters
    ----------
    sami_out : (sami2py.Model)
    exb : (float)
        vertical plasma drift at the apex of flux tube
    t_step : (int)
        time step for SAMI2
    """
    tube_list, d_time = eval_tubes(sami, exb, t_step)
    nf = len(tube_list)
    """
    calculate K by taking a simple gradient of N
    h1/N_e_1 are the values lower in altitude
    h2/N_e_2 are the values higher in altitude
    """
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

        #with all the requisite variables calculate the growth rate
        gam = calc_growth_rate(tube_list[ft])
        tube_list[ft].gamma = gam
    return tube_list, d_time

def exb_calc(coefficients, ve01, t):
    '''
    Parameters
    ----------
    coefficients : (array)
        10x2 array of Fourier coefficients
    ve01 : (float)
        flat offset for fourier function 0 by default
    t : (float)
        time in hours
    '''
    exb = ve01
    for i, term in enumerate(coefficients):
        a = term[0]
        b = term[1]
        exb += ((a * np.cos((i+1) * t * np.pi / 12))
              + (b * np.sin((i+1) * t * np.pi / 12)))
    return exb

def run_growth_calc(sami, coefficients=None, ve01=0):
    '''runs the growth rate calculation for a sami2 run. Requires external
       drift information until exb drifts from sami2 are an archived data prod.

    Parameters
    ----------
    sami : (sami2py.Model)
        sami2py model object
    coefficients : (array)
        10x2 array of fourier coefficients describe the vertical drift function
    ve01: (float)
        offset or 0th term of fourier fit
    '''
    time_steps = len(sami.ut)
    rtgr_sets = []
    if coefficients is None:
        coefficients = np.zeros((10, 2))
    if "lon0" in sami.attrs:
        lon0 = sami.lon0
    else:
        lon0 = np.mean(sami.glon)

    for i in range(time_steps):
        t = sami.ut[i]
        print(i)
        print(str(t))
        lt = t + lon0/15
        lt = lt % 24
        if coefficients is None:
            exb = sami['u1p'][i, 5, 80]
        else:
            exb = exb_calc(coefficients, ve01, lt)
        tube_list, t = rt_growth_rate(sami=sami, exb=exb, t_step=i)
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
