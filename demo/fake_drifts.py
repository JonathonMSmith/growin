import growin
import numpy as np
import matplotlib.pyplot as plt


def sinusoidal_drifts(amplitude=30):
    coef = np.array([[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])

    drift = []
    time = np.linspace(0, 24, 49)

    for t in time:
    
        d = growin.growth_rate.exb_calc(coef, 0 ,t)
        drift.append(d)
    
    drift = np.array(drift)
    drift = amplitude*drift
    return drift


def sinusoidal_pre_drifts(amplitude=30):
    drift = sinusoidal_drifts()
    drift[37] = 10
    drift[38] = 20
    drift[39] = np.mean([20, drift[40]])
    return drift


def fit_drift(drift):
    offset, fit_coef = growin.fourier_exb.fourier_fit(time, drift, 10)
    fit_drift = []
    for t in time:
    
        d = growin.growth_rate.exb_calc(fit_coef, 0 ,t)
        fit_drift.append(d)
    
    fit_drift = np.array(fit_drift)
    return fit_drift
