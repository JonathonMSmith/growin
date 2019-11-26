import numpy as np

def get_f10(year, day):
    f10 = np.loadtxt('f10.lst')
    for i, d in enumerate(f10[:, 1]):
        tmp_yr = f10[i, 0]
        yd = tmp_yr*1000 + d
        if yd > (year*1000 + day):
            f10 = np.mean([f10[i, 3], f10[i+1, 3]])
            return f10
