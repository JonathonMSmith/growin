"""functions and variables compute the rtgr using the sami3 model
"""
import os
import datetime
import numpy as np
import warnings
import xarray as xr
from . import growth_rate as gr
from . import utils

def get_growth(sami_filename, lon):
    '''get the sami instrument with growth rates calculated
       checks if there is an existing sami instrument with the appropriate tag
       and loads it. Otherwise it runs the growth rate calculation.

    Parameters
    ----------
    sami_filename : (str)
        full path and name of sami3 file
    lon : (int)
        geo longitude in degrees for SAMI run
    '''

    if os.path.isfile(sami_filename):
        sami = xr.load_dataset(sami_filename)
        return sami
    else:
        raise NameError("Invalid Filname, the file:" + sami_filename +
                        "does not exist")

    # TODO select just the longitude provided in the args
    sami = sami[lon]
    # TODO perpare the variable and coordinate names to be compatible
    sami.gamma = gr.run_growth_calc(sami)

    return sami

