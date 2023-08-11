#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Full license can be found in License.md
# -----------------------------------------------------------------------------

import numpy as np

""" 
Functions
-------------------------------------------------------------------------------
generate_path(tag, lon, year, day)
    Generates path to archived model runs based on input paramters.

set_archive_dir(path=None, store=None)
    Allows user to specify the location where the model outputs will be stored

get_unformatted_data(dat_dir, var_name, nz, nf, ni, nt, reshape=False)
    routine to interpret unformatted binary files created by the SAMI2 model
-------------------------------------------------------------------------------
"""


def generate_path(tag, lon=None, year=None, day=None, test=False, 
                  end_year=None):
    """Creates a path based on run tag, date, and longitude

    Parameters
    ----------
    tag : (string)
        specifies name of model run
    lon : (int)
        longitude of model run
    year : (int)
        year of model run
    day : (int)
        day of year of model run
    test : (bool)
        If True, use directory for test data.  If False, use archive_dir
        (default = False)

    Returns
    -------
    archive_path : (string)
        Complete path pointing to model archive for a given run
    """
    from os import path

    if not isinstance(tag, str):
        raise TypeError

    if test:
        from growin import test_data_dir
        top_directory = test_data_dir
    else:
        from growin import archive_dir
        top_directory = archive_dir

    # Check if top_directory is empty string, ie, user has not specified
    # a directory through set_archive_dir
    if top_directory:
        if end_year:
            str_fmt = '{year:4d}_{end_year:4d}/'
            archive_path = path.join(top_directory, tag,
                                     (str_fmt.format(year=year,
                                                     end_year=end_year)))
        else:
            str_fmt = 'lon{lon:03d}/{year:4d}_{day:03d}/'
            archive_path = path.join(top_directory, tag,
                                     (str_fmt.format(lon=lon,
                                                     year=year,
                                                     day=day)))
    else:
        raise NameError(''.join(('Archive Directory Not Specified: ',
                                 'Run growin.utils.set_archive_dir')))

    return archive_path


def set_archive_dir(path=None, store=True):
    """Set the top level directory growin uses to look for data and reload.

    Parameters
    ----------
    path : string
        valid path to directory growin uses to look for data
    store : bool
        if True, store data directory for future runs
    """
    import os
    import growin

    if os.path.isdir(path):
        if store:
            with open(growin.archive_path, 'w') as archive_file:
                archive_file.write(path)
        growin.archive_dir = path
    else:
        raise ValueError('Path does not lead to a valid directory.')

def fwhm(model, fac):
    """Return the full width half max in both altitude and local time for a
    growth rate calculation

    Parameters
    ----------
    model: xarray dataSet
    fac: int
    factor for time conversion to a 24 hour range
    """
    # first get the max val
    max_growth = model.gamma.max()
    
    # iterate out in four directions until the growth rate falls below max / 2
    rt_eq = model.gamma.assign_coords(lt=((model.gamma.ut / fac +
                                           model.gamma.lon.values / 15) % 24))

    idx = np.where(rt_eq == rt_eq.max())

    half_lts = np.where(rt_eq[0, :, idx[2]] >= max_growth / 2)
    half_alt = np.where(rt_eq[0, idx[1], :] >= max_growth / 2)

    lt_at_max = rt_eq['lt'][idx[1]]
    alt_at_max = rt_eq['alt'][idx[2]]

    lt_fwhm = abs(rt_eq['lt'][half_lts[0][-1]] - rt_eq['lt'][half_lts[0][0]])
    alt_fwhm = rt_eq['alt'][half_alt[0][-1]] - rt_eq['alt'][half_alt[0][0]]


    return max_growth, lt_at_max, alt_at_max, lt_fwhm, alt_fwhm
