#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2019, Jonathon Smith
# Full license can be found in License.md
# -----------------------------------------------------------------------------
"""
growin
-----------
"""

import os
import sys
import logging

__version__ = str('0.1a1')

# get home directory
home_dir = os.path.expanduser('~')
# get environment name
env_name = os.path.split(sys.prefix)[-1]
growin_dir = os.path.join(home_dir, '.growin', env_name)
# make sure a growin directory for model output exists
if not os.path.isdir(growin_dir):
    # create directory
    os.mkdir(growin_dir)
    print(''.join(('Created .growin directory in user home directory to',
                   'store settings.')))


archive_path = os.path.join(growin_dir, 'archive_path.txt')
if os.path.isfile(archive_path):
    # load up stored data path
    with open(archive_path, 'r') as f:
        archive_dir = f.readline()
else:
    # create file
    with open(archive_path, 'w+') as f:
        f.write('')
    archive_dir = ''
    print(''.join(('Run growin.utils.set_archive_dir to set the path to',
                   ' top-level directory for model outputs.')))

# load test_data directory
with open(os.path.join(growin_dir, 'test_data_path.txt'), 'r') as f:
    test_data_dir = f.readline()


# import main functions
try:
    from growin import _core, utils, fourier_exb, growth_rate
    from growin.fourier_exb import DriftInstrument
    from growin._core import get_growth_rates_survey, get_growth
except ImportError as errstr:
    logging.exception('problem importing growin: ' + str(errstr))
