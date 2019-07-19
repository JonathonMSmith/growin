#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2019, Jonathon Smith
# Full license can be found in License.md
# -----------------------------------------------------------------------------
"""
irfl
-----------
"""

import os
import logging

__version__ = str('0.1a1')

# get home directory
home_dir = os.path.expanduser('~')
# get environment name
if 'CONDA_DEFAULT_ENV' in os.environ:
    env_name = os.environ['CONDA_DEFAULT_ENV']
elif 'VIRTUAL_ENV' in os.environ:
    env_name = os.environ['VIRTUAL_ENV']
# set irfl directory path in home directory
irfl_dir = os.path.join(home_dir, '.irfl', env_name)
# make sure a irfl directory for model output exists
if not os.path.isdir(irfl_dir):
    # create directory
    os.mkdir(irfl_dir)
    print(''.join(('Created .irfl directory in user home directory to',
                   'store settings.')))


archive_path = os.path.join(irfl_dir, 'archive_path.txt')
if os.path.isfile(archive_path):
    # load up stored data path
    with open(archive_path, 'r') as f:
        archive_dir = f.readline()
else:
    # create file
    with open(archive_path, 'w+') as f:
        f.write('')
    archive_dir = ''
    print(''.join(('Run irfl.utils.set_archive_dir to set the path to',
                   ' top-level directory for model outputs.')))

# load test_data directory
with open(os.path.join(irfl_dir, 'test_data_path.txt'), 'r') as f:
    test_data_dir = f.readline()


# import main functions
try:
    from irfl import _core, utils, fourier_exb, growth_rate
    from irfl.fourier_exb import DriftInstrument
    from irfl._core import get_growth_rates_survey
except ImportError as errstr:
    logging.exception('problem importing irfl: ' + str(errstr))
