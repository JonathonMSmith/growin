#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2019, Jonathon Smith
# Full license can be found in License.md
# -----------------------------------------------------------------------------

from os import path, mkdir, makedirs
import sys
from setuptools import setup, find_packages


# Define a read function for using README for long_description
def read(fname):
    return open(path.join(path.dirname(__file__), fname)).read()


# generate path for fortran model files
here = path.abspath(path.dirname(__file__))
test_data_path = path.join(here, 'growin', 'tests', 'test_data')
# get environment name to create virtual environment specific archives

env_name = path.split(sys.prefix)[-1]
file_path = path.join(path.expanduser('~'), '.growin', env_name)

if not path.isdir(file_path):
    makedirs(file_path)
    print(''.join(('Created .growin directory in user home directory to',
                   'store settings.')))

with open(path.join(file_path, 'test_data_path.txt'), 'w+') as f:
    f.write(test_data_path)

setup()
