#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2019, Jonathon Smith
# Full license can be found in License.md
# -----------------------------------------------------------------------------

from os import path, mkdir, makedirs
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

setup(name='growin',
      version='0.1a2',
      url='gitlab.com/jklenzing/growin',
      author='Jeff Klenzing',
      author_email='jeffrey.klenzing@nasa.gov',
      description='Generate, read, and plot sami2 model runs',
      long_description=read('README.md'),
      packages=find_packages(),
      classifiers=[
          "Development Status :: 3 - Alpha",
          "Topic :: Scientific/Engineering :: Physics",
          "Intended Audience :: Science/Research",
          "License :: BSD",
          "Natural Language :: English",
          "Programming Language :: Python :: 2.7",
          "Programming Language :: Python :: 3.3",
          "Programming Language :: Python :: 3.4",
          "Programming Language :: Python :: 3.5",
          "Programming Language :: Python :: 3.6",
      ],
      include_package_data=True,
      zip_safe=False,
      test_suite='setup.growin_test_suite',
      )
