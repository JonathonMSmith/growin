<!--
                                                 /$$          
                                                |__/          
      /$$$$$$   /$$$$$$   /$$$$$$  /$$  /$$  /$$ /$$ /$$$$$$$ 
     /$$__  $$ /$$__  $$ /$$__  $$| $$ | $$ | $$| $$| $$__  $$
    | $$  \ $$| $$  \__/| $$  \ $$| $$ | $$ | $$| $$| $$  \ $$
    | $$  | $$| $$      | $$  | $$| $$ | $$ | $$| $$| $$  | $$
    |  $$$$$$$| $$      |  $$$$$$/|  $$$$$/$$$$/| $$| $$  | $$
     \____  $$|__/       \______/  \_____/\___/ |__/|__/  |__/
     /$$  \ $$                                                
    |  $$$$$$/                                                
     \______/                                                 
 -->   
<div align="center">
        <img height="0" width="0px">
        <img width="50%" src="/banner.png" alt="growin" title="growin"</img>
</div>


[![DOI](https://zenodo.org/badge/174235815.svg)](https://zenodo.org/badge/latestdoi/174235815)
[![Coverage Status](https://coveralls.io/repos/github/JonathonMSmith/growin/badge.svg?branch=master)](https://coveralls.io/github/JonathonMSmith/growin?branch=master)
[![Maintainability](https://api.codeclimate.com/v1/badges/c9a94135bded3475dea7/maintainability)](https://codeclimate.com/github/JonathonMSmith/growin/maintainability)
[![Documentation Status](https://readthedocs.org/projects/growin/badge/?version=latest)](https://growin.readthedocs.io/en/latest/?badge=latest)

#Overview

The flux-tube-integrated linear Rayleigh-Taylor instability (grow)th rate of 
(i)o(n)ospheric irregularities is often used in the study of equatorial spread 
F(ESF). The Rayleigh-Taylor instability is often accepted as the mechanism
responsible for ESF, therefore a computation of its growth rate is valuable. 
This package includes many common approximations used in the computation of the 
Rayleigh-Taylor growth rate.

#Installation

Python version 3.6 or greater is recommended

growin has many dependencies:
future
apexpy
netCDF4
numpy
matplotlib
pysat
sami2py
pyglow

and for best results try:
```
pip install future
pip install apexpy
pip install netCDF4
pip install numpy
pip install matplotlib
pip install pysat
pip install igrf12

git clone https://github.com/sami2py/sami2py.git
cd sami2py/python setup.py install
make -C sami2py/fortran compile

git clone https://github.com/JonathonMSmith/pyglow.git
cd pyglow
git checkout fejer_output
python setup.py install
```
Finally, clone and install growin

```
git clone https://github.com/JonathonMSmith/growin.git
cd growin/python setup.py install
```

Example
-------

In an interactive python shell, run:
```
import growin
```
growin and sami2py will remind you to set the top level directory that will hold the model output.
```
growin.utils.set_archive_dir(path=path)
sami2py.utils.set_archive_dir(path=path)
```
sami2py and growin will raise an error if this is not done before trying to run the model.

Now it's time to calculate some growth rates with a trivial example
```
growin.get_growth(tag='babys_first_growth_rate', day=42, year=2010, lon=81) 
```

this might take a while...
