# PECANS
Python Editable Chemical Atmospheric Numeric Solver

[![DOI](https://zenodo.org/badge/94652282.svg)](https://zenodo.org/badge/latestdoi/94652282)

## Getting started

To install, simply download or clone this repo to your computer and install the requirements
described in `requirements.txt` using `pip install -r requirements.txt` (called in the top
directory of this repo, where `requirements.txt` is). I recommend installing in a virtual
environment or Conda environment to avoid dependency conflicts with your base Python installation.

PECANS can also be installed as a package by running `python setup.py install` or `python setup.py develop`
in the top directory. This allows you to, for example, import PECANS into a Jupyter notebook.
The difference between the `install` and `develop` commands is that the latter will link to 
the copy of the repo downloaded, so that if you pull updates or otherwise change the code, those
changes will automatically be available to the package without needed to rerun the setup command.

## Documentation

Full documentation is available at https://pecans.readthedocs.io/en/stable/.

## Citing 

The DOI in the badge at the top of this document will reference the most recent version of PECANS.
Specific versions available are:

* [v0.1.0](https://doi.org/10.5281/zenodo.3386652)
