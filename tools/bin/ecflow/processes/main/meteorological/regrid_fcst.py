#!/usr/bin/env python
"""
regrid.py
usage: <python> <regrid_fcst.py> <configuration.cfg>

This script changes the grid and domain in accordance with the grid_file, using
remapcon.
Remapcon was selected to ensure that no precipitation was lost in the
regridding process.
"""
import os
import argparse
from cdo import Cdo

from tonic.io import read_config

cdo = Cdo()
# read in configuration file
parser = argparse.ArgumentParser(description='Regrid met data')
parser.add_argument('config_file', metavar='config_file',
                    help='configuration file')
args = parser.parse_args()
config_dict = read_config(args.config_file)
# set location of cdo binary
cdo.setCdo(config_dict['ECFLOW']['Cdo_Bin'])

# read in met location from config file
fcst_loc = config_dict['FORECAST']['metfcst_Loc']

# read in grid_file from config file
grid_file = config_dict['SUBDAILY']['GridFile']

# define model names for file name
modelnames = ['NCAR', 'NASA', 'GFDL', 'GFDL-FLOR', 'ENSMEAN', 'CMC1', 'CMC2',
              'CFSv2']

for model in modelnames:
        # in file
        orig_file = os.path.join(fcst_loc, '%s.nc' % (model))
        # out file
        regrid_file = os.path.join(fcst_loc, '%s.regrid.nc' % (model))

        # remove previous days file, cdo doesn't overwrite
        if os.path.isfile(regrid_file):
            os.remove(regrid_file)

        cdo.remapcon(grid_file, input=orig_file, output=regrid_file)
