#!/usr/bin/env python
"""
regrid.py
usage: <python> <regrid.py> <configuration.cfg>

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

# read in met location from config file
met_loc = config_dict['ECFLOW']['Met_Loc']
orig_met = config_dict['ECFLOW']['Orig_Met']
regridded_met = config_dict['ECFLOW']['Regridded_Met']

# read in grid_file from config file
grid_file = config_dict['SUBDAILY']['GridFile']

# in file
orig_file = os.path.join(met_loc, orig_met)
# out file
regrid_file = os.path.join(met_loc, regridded_met)

# remove previous days file, cdo doesn't overwrite
if os.path.isfile(regrid_file):
    os.remove(regrid_file)

cdo.remapcon(grid_file, input=orig_file, output=regrid_file)
