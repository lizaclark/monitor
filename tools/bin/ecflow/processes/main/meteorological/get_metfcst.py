#!/usr/bin/env python
"""
get_metfcst.py
usage: <python> <get_metfcst.py> <configuration.cfg>
This script downloads downscaled NMME meteorological forecast data through
xarray from
http://proxy.nkn.uidaho.edu/thredds/catalog/
    NWCSC_INTEGRATED_SCENARIOS_ALL_CLIMATE/bcsd-nmme/dailyForecasts/
    catalog.html
delivered through OPeNDAP. Because attributes are lost during download,
they are added back in. To start, we just download multi-model ensemble mean.
"""
import xarray as xr
import os
import argparse
import calendar
from collections import OrderedDict
from datetime import datetime, timedelta
import cf_units

from tonic.io import read_config
from monitor import model_tools

# read in configuration file
parser = argparse.ArgumentParser(description='Download met forecast data')
parser.add_argument('config_file', metavar='config_file',
                    help='configuration file')
args = parser.parse_args()
config_dict = read_config(args.config_file)

# read in meteorological data location
met_fcst_loc = config_dict['FORECAST']['metfcst_Loc']
old_config_file = config_dict['ECFLOW']['old_Config']
new_config_file = config_dict['ECFLOW']['new_Config']

# define variable names used when filling threads URL
# an abbreviation and a full name is needed
varnames = ['was', 'dps', 'tasmean', 'tasmax', 'tasmin', 'rsds',
            'pr', 'pet']
# define model names for file name
modelnames = ['GFDL']
# ['NCAR', 'NASA', 'GFDL', 'GFDL-FLOR', 'ENSMEAN', 'CMC1', 'CMC2',
#              'CFSv2']
new_units = {'was': 'm s-1', 'dps': 'degC', 'tasmean': 'degC',
             'tasmax': 'degC', 'tasmin': 'degC', 'rsds': 'W m-2',
             'pr': 'mm', 'pet': 'mm'}

# download metdata from http://thredds.northwestknowledge.net
met_dsets = dict()
for model in modelnames:
    for var in varnames:
        url = ('http://tds-proxy.nkn.uidaho.edu/thredds/dodsC/' +
               'NWCSC_INTEGRATED_SCENARIOS_ALL_CLIMATE/bcsd-nmme/' +
               'dailyForecasts/' +
               'bcsd_nmme_metdata_%s_forecast_%s_daily.nc' % (
                   model, var))
        print('Reading {0}'.format(url))
        ds = xr.open_dataset(url)
        if ds[var].attrs['units'] == 'F':
            ds[var].attrs['units'] = 'degF'
        units_in = cf_units.Unit(ds[var].attrs['units'])
        units_out = cf_units.Unit(new_units[var])
        # Perform units conversion
        units_in.convert(ds[var].values[:], units_out, inplace=True)
        ds[var].attrs['units'] = new_units[var]
        outfile = os.path.join(met_fcst_loc, '%s_%s.nc' % (model, var))
        print('writing {0}'.format(outfile))
        ds.to_netcdf(outfile, mode='w', format='NETCDF4')
