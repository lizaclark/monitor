#!/usr/bin/env python
"""
create_metsim_fcst_cfg.py
usage: <python> <create_metsim_fcst_cfg.py> <configuration.cfg>
This script writes a MetSim configuration file and runs MetSim for
NMME forecast data.
"""
import os
import argparse
import subprocess
import xarray as xr
import pandas as pd

from monitor import model_tools
from tonic.io import read_config

# read in configuration file
parser = argparse.ArgumentParser(description='Create configuration file ' +
                                 'for MetSim')
parser.add_argument('config_file', metavar='config_file',
                    help='configuration file')
args = parser.parse_args()
config_dict = read_config(args.config_file)

state = config_dict['SUBDAILY']['Fcst_Met_State']
old_config_file = config_dict['FORECAST']['MetSim_Template']
new_config_file = config_dict['FORECAST']['MetSim_Cfg']
start_date = config_dict['SUBDAILY']['Subd_Met_End_Date']
out_dir = config_dict['FORECAST']['metfcst_Loc']
out_state = os.path.join(out_dir, 'fcst_state_out.nc')
modelvars = ['ENSMEAN']
for model in modelvars:
    forecast = os.path.join(config_dict['FORECAST']['metfcst_Loc'],
                            '{}.nc'.format(model))
    xds = xr.open_dataset(forecast)
    end_date = pd.to_datetime(str(xds.time.values[-1])).strftime('%Y/%m/%d')
    kwargs = {'STARTDATE': start_date.replace('-', '/'), 'ENDDATE': end_date,
              'FORCING': forecast, 'DOMAIN': config_dict['SUBDAILY']['GridFile'],
              'INSTATE': state, 'OUTDIR': out_dir, 'OUTSTATE': out_state}
    model_tools.replace_var_pythonic_config(
        old_config_file, new_config_file, header=None, **kwargs)
    subprocess.check_call(['ms', new_config_file])
