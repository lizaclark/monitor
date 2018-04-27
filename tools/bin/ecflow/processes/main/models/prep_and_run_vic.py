#!/usr/bin/env python
"""
create_global_file.py
usage: <python> <create_global_file.py> <configuration.cfg>

This script creates a global file from the template with the correct
 model start and end dates.
"""
import argparse
import subprocess
from tonic.io import read_config
from monitor import model_tools
import os
from dateutil.parser import parse
######### ----------------------------------------###########

# read in configuration file
parser = argparse.ArgumentParser(description='Reorder dimensions')
parser.add_argument('config_file', metavar='config_file',
                    help='configuration file')
parser.add_argument('time_horizon_type', help='MONITOR, MED_FCST, or ' +
                    'SEAS_FCST. Should correspond to section header in ' +
                    'config_file')
args = parser.parse_args()
config_dict = read_config(args.config_file)
section = args.time_horizon_type

# read in the desired variables from the config file
global_template = config_dict['DOMAIN']['GlobalFileTemplate']
global_file_path = config_dict[section]['GlobalFilePath']
state_path = config_dict[section]['StatePath']
in_state_path = config_dict[section]['InStatePath']
result_path = config_dict[section]['OutputDirRoot']
forc_path = config_dict[section]['Subd_Out_Dir']
# get important dates
vic_start_date = config_dict[section]['Start_Date']
vic_end_date = config_dict[section]['End_Date']
vic_save_state = config_dict[section]['vic_save_state']

forcing_prefix = os.path.join(
    forc_path, 'forcing_{0}-{1}.'.format(vic_start_date.replace('-', ''),
                                   vic_end_date.replace('-', ''))) 

# parse out year, month, and day from the model dates, which have
# the form YYYY-MM-DD
start = parse(vic_start_date)
end = parse(vic_end_date)
save_state = parse(vic_save_state)

# generate the path to the initial state file
in_state = os.path.join(in_state_path, 'state.%s%s%s_00000.nc' % (
    vic_start_date[:4], vic_start_date[5:7], vic_start_date[8:10]))

kwargs = {
    'Start_Year': start.year,
    'Start_Month': start.month,
    'Start_Day': start.day,
    'End_Year': end.year,
    'End_Month': end.month,
    'End_Day': end.day,
    'State_Year': save_state.year,
    'State_Month': save_state.month,
    'State_Day': save_state.day,
    'In_State': in_state,
    'State_Path': state_path,
    'Result_Path': result_path,
    'Forcing_Prefix': forcing_prefix}


model_tools.replace_var_pythonic_config(
    global_template, global_file_path, header=None, **kwargs)

subprocess.run(['mpi.hydra', '-np', '15', config_dict['ECFLOW']['Executable'], '-g', global_file_path])
