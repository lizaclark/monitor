#!/usr/bin/env python
"""
get_metdata.py
usage: <python> <get_metdata.py> <configuration.cfg>
This script downloads historic meteorological data through xarray from
http://thredds.northwestknowledge.net:8080/thredds/reacch_climate_MET_catalog.html
delivered through OPeNDAP. Because attributes are lost during download,
they are added back in.
The script also reorders dimensions of the downloaded files, regrids the data
to a 1/16th degree grid and combines all years for each parameter.
"""
import xarray as xr
import os
import argparse
import calendar
from collections import OrderedDict
from datetime import datetime, timedelta
from tonic.io import read_config
from monitor import model_tools

import numpy as np
from nco import Nco
nco = Nco()
from cdo import Cdo
cdo = Cdo()

######### ----------------------------------------###########

# read in configuration file
parser = argparse.ArgumentParser(description='Download met data')
parser.add_argument('config_file', metavar='config_file',
                    help='configuration file')
args = parser.parse_args()
config_dict = read_config(args.config_file)

downloads_loc = config_dict['METEOROLOGICAL']['Downloads_Loc']
reorder_loc = config_dict['METEOROLOGICAL']['Reorder_Loc']
regrid_loc = config_dict['METEOROLOGICAL']['Regrid_Loc']
grid_file = config_dict['METEOROLOGICAL']['Grid_File']

start_year = config_dict['METEOROLOGICAL']['Start_Year']
end_year = config_dict['METEOROLOGICAL']['End_Year']

# set a variable name for the number of lats and lons
num_lat = 584
num_lon = 1385

# define variable names used when filling threads URL
varnames = [('pr', 'precipitation_amount'), ('tmmn', 'air_temperature'),
            ('tmmx', 'air_temperature'), ('vs', 'wind_speed'),
            ('srad', 'surface_downwelling_shortwave_flux_in_air'), ('sph', 'specific_humidity')]
"""
# create attribute dictionaries

# global attributes
datestring = datetime.now()
today_date = datestring.strftime('%d %B %Y')

globe_attrs = OrderedDict()
globe_attrs[
    'author'] = "John Abatzoglou - University of Idaho, jabatzoglou@uidaho.edu"
globe_attrs['date'] = today_date
globe_attrs[
    'note1'] = "The projection information for this file is: GCS WGS 1984."
globe_attrs['note2'] = ("Citation: Abatzoglou, J.T., 2013, Development of gridded surface " +
                        "meteorological data for ecological applications and modeling, " +
                        "International Journal of Climatology, DOI: 10.1002/joc.3413")
globe_attrs['last_permanent_slice'] = "50"
globe_attrs['note3'] = ("Data in slices after last_permanent_slice (1-based) are " +
                        "considered provisional and subject to change with subsequent updates")

# latitude attributes
lat_attrs = OrderedDict()
lat_attrs['units'] = "degrees_north"
lat_attrs['description'] = "latitude"

# longitude attributes
lon_attrs = OrderedDict()
lon_attrs['units'] = "degrees_east"
lon_attrs['description'] = "longitude"

# time attributes
day_attrs = OrderedDict()
day_attrs['units'] = "days since 1900-01-01 00:00:00"
day_attrs['calendar'] = "gregorian"
day_attrs['description'] = "days since 1900-01-01"

esri_str = ("GEOGCS[\\\"GCS_WGS_1984\\\",DATUM" +
            "[\\\"D_WGS_1984\\\",SPHEROID\\\"WGS_1984\\\"," +
            "6378137.0,298.257223563]],PRIMEM[\\\"Greenwich\\\",0.0]," +
            "UNIT[\\\"Degree\\\",0.0174532925199433]]")

# parameter attributes
# precipitation
pr_attrs = OrderedDict()
pr_attrs['units'] = "mm"
pr_attrs['description'] = "Daily Accumulation Precipitation"
pr_attrs['_FillValue'] = -32767.
pr_attrs['esri_pe_string'] = esri_str
pr_attrs['coordinates'] = "lon lat"
pr_attrs['cell_methods'] = "time: sum(intervals: 24 hours)"
pr_attrs['missing_value'] = -32767.

# minimum temperature
tmmn_attrs = OrderedDict()
tmmn_attrs['units'] = "K"
tmmn_attrs['description'] = "Daily Minimum Temperature"
tmmn_attrs['_FillValue'] = -32767.
tmmn_attrs['esri_pe_string'] = esri_str
tmmn_attrs['coordinates'] = "lon lat"
tmmn_attrs['cell_methods'] = "time: sum(interval: 24 hours)"
tmmn_attrs['height'] = "2 m"
tmmn_attrs['missing_value'] = -32767.

# maximum temperature
tmmx_attrs = OrderedDict()
tmmx_attrs['units'] = "K"
tmmx_attrs['description'] = "Daily Maximum Temperature"
tmmx_attrs['_FillValue'] = -32767.
tmmx_attrs['esri_pe_string'] = esri_str
tmmx_attrs['coordinates'] = "lon lat"
tmmx_attrs['cell_methods'] = "time: sum(interval: 24 hours)"
tmmx_attrs['height'] = "2 m"
tmmx_attrs['missing_value'] = -32767.

# wind speed
vs_attrs = OrderedDict()
vs_attrs['units'] = "m/s"
vs_attrs['description'] = "Daily Mean Wind Speed"
vs_attrs['_FillValue'] = -32767.
vs_attrs['esri_pe_string'] = esri_str
vs_attrs['coordinates'] = "lon lat"
vs_attrs['height'] = "10 m"
vs_attrs['missing_value'] = -32767.

# shortwave radiation
srad_attrs = OrderedDict()
srad_attrs['units'] = "W m-2"
srad_attrs['description'] = "Daily Mean Downward Shortwave Radiation At Surface"
srad_attrs['_FillValue'] = -32767.
srad_attrs['esri_pe_string'] = esri_str
srad_attrs['coordinates'] = "lon lat"
srad_attrs['missing_value'] = -32767.

# specific humidity
sph_attrs = OrderedDict()
sph_attrs['units'] = "kg/kg"
sph_attrs['description'] = "Daily Mean Specific Humidity"
sph_attrs['_FillValue'] = -32767.
sph_attrs['esri_pe_string'] = esri_str
sph_attrs['coordinates'] = "lon lat"
sph_attrs['height'] = "2 m"
sph_attrs['missing_value'] = -32767.

# download metdata from http://thredds.northwestknowledge.net
met_dsets = dict()
num_startofyear = 0

#for year in range(start_year, end_year+1):
for year in range(2015,2016):
    for var in varnames:
        if calendar.isleap(int(year)):
            num_endofyear = 365
        else:
            num_endofyear = 364

        url = ("http://thredds.northwestknowledge.net:8080" +
               "/thredds/dodsC/MET/%s/%s_%s.nc?lon[0:1:%s]," % (var[0], var[0], year, num_lon) +
               "lat[0:1:%s],day[%s:1:%s]," % (num_lat, num_startofyear, num_endofyear) +
               "%s[%s:1:%s]" % (var[1], num_startofyear, num_endofyear) +
               "[0:1:%s][0:1:%s]" % (num_lon, num_lat))

        # download data and add general attributes
        ds = xr.open_dataset(url)
        ds.lat.attrs = lat_attrs
        ds.lon.attrs = lon_attrs
        ds.day.attrs = day_attrs
        ds.attrs = globe_attrs

        # place data in dict
        met_dsets[os.path.join('%s_%s' % (var[0], year))] = ds

# add variable specific attributes and save as netcdf
#for year in range(start_year, end_year):
for year in range(2015, 2016):
    met_dsets['pr_%s' % (year)].precipitation_amount.attrs = pr_attrs
    met_dsets['pr_%s' % (year)].to_netcdf(os.path.join(downloads_loc, 'pr_%s.nc' % (year)),
                          mode='w', format='NETCDF4')

    met_dsets['tmmn_%s' % (year)].air_temperature.attrs = tmmn_attrs
    met_dsets['tmmn_%s' % (year)].to_netcdf(os.path.join(downloads_loc, 'tmmn_%s.nc' % (year)),
                            mode='w', format='NETCDF4')

    met_dsets['tmmx_%s' % (year)].air_temperature.attrs = tmmx_attrs
    met_dsets['tmmx_%s' % (year)].to_netcdf(os.path.join(downloads_loc, 'tmmx_%s.nc' % (year)),
                            mode='w', format='NETCDF4')

    met_dsets['vs_%s' % (year)].wind_speed.attrs = vs_attrs
    met_dsets['vs_%s' % (year)].to_netcdf(os.path.join(downloads_loc, 'vs_%s.nc' % (year)),
                          mode='w', format='NETCDF4')

    met_dsets['srad_%s' %
        (year)].surface_downwelling_shortwave_flux_in_air.attrs = srad_attrs
    met_dsets['srad_%s' % (year)].to_netcdf(os.path.join(downloads_loc, 'srad_%s.nc' % (year)),
                            mode='w', format='NETCDF4')

    met_dsets['sph_%s' % (year)].specific_humidity.attrs = sph_attrs
    met_dsets['sph_%s' % (year)].to_netcdf(os.path.join(downloads_loc, 'sph_%s.nc' % (year)),
                           mode='w', format='NETCDF4')

# ncpdq to reorder dimensions
for year in range(start_year, end_year + 1):
    for var in varnames:
        # in file
        nc_file = os.path.join(downloads_loc, '%s_%s.nc' % (var[0], year))
        # out file
        reorder_file = os.path.join(reorder_loc, '%s_%s.reorder.nc' % (var[0], year))

        # remove a previous file, ncpdq doesn't overwrite
        if os.path.isfile(reorder_file):
            os.remove(reorder_file)

        nco.ncpdq(input=nc_file, output=reorder_file, arrange='day,lat,lon')

# cdo remapcon to change grid and domain

for year in range(start_year, end_year + 1):
    for var in varnames:
        # in file
        reorder_file = os.path.join(
            reorder_loc, '%s_%s.reorder.nc' % (var[0], year))
        # out file
        regrid_file = os.path.join(
            regrid_loc, '%s_%s.regrid.nc' % (var[0], year))

        # remove previous days file, cdo doesn't overwrite
        if os.path.isfile(regrid_file):
            os.remove(regrid_file)

        cdo.remapcon(grid_file, input=reorder_file, output=regrid_file)
"""
# concatenate files for each variable
for var in varnames:
    list_of_filenames = []
    for year in range(start_year, end_year + 1):
        list_of_filenames.append(os.path.join(
            regrid_loc, '%s_%s.regrid.nc' % (var[0], year)))
    out_path = os.path.join(regrid_loc, '%s_%s_%s.nc' %
                            (var[0], start_year, end_year))
    # remove a previous file, cdo doesn't overwrite
    if os.path.isfile(out_path):
        os.remove(out_path)
    # concatenate all years
    cdo.mergetime(input=" ".join(list_of_filenames),
                            output=out_path)

