# Certain files have very small mistmaches (e.g. 7.6293945e-06) in the lat/lon coordinates
# This is rare (typically 1 file only) and doesn't occur for every lat/lon value either, 
# so we assume this is a processing error in the source data, and that we can therefore 
# get away with overwriting those values.

import os
import pandas as pd
from pathlib import Path
import shutil
import sys
import xarray as xr
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

# --- Command line arguments
if len(sys.argv) != 2:
    print("Usage: python 6e_merge_daymet_vars_to_single_year <array_index>")
    sys.exit(1)
else:
    ix = int(sys.argv[1]) # Basin index in cs_meta

# --- Config handling
# Specify where the config file can be found
config_file = '../0_config/config.txt'

# Get the required info from the config file
data_path = cs.read_from_config(config_file,'data_path')

# CAMELS-spat metadata
cs_meta_path = cs.read_from_config(config_file,'cs_basin_path')
cs_meta_name = cs.read_from_config(config_file,'cs_meta_name')
cs_unusable_name = cs.read_from_config(config_file,'cs_unusable_name')

# Basin folder
cs_basin_folder = cs.read_from_config(config_file, 'cs_basin_path')
basins_path = Path(data_path) / cs_basin_folder

# --- Data loading
# CAMELS-spat metadata file
cs_meta_path = Path(data_path) / cs_meta_path
cs_meta = pd.read_csv(cs_meta_path / cs_meta_name)

# Open list of unusable stations; Enforce reading IDs as string to keep leading 0's
cs_unusable = pd.read_csv(cs_meta_path / cs_unusable_name, dtype={'Station_id': object})

# --- Processing
# Define cs_meta row
row = cs_meta.iloc[ix] # needs to be between 0  and 1697
basin_id = row['Country'] + '_' + row['Station_id']

# Check if we need to run for this station at all
missing = cs.flow_obs_unavailable(cs_unusable, row.Country, row.Station_id)
if 'iv' in missing and 'dv' in missing: 
    print(f'No flow observations for basin {basin_id}. Exiting.')
    sys.exit(0) # with next station, because we have no observations at all for this station. Error code 0: clean exit, no problems
print('--- Now running basin {}. {}'.format(ix, basin_id))

# Find the Daymet files
raw_fold, _, _ = cs.prepare_forcing_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder) # Returns folders only, not file names
daymet_files = list(raw_fold.glob('daymet_*.nc')) 
daymet_files.sort()

# Open the first and extract lat and lon values
data = xr.open_dataset(daymet_files[0])
src_lat = data['lat'].values
src_lon = data['lon'].values
src_x = data['x'].values
src_y = data['y'].values
data.close()

# Loop over the rest
for file in daymet_files[1:]:

    # Reset the flags
    flag_lat = False
    flag_lon = False
    flag_x = False
    flag_y = False

    # Check if we need to overwrite anything, and report what the differences are
    data = xr.open_dataset(file)

    lat_diff = abs(data['lat'].values - src_lat)
    if (lat_diff > 0).any():
        flag_lat = True
        num_dif = (lat_diff > 0).sum()
        siz_dif = lat_diff[lat_diff > 0].mean()
        print(f'\nLatitude differences found in {file}')
        print(f' - Number of different values:  {num_dif}')
        print(f' - Average size of differences: {siz_dif}')

    lon_diff = abs(data['lon'].values - src_lon)
    if (lon_diff > 0).any():
        flag_lon = True
        num_dif = (lon_diff > 0).sum()
        siz_dif = lon_diff[lon_diff > 0].mean()
        print(f'\nLongitude differences found in {file}')
        print(f' - Number of different values:  {num_dif}')
        print(f' - Average size of differences: {siz_dif}')

    x_diff = abs(data['x'].values - src_x)
    if (x_diff > 0).any():
        flag_x = True
        num_dif = (x_diff > 0).sum()
        siz_dif = x_diff[x_diff > 0].mean()
        print(f'\nX differences found in {file}')
        print(f' - Number of different values:  {num_dif}')
        print(f' - Average size of differences: {siz_dif}')

    y_diff = abs(data['y'].values - src_y)
    if (y_diff > 0).any():
        flag_y = True
        num_dif = (y_diff > 0).sum()
        siz_dif = y_diff[y_diff > 0].mean()
        print(f'\nY differences found in {file}')
        print(f' - Number of different values:  {num_dif}')
        print(f' - Average size of differences: {siz_dif}')

    data.close()

    # If flagged, create a temporary copy of the file and use that to overwrite data
    if flag_lat or flag_lon or flag_x or flag_y:
        temp_file = file.with_name(file.name + '.tmp')
        shutil.copy(file, temp_file)
        data = xr.open_dataset(temp_file)
        if flag_lat:
            data['lat'].values = src_lat
        if flag_lon:
            data['lon'].values = src_lon
        if flag_x:
            data['x'].values = src_x
        if flag_y:
            data['y'].values = src_y
        data.to_netcdf(file)
        data.close()
        os.remove(temp_file)

