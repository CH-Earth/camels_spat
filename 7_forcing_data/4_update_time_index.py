# UTC to LST
# We want all the forcing to be in local standard time to match the streamflow observations.

import glob
import os
import sys
import netCDF4 as nc4
import pandas as pd
import shutil
import xarray as xr
from pathlib import Path
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

# --- Reruns 2024-05-21
# These fix various small errors discovered during data use
rerun_file = Path('/globalhome/wmk934/HPC/camels_spat/7_forcing_data/forcing_check_logs/reruns_20240516.csv')
reruns = pd.read_csv(rerun_file)
# --- Reruns 2024-05-21

## Config handling
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

## Data loading
# CAMELS-spat metadata file
cs_meta_path = Path(data_path) / cs_meta_path
cs_meta = pd.read_csv(cs_meta_path / cs_meta_name)
cs_unusable = pd.read_csv(cs_meta_path / cs_unusable_name,  dtype={'Station_id': object}) # Enforce reading IDs as string to keep leading 0's

## Processing
debug_message = f'\n!!! CHECK DEBUGGING STATUS: \n- Full run\n'

# Get the array index from the command-line argument
if len(sys.argv) != 2:
    print("Usage: python 4_update_time_index.py <array_index>")
    sys.exit(1)
else:
    ix = int(sys.argv[1])

# Define cs_meta row
row = cs_meta.iloc[ix] # needs to be between 0  and 1697

# --- Reruns 2024-05-21
this_basin = row['Country'] + '_' + row['Station_id']
if this_basin not in reruns['basin'].values:
    print(f'No reruns for basin {this_basin}. Exiting.')
    sys.exit(0) # with next station, because we have no reruns for this station. Error code 0: clean exit, no problems
else:
    print(f'Running reruns for basin {this_basin}.')
# --- Reruns 2024-05-21

# Get forcing paths
basin_id, _, _, _, _ = cs.prepare_delineation_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder)
raw_fold, lump_fold, dist_fold = cs.prepare_forcing_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder) # Returns folders only, not file names
print('--- Now running basin {}. {}'.format(ix, basin_id))

# Create a backup folder at the same level as raw_fold, lump_fold, dist_fold
backup_fold = raw_fold.parent / 'backup_before_time_update'
backup_fold.mkdir(exist_ok=True)

# Check if we need to run downloads for this station at all
missing = cs.flow_obs_unavailable(cs_unusable, row.Country, row.Station_id)
if 'iv' in missing and 'dv' in missing: 
    sys.exit(0) # graceful exit because we don't need to process this basin and don;'t want to pollute the SLURM summary with unhelpful job failed messages
    
# Find the files
raw_files = sorted(glob.glob(str(raw_fold/'*.nc'))) # list
lump_files = sorted(glob.glob(str(lump_fold/'*.nc'))) # list
dist_files = sorted(glob.glob(str(dist_fold/'*.nc'))) # list
all_files = raw_files + lump_files + dist_files

# Rerun 20240521 - ONLY do the lump and dist EM_Earth files
# We need to do this because during remapping I used the wrong EM-Earth filename
#  and then proceeded with the rest of the workflow (3b_add_wind_dir, 4_update_time_index,
#  5_add_time_bounds) without noticing the error. The lump and dist files are the only ones
#  that need to be updated, because the raw files were processed correctly.
all_files = lump_files + dist_files
all_files = [file for file in all_files if 'EM_Earth' in file]
# Rerun 20240521

# Find LST
# We can simply use dv_flow_obs_timezone here because we already know the USA gauges
#  show consistent LSTs for IV and DV observations, and for CAN we only have DV LST
#  anyway.
lst = row['dv_flow_obs_timezone']
if lst == 'NST':
    print(f'NST found for {row.Station_id}, switching to AST')
    lst = 'AST' # set Newfoundland Standard Time (UTC-3h30) to Atlantic Standard Time (UTC-4h), because we only have forcing at whole hours
utc = cs.tz_abbreviation_to_utc(lst) # e.g. 'UTC-04'
offset = cs.relative_utc_to_float_offset_in_hours(utc) # e.g. -4.0
    
# Open files with xarray, update time values with pandas and replace in file
for file in all_files:
    if 'invariant' in file:
        continue # skip the ERA5 invariants that don't have a time dimension
    print(f'Processing {file}')
    
    # create a backup copy of the file in case processing ends halfway through a file
    backup_file = file.replace('.nc','_BACKUP.nc')
    shutil.copy2(file,backup_file)

    # Rerun 20240521 - create an extra backup for the NST case
    # This seems safer because in case we need to re-do the time update
    if row['dv_flow_obs_timezone'] == 'NST':
        file_name = os.path.basename(file)
        backup_file2 = backup_fold / file_name
        shutil.copy2(file,backup_file2)
    # Rerun 20240521 - create an extra backup for the NST case
    
    # Update times
    with nc4.Dataset(file, 'a') as f:
        time_variable = f.variables['time']
        time_variable[:] = time_variable[:] + offset
    
    # Remove the backup
    os.remove(backup_file)

print(debug_message)


