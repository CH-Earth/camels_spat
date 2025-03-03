# RDRS UTC to LST
# We want all the forcing to be in local standard time to match the streamflow observations.

import glob
import os
import sys
import netCDF4 as nc4
import pandas as pd
import shutil
import tarfile
import time
import xarray as xr
from pathlib import Path
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

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
    print("Usage: python 7g_update_time_index.py <array_index>")
    sys.exit(1)
else:
    ix = int(sys.argv[1])

# Define cs_meta row
row = cs_meta.iloc[ix] # needs to be between 0  and 1697

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

# Clean-up if needed
if (raw_fold/'rdrs_month').exists():
    # Move the RDRS raw files into the main raw folder
    raw_files = sorted(glob.glob(str(raw_fold/'rdrs_month'/'RDRS*.nc'))) # list
    for raw_file in raw_files:
        shutil.copy2(raw_file,raw_fold/Path(raw_file).name)

    # Now tar the rdrs_month folder to save space but keep the backup
    tar_month = raw_fold / 'rdrs_month_tar' 
    tar_month.mkdir(exist_ok=True, parents=True)
    tar_file = tar_month / f'{basin_id}_rdrs_month.tar.gz'
    with tarfile.open(tar_file, "w:gz") as tar:
        for file in raw_files:
            tar.add(file, arcname=Path(file).name)

    # Remove the rdrs_month folder
    shutil.rmtree(raw_fold/'rdrs_month')

# Find the files now we've moved the whole thing
raw_files = sorted(glob.glob(str(raw_fold/'RDRS*.nc'))) # list
lump_files = sorted(glob.glob(str(lump_fold/'RDRS*.nc'))) # list
dist_files = sorted(glob.glob(str(dist_fold/'RDRS*.nc'))) # list
all_files = raw_files + lump_files + dist_files

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
   
    # create a backup copy of the file in case processing ends halfway through a file
    backup_file = file.replace('.nc','_BACKUP.nc')
    shutil.copy2(file,backup_file)
    
    # Update times
    with nc4.Dataset(file, 'a') as f:
        time_variable = f.variables['time']
        time_variable[:] = time_variable[:] + offset

        # Track what we did
        history = f' On {time.ctime(time.time())}: updated time index from UTC to {lst}.'
        if 'History' in f.ncattrs():
            current_history = f.getncattr('History')
            new_history = f'{current_history}{history}'
            f.setncattr('History', new_history)
        elif 'history' in f.ncattrs():
            current_history = f.getncattr('history')
            new_history = f'{current_history}{history}'
            f.setncattr('history', new_history)
    
    # Remove the backup
    os.remove(backup_file)

print(debug_message)


