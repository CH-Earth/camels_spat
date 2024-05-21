# Add time_bnds to forcing

import glob
import sys
import pandas as pd
from pathlib import Path
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

## Config handling
# Specify where the config file can be found
config_file = '../0_config/config.txt'

# Get the required info from the config file
data_path     = cs.read_from_config(config_file,'data_path')

# CAMELS-spat metadata
cs_meta_path  = cs.read_from_config(config_file,'cs_basin_path')
cs_meta_name  = cs.read_from_config(config_file,'cs_meta_name')
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
# Get the array index from the command-line argument
if len(sys.argv) != 2:
    print("Usage: python 5_add_time_bounds.py <array_index>")
    sys.exit(1)
else:
    ix = int(sys.argv[1])

# Get the appropriate row
row = cs_meta.iloc[ix]

# Check if we need to run downloads for this station at all
missing = cs.flow_obs_unavailable(cs_unusable, row.Country, row.Station_id)
if 'iv' in missing and 'dv' in missing: 
    sys.exit(0) # graceful exit, keeps slurm logs clean

# Get forcing paths
basin_id, _, _, _, _ = cs.prepare_delineation_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder)
raw_fold, lump_fold, dist_fold = cs.prepare_forcing_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder) # Returns folders only, not file names
print('--- Now running basin {}. {}'.format(ix, basin_id))

# Find the files
raw_files = sorted(glob.glob(str(raw_fold/'*.nc'))) # list
lump_files = sorted(glob.glob(str(lump_fold/'*.nc'))) # list
dist_files = sorted(glob.glob(str(dist_fold/'*.nc'))) # list
all_files = raw_files + lump_files + dist_files

# Get LST for this station
LST = row['dv_flow_obs_timezone']
if LST == 'NST':
    print(f'NST found for {row.Station_id}, switching to AST')
    LST = 'AST' # set Newfoundland Standard Time (UTC-3h30) to Atlantic Standard Time (UTC-4h), because we only have forcing at whole hours

# Open files, add time_bnds, and close
for file in all_files:
    
    # ERA5_2023-01-01_invariants.nc' doesn't have a time variable and 
    #  thus doesn't need a time_bnds addition
    if 'invariant' in file:
        continue 

    # Add time_bnds specific to how each data set treats timestamps
    if 'era5' in file.lower():
        cs.add_time_bnds(file,'era5',LST)
    elif 'earth' in file.lower():
        cs.add_time_bnds(file,'em-earth',LST)

print(f'Completed {basin_id}')