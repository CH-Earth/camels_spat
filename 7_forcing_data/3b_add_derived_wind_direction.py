# Add wind direction
import glob
import sys
import netCDF4 as nc4
import pandas as pd
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

# Open list of unusable stations; Enforce reading IDs as string to keep leading 0's
cs_unusable = pd.read_csv(cs_meta_path / cs_unusable_name, dtype={'Station_id': object})

## Processing
debug_message = f'\n!!! CHECK DEBUGGING STATUS: \n- Full run in progress\n'

# Get the array index from the command-line argument
if len(sys.argv) != 2:
    print("Usage: python 3b_add_derived_wind_direction.py <array_index>")
    sys.exit(1)
else:
    ix = int(sys.argv[1])

# Define cs_meta row
row = cs_meta.iloc[ix] # needs to be between 0  and 1697

print(debug_message)
    
# Get shapefile path to determine download coordinates, and forcing destination path
basin_id, shp_lump_path, shp_dist_path, _, _ = cs.prepare_delineation_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder)
raw_fold, lump_fold, dist_fold = cs.prepare_forcing_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder) # Returns folders only, not file names
print('--- Now running basin {}. {}'.format(ix, basin_id))
    
# Check if we need to run downloads for this station at all
missing = cs.flow_obs_unavailable(cs_unusable, row.Country, row.Station_id)
if 'iv' in missing and 'dv' in missing: 
    sys.exit(0) # gracefully exit run, because we have no observations at all for this station

# Find the files
era5_merged_files = sorted(glob.glob(str(raw_fold/'ERA5_[0-9][0-9][0-9][0-9]-[0-9][0-9].nc'))) # list
era5_lump_files = sorted(glob.glob(str(lump_fold/'ERA5_lumped_*.nc'))) # list
era5_dist_files = sorted(glob.glob(str(dist_fold/'ERA5_dist_*.nc'))) # list
era5_all_files = era5_merged_files + era5_lump_files + era5_dist_files

# Loop over the files and add new variables
for file in era5_all_files:
    print(f'Processing {file}')
    with nc4.Dataset(file, 'r+') as f:

        # Resume after interrupts

        # Add wind direction, function argument 'dims' toggles a switch away from
        #  default behavior in make_nc_variable() which assumes lat & lon dims exist
        if 'hru' in f.dimensions:
            f = cs.derive_wind_direction(f,dims='hru')
        else:
            f = cs.derive_wind_direction(f)

print(debug_message)