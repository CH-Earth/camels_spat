# Add derived variables
import glob
import sys
import pandas as pd
from pathlib import Path
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

# --- Reruns 2024-05-19
# These fix various small errors discovered during data use
rerun_file = Path('/globalhome/wmk934/HPC/camels_spat/7_forcing_data/forcing_check_logs/reruns_20240516.csv')
reruns = pd.read_csv(rerun_file)
# --- Reruns 2024-05-19

# --- Config Handling
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
debug_message = f'\n!!! CHECK DEBUGGING STATUS: \n- Full run in progress\n'

# Get the array index from the command-line argument
if len(sys.argv) != 2:
    print("Usage: python 1e_add_derived_variables.py <array_index>")
    sys.exit(1)
else:
    ix = int(sys.argv[1])

# Define cs_meta row
row = cs_meta.iloc[ix] # needs to be between 0  and 1697

# --- Reruns 2024-05-19
this_basin = row['Country'] + '_' + row['Station_id']
if this_basin not in reruns['basin'].values:
    print(f'No reruns for basin {this_basin}. Exiting.')
    sys.exit(0) # with next station, because we have no reruns for this station. Error code 0: clean exit, no problems
else:
    print(f'Running reruns for basin {this_basin}.')
# --- Reruns 2024-05-19

print(debug_message)
    
# Get shapefile path to determine download coordinates, and forcing destination path
basin_id, shp_lump_path, shp_dist_path, _, _ = cs.prepare_delineation_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder)
raw_fold, _, _ = cs.prepare_forcing_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder) # Returns folders only, not file names
print('--- Now running basin {}. {}'.format(ix, basin_id))

# Check if we need to run downloads for this station at all
missing = cs.flow_obs_unavailable(cs_unusable, row.Country, row.Station_id)
if 'iv' in missing and 'dv' in missing: 
    sys.exit(0) # gracefully exit run, because we have no observations at all for this station

# Find the files
era5_merged_files = sorted(glob.glob(str(raw_fold/'ERA5_[0-9][0-9][0-9][0-9]-[0-9][0-9].nc'))) # list

# Loop over the files and add new variables
for era5_merged_file in era5_merged_files:
    cs.add_derived_variables(era5_merged_file)

print(debug_message)