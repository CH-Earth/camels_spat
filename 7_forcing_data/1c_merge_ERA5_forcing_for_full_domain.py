# Merge ERA5 data for full domain
# We have ERA5 data on a monthly basis, but spread over two files for different variables. Here we merge those into one monthly file for further processing.

# Code based on: https://github.com/CH-Earth/CWARHM/tree/main/3a_forcing/2_merge_forcing

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
temp_folder = cs.read_from_config(config_file, 'temp_path') # Should exist at this point
temp_folder = '/gpfs/tp/gwf/gwf_cmt/wknoben/camels_spat/temp/ERA5/1_download' # re-runs 2024-05-6
## Processing
# Find the files and sort them according to date
surface_files  = sorted( glob.glob(temp_folder + '/ERA5_*_surface_variables.nc') )
pressure_files = sorted( glob.glob(temp_folder + '/ERA5_*_pressure_variables.nc') )

debug_message = f'!!! Warning: Check debugging status:\n-Full run in progress'

# Loop over the files, check that they match temporally, and merge
print(debug_message)
for surface_file, pressure_file in zip(surface_files,pressure_files):
        
    # Check that dimensions match
    flag,msg = cs.compare_era5_netcdf_dimensions(surface_file, pressure_file)
    if not flag: 
        print(msg)
        continue # Print error and skip to next set of files if dimensions for the current files don't match

    # Merge
    merged_file = surface_file.replace('_surface_variables','')
    cs.merge_era5_surface_and_pressure_files(surface_file, pressure_file, merged_file)
    print(f'Merged into {merged_file}')

print(debug_message)