# Merge EM-Earth data for full domain
import os
import re
import glob
import sys
import pandas as pd
from pathlib import Path
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

# --- Config handling
# Specify where the config file can be found
config_file = '../0_config/config.txt'

# Get the required info from the config file
temp_folder = cs.read_from_config(config_file, 'temp_path') # Should exist at this point
temp_folder = '/gpfs/tp/gwf/gwf_cmt/wknoben/camels_spat/temp_em_earth' # Overwrite with actual location, so we can run a few things side-by-side

# --- Processing
# Find the EM-Earth files and ensure they are sorted
em_earth_fold = temp_folder + '/EM_Earth_v1/deterministic_hourly'
p_files = sorted(glob.glob(str(em_earth_fold + '/prcp/NorthAmerica/*.nc'))) # list
t_files = sorted(glob.glob(str(em_earth_fold + '/tmean/NorthAmerica/*.nc'))) # list

# Make a folder for the merged files
merged_path = Path(em_earth_fold) / 'merged'
merged_path.mkdir(exist_ok=True, parents=True)

# --- Prepare file name replacement so that we have consistent notations between ERA5 files (YYYY-MM) and EM-Earth
# Define a regular expression pattern to match the format 'YYYYMM'
pattern = r'(\d{4})(\d{2})'

# Define a replacement pattern with groups to insert '-'
replacement = r'\1-\2'

debug_message = f'!!! Warning: Check debugging status:\n-Running with fixed EM-Earth file path\n-Full run in progress'

# Loop over the files, check that they match temporally, and merge
print(debug_message)
for p_file, t_file in zip(p_files,t_files):
        
    # Check that dimensions match
    flag,msg = cs.compare_em_earth_netcdf_dimensions(p_file, t_file)
    if not flag: 
        print(msg)
        continue # Print error and skip to next set of files if dimensions for the current files don't match

    # Merge
    file_name = re.sub(pattern, replacement, os.path.basename(p_file)) # Use re.sub to replace YYYYMM format with YYYY-MM
    merged_file = merged_path / file_name
    cs.merge_em_earth_prcp_and_tmean_files(p_file,t_file,merged_file)

print(debug_message)