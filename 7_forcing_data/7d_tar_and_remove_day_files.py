import pandas as pd
from pathlib import Path
import shutil
import sys
import tarfile
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

# --- Command line arguments
if len(sys.argv) != 2:
    print("Usage: python 7c_merge_rdrs_to_month.py <array_index>")
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

# Define where the RDRS folder with daily data is
raw_fold, _, _ = cs.prepare_forcing_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder) # Returns folders only, not file names
temp_day_folder = raw_fold / 'rdrs_day'

# Find the yearly folders in the RDRS temp directory
year_folders = [f.name for f in temp_day_folder.iterdir() if f.is_dir()]

# --- Clean up
temp_tar_folder = raw_fold / 'rdrs_day_tar'
temp_tar_folder.mkdir(exist_ok=True, parents=True)
temp_tar_file = temp_tar_folder / f'{basin_id}_rdrs_day.tar.gz'

# find paths to all files
day_files = []
for year in year_folders:
    files = [f for f in (temp_day_folder / year).iterdir() if f.is_file()]
    day_files.extend(files)

# Tar the files
with tarfile.open(temp_tar_file, 'w:gz') as tar:
    for file in day_files:
        tar.add(temp_day_folder / year / file, arcname=file)

# Remove the individual day files
shutil.rmtree(temp_day_folder)