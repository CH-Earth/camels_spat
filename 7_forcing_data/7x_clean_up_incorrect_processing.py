import glob
import os
import pandas as pd
from pathlib import Path
import shutil
import sys
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

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

# --- Loop over the basins and remove any old RDRS files
for ix, row in cs_meta.iterrows():
    basin_id = row['Country'] + '_' + row['Station_id']
    
    # Check if we need to run for this station at all
    missing = cs.flow_obs_unavailable(cs_unusable, row.Country, row.Station_id)
    if 'iv' in missing and 'dv' in missing: 
        print(f'No flow observations for basin {basin_id}. Exiting.')
        continue # with next station, because we have no observations at all for this station. Error code 0: clean exit, no problems
    print('--- Now running basin {}. {}'.format(ix, basin_id))
    
    # Define where the raw forcing folder is
    raw_fold, _, _ = cs.prepare_forcing_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder) # Returns folders only, not file names

    # Remove everything RDRS-related
    # rdrs_day_tar: folder with backed up daily files
    shutil.rmtree(raw_fold / 'rdrs_day_tar', ignore_errors=True)

    # rdrs_month: folder with monthly files
    shutil.rmtree(raw_fold / 'rdrs_month', ignore_errors=True)

    # rdrs_backup_before_units_pet: folder with backed up daily files before unit conversion
    shutil.rmtree(raw_fold / 'rdrs_backup_before_units_pet', ignore_errors=True)

    # "RDRS_*.nc" files in the raw folder
    rdrs_files = glob.glob(str(raw_fold / 'RDRS_*.nc'))
    for f in rdrs_files:
        os.remove(f)