import glob
import os
import sys
import pandas as pd
import xarray as xr
from datetime import datetime
from pathlib import Path
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs
config_file = '../0_config/config.txt'

# --- Reruns 2024-05-21
# These fix various small errors discovered during data use
rerun_file = Path('/globalhome/wmk934/HPC/camels_spat/7_forcing_data/forcing_check_logs/reruns_20240516.csv')
reruns = pd.read_csv(rerun_file)
# --- Reruns 2024-05-21

# Get the required info from the config file
data_path = cs.read_from_config(config_file,'data_path')

# CAMELS-spat metadata
cs_meta_path = cs.read_from_config(config_file,'cs_basin_path')
cs_meta_name = cs.read_from_config(config_file,'cs_meta_name')
cs_unusable_name = cs.read_from_config(config_file,'cs_unusable_name')

# Basin folder
cs_basin_folder = cs.read_from_config(config_file, 'cs_basin_path')
basins_path = Path(data_path) / cs_basin_folder
cs_meta_path = Path(data_path) / cs_meta_path
cs_meta = pd.read_csv(cs_meta_path / cs_meta_name)

# Open list of unusable stations; Enforce reading IDs as string to keep leading 0's
cs_unusable = pd.read_csv(cs_meta_path / cs_unusable_name, dtype={'Station_id': object})

# loop
for ix, row in cs_meta.iterrows():
    
    # --- Reruns 2024-05-21
    this_basin = row['Country'] + '_' + row['Station_id']
    if this_basin not in reruns['basin'].values:
        #print(f'No reruns for basin {this_basin}. Exiting.')
        continue # with next station, because we have no reruns for this station. Error code 0: clean exit, no problems
    else:
        print(f'Fixing EM-Earth for basin {this_basin}.')
    # --- Reruns 2024-05-21

    # Check if we need to run downloads for this station at all
    missing = cs.flow_obs_unavailable(cs_unusable, row.Country, row.Station_id)
    if 'iv' in missing and 'dv' in missing: 
        print(f'No flow observations for basin {this_basin}. Exiting.')
        continue # with next station, because we have no observations at all for this station

    # Get forcing destination path
    raw_fold, lump_fold, dist_fold = cs.prepare_forcing_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder) # Returns folders only, not file names

    # Get EM-Earth files
    em_earth_files = sorted(glob.glob(str(raw_fold/'EM_Earth_[0-9][0-9][0-9][0-9]-[0-9][0-9].nc'))) # list

    # Loop through the files and remove the timebnds info
    for file in em_earth_files:
        if 'raw' in file and 'EM_Earth' in file:
            ds = xr.open_dataset(file)
            if 'time_bnds()' in ds.History:
                ds = ds.drop('time_bnds')
                ds = ds.drop_dims('nbnds')
                ds.to_netcdf(file)
            ds.close()

    print(f'Done fixing time_bnds in {this_basin}')
