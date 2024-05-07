import glob
import os
import sys
import pandas as pd
from datetime import datetime
from pathlib import Path
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs
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
cs_meta_path = Path(data_path) / cs_meta_path
cs_meta = pd.read_csv(cs_meta_path / cs_meta_name)

# Open list of unusable stations; Enforce reading IDs as string to keep leading 0's
cs_unusable = pd.read_csv(cs_meta_path / cs_unusable_name, dtype={'Station_id': object})

# loop
for ix, row in cs_meta.iterrows():

    # Check if we need to run downloads for this station at all
    missing = cs.flow_obs_unavailable(cs_unusable, row.Country, row.Station_id)
    if 'iv' in missing and 'dv' in missing: 
        print(f'No flow observations for basin {basin_id}. Exiting.')
        continue # with next station, because we have no observations at all for this station
    
    # Get shapefile path to determine download coordinates, and forcing destination path
    basin_id, shp_lump_path, shp_dist_path, _, _ = cs.prepare_delineation_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder)
    raw_fold, lump_fold, dist_fold = cs.prepare_forcing_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder) # Returns folders only, not file names
    
    # Find the ERA5 and EM-Earth files we already have in the [basin] -> [forcing] -> [lump folder] & [dist_folder]
    lump_eme_files = sorted(glob.glob(str(lump_fold/'EM-Earth_*.nc'))) # list
    dist_eme_files = sorted(glob.glob(str(dist_fold/'EM-Earth_*.nc'))) # list
    lump_era_files = sorted(glob.glob(str(lump_fold/'ERA5_*.nc'))) # list
    dist_era_files = sorted(glob.glob(str(dist_fold/'ERA5_*.nc'))) # list
    
    # Remove the last file in each list, in case something has gone wrong with writing that
    print(f'Cleaning up {basin_id}')
    for list in [lump_eme_files, dist_eme_files, lump_era_files, dist_era_files]:
        if len(list) == 0:
            continue # because no files in this category have been remapped yet
        else:
            file = list[-1]
            print(f'   Removing {file}')
            os.remove(file)

print(f'Done removing files')