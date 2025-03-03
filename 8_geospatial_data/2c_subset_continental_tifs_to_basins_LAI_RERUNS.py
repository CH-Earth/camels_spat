# Subset GeoTIFFs to basin shapefile outlines
import glob
import os
import pandas as pd
import sys
from pathlib import Path
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

# --- Config handling
# Specify where the config file can be found
config_file = '../0_config/config.txt'

# Get the required info from the config file
data_path            = cs.read_from_config(config_file,'data_path')
geospatial_temp_path = cs.read_from_config(config_file,'geospatial_temp_path')

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

# --- Find data files
data_folder = Path(data_path) / geospatial_temp_path
data_files = glob.glob( str(data_folder/'lai/raw/*.tif')) # E.g., geospatial_temp/lai

# --- Processing
# Get the array index from the command-line argument
if len(sys.argv) != 2:
    print("Usage: python 2c_subset_continental_tifs_to_basins_LAI_RERUNS.py <array_index>")
    sys.exit(1)
else:
    ix = int(sys.argv[1])

# Get metadata
row = cs_meta.iloc[ix]

# Get shapefile path to determine download coordinates, and forcing destination path
basin_id, shp_lump_path, _, _, _ = cs.prepare_delineation_outputs(cs_meta, ix, basins_path)
print(f'Processing LAI GeoTIFFs for {basin_id}')

# Loop over the files we want to subset
for file in data_files:
    print(f' - {file}')

    # Get the relative path compared to download folder; 
    # In other words, find which folders we want to create for the basin
    relative_path = Path(file).relative_to(data_folder)
    des_folder = basins_path / 'basin_data' / basin_id / 'geospatial' / os.path.dirname(relative_path)
    
    # Subset the file
    # 'buffer' adds a small, data-set dependent, buffer around the shapefile to ensure full coverage
    cs.subset_geotiff_to_shapefile(file,shp_lump_path,des_folder, buffer=True)
