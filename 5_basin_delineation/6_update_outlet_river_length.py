# Load river shapefiles and calculate the length of each river segment. Save the updated shapefiles.
import sys
import geopandas as gpd
import pandas as pd
from pathlib import Path
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

# --- CONFIG HANDLING
# Specify where the config file can be found
config_file = '../0_config/config.txt'

# Get the required info from the config file
data_path     = cs.read_from_config(config_file,'data_path')

# CAMELS-spat metadata
cs_meta_path  = cs.read_from_config(config_file,'cs_basin_path')
cs_meta_name  = cs.read_from_config(config_file,'cs_meta_name')

# Basin folder
cs_basin_folder = cs.read_from_config(config_file, 'cs_basin_path')

# Equal area projection
equal_area_proj = cs.read_from_config(config_file, 'equal_area_crs')

# --- DATA LOADING
# CAMELS-spat metadata file
cs_meta_path = Path(data_path) / cs_meta_path
cs_meta = pd.read_csv(cs_meta_path / cs_meta_name)

# --- PROCESSING
for i,row in cs_meta.iterrows():
    gauge_id, _, shp_dist_path, _, _ = cs.prepare_delineation_outputs(cs_meta, i, Path(data_path)/cs_basin_folder)
    # Load the river shapefile
    shp_riv = gpd.read_file(str(shp_dist_path).format('river'))  
    # Calculate the length of each river segment
    eap_riv = shp_riv.to_crs(equal_area_proj)
    eap_lkm = eap_riv['geometry'].length / 1000
    shp_riv['new_len_km'] = eap_lkm
    max_diff = (shp_riv['lengthkm'] - shp_riv['new_len_km']).max()
    # Save the updated shapefiles
    shp_riv.to_file(str(shp_dist_path).format('river'))
    print(f'{i:04}. {gauge_id}: Updated river shapefile lengths. Max difference: {max_diff:.2f} km')