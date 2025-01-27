# Calculate basin attributes
# Accepts basin ID as command line input

import numpy as np
import os
import pandas as pd
import shutil
import sys
from pathlib import Path
sys.path.append(str(Path().absolute().parent))
from python_cs_functions import config as cs, attributes as csa
from python_cs_functions.delineate import prepare_delineation_outputs

# --- Config handling
# Specify where the config file can be found
config_file = '../0_config/config.txt'

# Get the required info from the config file
data_path            = cs.read_from_config(config_file,'data_path')

# CAMELS-spat metadata
cs_meta_path = cs.read_from_config(config_file,'cs_basin_path')
cs_meta_name = cs.read_from_config(config_file,'cs_meta_name')
cs_unusable_name = cs.read_from_config(config_file,'cs_unusable_name')

# Basin folder
cs_basin_folder = cs.read_from_config(config_file, 'cs_basin_path')
basins_path = Path(data_path) / cs_basin_folder

# Get the attribute folder
att_folder = cs.read_from_config(config_file, 'att_path')
att_path = basins_path / att_folder
att_path.mkdir(parents=True, exist_ok=True)

# --- Data loading
# CAMELS-spat metadata file
cs_meta_path = Path(data_path) / cs_meta_path
cs_meta = pd.read_csv(cs_meta_path / cs_meta_name)

# Open list of unusable stations; Enforce reading IDs as string to keep leading 0's
cs_unusable = pd.read_csv(cs_meta_path / cs_unusable_name, dtype={'Station_id': object})

# --- Command line arguments
# Get the array index from the command-line argument
if len(sys.argv) != 2:
    print("Usage: python 1_calculate_attributes.py <array_index>")
    sys.exit(1)
else:
    ix = int(sys.argv[1])

# Get metadata
row = cs_meta.iloc[ix]

# Skip the basins for which we don't have daily flow obs
if (row['Station_id'] == cs_unusable['Station_id']).any():
    ixs = np.where(row['Station_id'] == cs_unusable['Station_id'])
    unusable = cs_unusable.iloc[ixs]    
    for i,r in unusable.iterrows():
        if r['Missing'] == 'dv':
            print('no daily streamflow observations available for this basin.')
            sys.exit(0)
        
# --- Processing
# Names of the data folders with catchment maps
data_subfolders = ['rdrs', 'worldclim', 'hydrology', 'lai', 'forest_height', 
                   'glclu2019', 'modis_land', 'lgrip30', 'merit', 'hydrolakes', 
                   'pelletier', 'soilgrids', 'glhymps'] # 

# Every attribute needs a list, so that we can efficiently construct a dataframe later
l_gauges = [] # station ID

# Get the paths
basin_id, shp_lump_path, shp_dist_path, _, _ = prepare_delineation_outputs(cs_meta, ix, basins_path)
geo_folder = basins_path / 'basin_data' / basin_id / 'geospatial'
met_folder = basins_path / 'basin_data' / basin_id / 'forcing'
hyd_folder = basins_path / 'basin_data' / basin_id / 'observations'

# Check if we have already completed this basin
att_file = f'attributes_{basin_id}.csv'
#if os.path.exists(att_path/att_file):
#    print(f'{att_file} exists. Skipping.')
#    sys.exit(0)

# Data storage
l_gauges.append(basin_id) # Update the Station list
l_values = [] # Initialize an empty list where we'll store this basin's attributes
l_index = [] # Initialize an empty list where we'll store the attribute descriptions

# Make a temporary folder for storage
temp_path = basins_path / 'basin_data' / basin_id / 'temp'
temp_path.mkdir(exist_ok=True, parents=True)

# Define the shapefiles
shp = str(shp_lump_path) # because zonalstats wants a file path, not a geodataframe
riv = str(shp_dist_path).format('river') # For topographic attributes

# Data-specific processing
print(f'Processing geospatial data into attributes for {basin_id}')
for dataset in data_subfolders:
    print(f' - processing {dataset}')
    
    ## CLIMATE
    if dataset == 'era5':
        l_values, l_index, ds_precip, ds_era5 = csa.attributes_from_era5(met_folder, shp, 'era5', l_values, l_index)
    if dataset == 'rdrs':
        l_values, l_index, ds_precip, ds_rdrs = csa.attributes_from_rdrs(met_folder, shp, dataset, l_values, l_index)
    if dataset == 'worldclim':
        csa.oudin_pet_from_worldclim(geo_folder, dataset) # Get an extra PET estimate to sanity check ERA5 outcomes
        csa.aridity_and_fraction_snow_from_worldclim(geo_folder, dataset) # Get monthly aridity and fraction snow maps
        l_values, l_index = csa.attributes_from_worldclim(geo_folder, dataset, shp, l_values, l_index)
    
    ## LAND COVER
    if dataset == 'forest_height':
        l_values, l_index = csa.attributes_from_forest_height(geo_folder, dataset, shp, l_values, l_index)
    if dataset == 'lai':
        l_values, l_index = csa.attributes_from_lai(geo_folder, dataset, temp_path, shp, l_values, l_index)
    if dataset == 'glclu2019':
        l_values, l_index = csa.attributes_from_glclu2019(geo_folder, dataset, shp, l_values, l_index)
    if dataset == 'modis_land':
        l_values, l_index = csa.attributes_from_modis_land(geo_folder, dataset, shp, l_values, l_index)
    if dataset == 'lgrip30':
        l_values, l_index = csa.attributes_from_lgrip30(geo_folder, dataset, shp, l_values, l_index)
    
    ## TOPOGRAPHY
    if dataset == 'merit':
        l_values, l_index = csa.attributes_from_merit(geo_folder, dataset, shp, riv, row, l_values, l_index)
    
    ## OPENWATER
    if dataset == 'hydrolakes':
        l_values, l_index = csa.attributes_from_hydrolakes(geo_folder, dataset, l_values, l_index)
    if dataset == 'hydrology':
        l_values, l_index = csa.attributes_from_streamflow(hyd_folder, dataset, basin_id, ds_precip, row, l_values, l_index)
    
    ## SOIL
    if dataset == 'pelletier':
        l_values, l_index = csa.attributes_from_pelletier(geo_folder, dataset, shp, l_values, l_index)
    if dataset == 'soilgrids':
        l_values, l_index = csa.attributes_from_soilgrids(geo_folder, dataset, shp, l_values, l_index)
    
    ## GEOLOGY
    if dataset == 'glhymps':
        l_values, l_index = csa.attributes_from_glhymps(geo_folder, dataset, l_values, l_index)

# Remove the temporary folder
shutil.rmtree(temp_path)

# --- Make the dataframe
# Add a fake second station, because otherwise we can't make a dataframe
l_gauges = [basin_id,'tmp_column']

# Make the dataframe
input_dict = dict(zip(l_gauges, [l_values,l_values]))
df = pd.DataFrame(input_dict)

# Set the index
multi_index = pd.MultiIndex.from_tuples(l_index, names=['Category', 'Attribute', 'Unit', 'Source'])
df.index = multi_index

# Drop the fake extra column
df = df.drop(columns=['tmp_column'], axis=1)

# Save to file
att_file = f'attributes_{basin_id}.csv'
df.to_csv(att_path/att_file)
        