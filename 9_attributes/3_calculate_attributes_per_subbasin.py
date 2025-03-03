# Calculate basin attributes for each subbasin
# Accepts basin ID as command line input

import geopandas as gpd
import numpy as np
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

# Get the temporary data folder
cs_temp_folder = cs.read_from_config(config_file, 'temp_path')
temp_path = Path(cs_temp_folder)
temp_path.mkdir(exist_ok=True, parents=True)

# Get the attribute folder
att_folder = cs.read_from_config(config_file, 'att_path')
att_path = basins_path / att_folder / 'distributed'
att_path.mkdir(parents=True, exist_ok=True)

# Get the CRS we want to do area calculations in (possibly helpful for GLHYMPS and HydroLAKES)
ea_crs = cs.read_from_config(config_file, 'equal_area_crs')


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
# Specifically set the case
case = 'distributed'

# Data sets to use: Removed 'hydrology', because we don't have gauges for every subbasin
data_subfolders = ['rdrs', 'worldclim', 'lai', 'forest_height', 'glclu2019', 'modis_land', 'lgrip30', 'merit', 'hydrolakes', 'pelletier', 'soilgrids', 'glhymps']

# Get the paths
basin_id, _, shp_dist_path, _, _ = prepare_delineation_outputs(cs_meta, ix, basins_path)
geo_folder = basins_path / 'basin_data' / basin_id / 'geospatial'
met_folder = basins_path / 'basin_data' / basin_id / 'forcing'
hyd_folder = basins_path / 'basin_data' / basin_id / 'observations'

# Define the shapefiles
shp = str(shp_dist_path).format('basin') # because zonalstats wants a file path, not a geodataframe
riv = str(shp_dist_path).format('river') # For topographic attributes

# Make a temporary folder for storage
temp_path = temp_path / basin_id 
temp_path.mkdir(exist_ok=True, parents=True)

# Load the shapefile to get the sub-basin order for geotiffs and shapefiles
gdf = gpd.read_file(shp)

# Data storage
l_comids_geo = gdf['COMID'].to_list() # Store the polygon IDs into a list, we'll later use this as columns in the attribute dataframes
l_values_geo = [[] for _ in range(len(l_comids_geo))] # Initialize an empty list where we'll store this basin's attributes - nested lists for each subbasin
l_index_geo = [] # Initialize an empty list where we'll store the attribute descriptions - this will be our dataframe index

l_values_lakes = [[] for _ in range(len(l_comids_geo))] # HydroLAKES open water bodies shapefile
l_index_lakes = []
l_values_glhymps = [[] for _ in range(len(l_comids_geo))] # GLHYMPS geology shapefile
l_index_glhymps = []

l_values_met = [] # RDRS meteorological netcdf data - creating a nested list with subbasins is done differently for RDRS, so this line is supposed to look different from the other l_values_[x] lists
l_index_met = []

# Data-specific processing
print(f'Processing geospatial data into attributes for {basin_id}')
for dataset in data_subfolders:
    print(f' - processing {dataset}')

    ## CLIMATE
    if dataset == 'rdrs':
        #l_values, l_index, ds_precip, ds_rdrs = csa.attributes_from_rdrs(met_folder, shp, dataset, l_values, l_index)
        l_values_met, l_index_met, l_comids_met, _ = csa.attributes_from_rdrs(met_folder, shp, dataset, l_values_met, l_index_met)
    if dataset == 'worldclim':
        csa.oudin_pet_from_worldclim(geo_folder, dataset) # Get an extra PET estimate to sanity check RDRS outcomes
        csa.aridity_and_fraction_snow_from_worldclim(geo_folder, dataset) # Get monthly aridity and fraction snow maps
        l_values_geo, l_index_geo = csa.attributes_from_worldclim(geo_folder, dataset, shp, l_values_geo, l_index_geo, case=case)

    ## LAND COVER
    if dataset == 'forest_height':
        l_values_geo, l_index_geo = csa.attributes_from_forest_height(geo_folder, dataset, shp, l_values_geo, l_index_geo, case=case)
    if dataset == 'lai':
        l_values_geo, l_index_geo = csa.attributes_from_lai(geo_folder, dataset, temp_path, shp, l_values_geo, l_index_geo, case=case)
    if dataset == 'glclu2019':
        l_values_geo, l_index_geo = csa.attributes_from_glclu2019(geo_folder, dataset, shp, l_values_geo, l_index_geo, case=case)
    if dataset == 'modis_land':
        l_values_geo, l_index_geo = csa.attributes_from_modis_land(geo_folder, dataset, shp, l_values_geo, l_index_geo, case=case)
    if dataset == 'lgrip30':
        l_values_geo, l_index_geo = csa.attributes_from_lgrip30(geo_folder, dataset, shp, l_values_geo, l_index_geo, case=case)

    ## TOPOGRAPHY
    if dataset == 'merit':
        l_values_geo, l_index_geo, l_comids_merit = csa.attributes_from_merit(geo_folder, dataset, shp, riv, row, l_values_geo, l_index_geo, equal_area_crs=ea_crs, case=case)
        assert (l_comids_geo == l_comids_merit).all(), f"mismatch between COMID orders determined before and inside attributes_from_merit()"

    ## OPENWATER
    if dataset == 'hydrolakes':
        #l_values, l_index = csa.attributes_from_hydrolakes(geo_folder, dataset, l_values, l_index)
        l_values_lakes, l_index_lakes, l_comids_lakes = csa.attributes_from_hydrolakes(geo_folder, dataset, shp, ea_crs, 
                                                                                       l_values_lakes, l_index_lakes)

    ## SOIL
    if dataset == 'pelletier':
        l_values_geo, l_index_geo = csa.attributes_from_pelletier(geo_folder, dataset, shp, l_values_geo, l_index_geo, case=case)
    if dataset == 'soilgrids':
        l_values_geo, l_index_geo = csa.attributes_from_soilgrids(geo_folder, dataset, shp, l_values_geo, l_index_geo, case=case)

    ## GEOLOGY
    if dataset == 'glhymps':
        #l_values, l_index = csa.attributes_from_glhymps(geo_folder, dataset, l_values, l_index)
        l_values_glhymps, l_index_glhymps, l_comids_glhymps = csa.attributes_from_glhymps(geo_folder, dataset, shp, 
                                                                                          l_values_glhymps, l_index_glhymps, 
                                                                                          equal_area_crs=ea_crs)

## MERGE THE DATAFRAMES
# Make the individual dataframes
# - RDRS
att_df_met = pd.DataFrame(data = l_values_met, index = l_comids_met).transpose()
multi_index = pd.MultiIndex.from_tuples(l_index_met, names=['Category', 'Attribute', 'Unit', 'Source'])
att_df_met.index = multi_index

# - geotiffs
att_df_geo = pd.DataFrame(data = l_values_geo, index = l_comids_geo).transpose()
multi_index = pd.MultiIndex.from_tuples(l_index_geo, names=['Category', 'Attribute', 'Unit', 'Source'])
att_df_geo.index = multi_index

# - HydroLAKES
att_df_lakes = pd.DataFrame(data = l_values_lakes, index = l_comids_lakes).transpose()
multi_index = pd.MultiIndex.from_tuples(l_index_lakes, names=['Category', 'Attribute', 'Unit', 'Source'])
att_df_lakes.index = multi_index

# - GLHYMPS
att_df_glhymps = pd.DataFrame(data = l_values_glhymps, index = l_comids_glhymps).transpose()
multi_index = pd.MultiIndex.from_tuples(l_index_glhymps, names=['Category', 'Attribute', 'Unit', 'Source'])
att_df_glhymps.index = multi_index

# Check we have the same columns in all dfs
geo_columns = att_df_geo.columns.unique().sort_values()
met_columns = att_df_met.columns.unique().sort_values()
lak_columns = att_df_lakes.columns.unique().sort_values()
glh_columns = att_df_glhymps.columns.unique().sort_values()
assert (geo_columns == met_columns).all(), f"COMID mismatches between meteo and geospatial attribute dataframes"
assert (geo_columns == lak_columns).all(), f"COMID mismatches between hydrolakes and geospatial attribute dataframes"
assert (geo_columns == glh_columns).all(), f"COMID mismatches between glhymps and geospatial attribute dataframes"

# Merge
att_df = pd.concat([att_df_met,att_df_geo,att_df_lakes,att_df_glhymps])

# Save
att_file = f'attributes_{basin_id}.csv'
att_df.to_csv(att_path/att_file)

# Remove the temporary folder
shutil.rmtree(temp_path)









