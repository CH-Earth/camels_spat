# Loops over soilgrids data and calculates porosity and conducivity

# NOTE: we're running this AFTER we've already subsetted the continental maps to basins
# In a better world, we'd do this for the big map and subset later but it is what it is.

import numpy as np
import os
import pandas as pd
import rasterio
import sys
from pathlib import Path
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

# --- Processing
# Get the array index from the command-line argument
if len(sys.argv) != 2:
    print("Usage: python 3a_derive_soil_properties.py <array_index>")
    sys.exit(1)
else:
    ix = int(sys.argv[1])

# Define cs_meta row
row = cs_meta.iloc[ix] # needs to be between 0  and 1697

# Get shapefile path to determine download coordinates, and forcing destination path
basin_id, _, _, _, _ = cs.prepare_delineation_outputs(cs_meta, ix, basins_path)
print(f' - Processing {basin_id}')

# Define the soilgrids folders
soil_folder = basins_path / 'basin_data' / basin_id / 'geospatial' / 'soilgrids' / 'raw'
sand_folder = soil_folder / 'sand'
clay_folder = soil_folder / 'clay'

# Define the output folders
porosity_folder = soil_folder / 'porosity'
conductivity_folder = soil_folder / 'conductivity'
porosity_folder.mkdir(parents=True, exist_ok=True)
conductivity_folder.mkdir(parents=True, exist_ok=True)

# --- Processing
depths = ['0-5cm','5-15cm','15-30cm','30-60cm','60-100cm','100-200cm']
layer = 'mean'

for depth in depths:

    # Load the sand and clay geotiffs
    sand_file = sand_folder / f'sand_{depth}_{layer}.tif'
    clay_file = clay_folder / f'clay_{depth}_{layer}.tif'

    with rasterio.open(sand_file) as src:
        sand = np.ma.masked_equal(src.read(), src.nodata)
    with rasterio.open(clay_file) as src:
        clay = np.ma.masked_equal(src.read(), src.nodata)
        profile = src.profile
        profile.update(
            dtype=rasterio.float32,
            count=1,
            compress='lzw',
            nodata=src.nodata
        )

    # Convert the permille values to percent
    sand = sand / 10
    clay = clay / 10

    # Calculate porosity
    theta_s = 50.5 -0.142 * sand - 0.037 * clay # mean Theta_s (Table 4, Cosby et al., 1984), [%]
    porosity = theta_s / 100 # Convert to fraction to match CAMELS

    # Calculate conductivity
    log_ks = -0.6 + 0.0126 * sand - 0.0064 * clay # mean log Ks (Table 4, Cosby et al., 1984), [inch hr-1]
    conductivity = 2.54 * np.power(10,log_ks) # Convert to cm hr-1 to match CAMELS; can't use 10**log_ks because that will be applied to both data AND mask value

    # Save the results to new geotiffs
    porosity_file = porosity_folder / f'porosity_{depth}_{layer}.tif'
    conductivity_file = conductivity_folder / f'conductivity_{depth}_{layer}.tif'

    with rasterio.open(porosity_file, 'w', **profile) as dst:
        dst.write(porosity.filled().astype(rasterio.float32))
    
    with rasterio.open(conductivity_file, 'w', **profile) as dst:
        dst.write(conductivity.filled().astype(rasterio.float32))
