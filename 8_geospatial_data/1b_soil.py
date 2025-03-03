import os
import shutil
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
soil_path            = cs.read_from_config(config_file,'soil_path')
soil_url             = cs.read_from_config(config_file,'soil_url')
download_area        = cs.read_from_config(config_file,'geospatial_area')

# --- Convert geospatial coordinates to specific subsetting coordinates for this data product
subset_window = cs.geospatial_coordinates_to_download_coordinates(download_area, 'soilgrids')

# --- Initial folder settings
download_folder = Path(data_path) / geospatial_temp_path / 'soilgrids' / 'download' # downloads go into subfolders here
raw_folder =  Path(data_path) / geospatial_temp_path / 'soilgrids' / 'raw' # final outputs go into subfolders here
download_folder.mkdir(parents=True, exist_ok=True)
raw_folder.mkdir(parents=True, exist_ok=True)

# --- Processing
# Get the array index from the command-line argument
if len(sys.argv) != 2:
    print("Usage: python 3a_average_focring_to_hrus.py <array_index>")
    sys.exit(1)
else:
    ix = int(sys.argv[1])

# Create the download lists
fields = ['bdod','cfvo','clay','sand','silt','soc'] # wrb and landmask are have different subfolders
depths = ['0-5cm','5-15cm','15-30cm','30-60cm','60-100cm','100-200cm']
layers = ['mean','Q0.05','Q0.5','Q0.95','uncertainty']

all_fields = []
all_depths = []
all_layers = []
for field in fields:
    for depth in depths:
        for layer in layers:
            all_fields.append(field)
            all_depths.append(depth)
            all_layers.append(layer)

# Select the one we're after, using the array job ID
field = all_fields[ix]
depth = all_depths[ix]
layer = all_layers[ix]

# Prepare the folders and file names
product = f'{field}_{depth}_{layer}'
download_destination = download_folder/product
output_folder = raw_folder/field
output_file = f'{product}.tif'

# Ensure we don't duplicate runs
#if os.path.isfile(output_folder/output_file):
#    print(f'{output_folder/output_file} already exists.')
#    sys.exit(0)

# Download the data
download_url = f'{soil_url}{field}/{product}/' # trailing '/' is critical or urljoin will ignore the final part later
cs.download_all_soilgrids_tiles_into_folder(download_url, download_destination)

# Process individual tiles into a single GeoTIFF of the domain of interest
cs.process_soilgrids_tiles_into_single_geotiff(download_destination, 
                                            output_folder, output_file, 
                                            to_crs='EPSG:4326', subset_window=subset_window)