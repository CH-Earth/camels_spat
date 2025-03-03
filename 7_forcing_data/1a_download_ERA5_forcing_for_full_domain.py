import sys
import concurrent.futures
from datetime import datetime, timedelta
import pandas as pd
from pathlib import Path
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

#### Config handling
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

# Merged shape
merged_shape_folder = cs.read_from_config(config_file, 'merged_shp_dir')
merged_shape_name   = cs.read_from_config(config_file, 'merged_shp_name')

# Temporary download path
temp_folder = Path( cs.read_from_config(config_file, 'temp_path') )

# Specify the paths
merged_shape_path = Path(data_path) / cs_basin_folder / merged_shape_folder

temp_folder.mkdir(parents=True, exist_ok=True)


#### Data loading
# CAMELS-spat metadata file
cs_meta_path = Path(data_path) / cs_meta_path
cs_meta = pd.read_csv(cs_meta_path / cs_meta_name)

# Open list of unusable stations; Enforce reading IDs as string to keep leading 0's
cs_unusable = pd.read_csv(cs_meta_path / cs_unusable_name, dtype={'Station_id': object})

#### Processing
## Define the spatial extent of the whole domain
# Find bounding box
bounds = cs.find_shapefile_bounds(merged_shape_path/merged_shape_name)

# Convert bounding box to download coordinates
coords_era5, _, _ = cs.find_download_coords_from_bounds(bounds, target='ERA5')
print(f'Spatial extent of download domain \nBounds : {bounds} \nCoords : {coords_era5}')

## Define the temporal extent of the whole domain
# Start with some dates that we will overwrite with the true values in the next cell
start_date = '2023-12-31' 
final_date = '1900-01-01'

for ix,row in cs_meta.iterrows():
    
    # Check if any flow observations are missing
    missing = cs.flow_obs_unavailable(cs_unusable, row.Country, row.Station_id)
    if 'iv' in missing and 'dv' in missing: 
        continue # with next station, because we have no observations at all for this station
    
    # From meta-data, get download period
    times_flow = cs.find_flow_obs_times_from_metadata(row, missing)
    times_era5 = cs.round_flow_obs_to_days(times_flow)
    
    # Compare these times with what we already have, and update if we need to extend the download window on either side
    start_date = min(start_date, times_era5[0])
    final_date = max(final_date, times_era5[1])

# Manually fix the start date at a time before which we have the data already
#final_date = '2020-11-01'

# Convert to datetimes
start_date = datetime.strptime(start_date, '%Y-%m-%d')
final_date = datetime.strptime(final_date, '%Y-%m-%d')
print(f'Temporal extent of download domain \nStart : {start_date} \nEnd   : {final_date}')

#### Download the data
# Get the time-invariant data
cs.download_era5_time_invariant_data_to_netcdf(coords_era5, temp_folder/'ERA5'/'ERA5_2023-01-01_invariants.nc')

# Convert start and end dates into two lists of start and end dates, that we'll iterate over
start_list,end_list = cs.convert_start_and_end_dates_to_era5_download_lists(start_date,final_date)

# Function to call inside the parallel loop
def download_era5(start, end):
    
    # Convert to relevant
    yyyy_mm = start.strftime('%Y-%m') # filename
    start   = start.strftime('%Y-%m-%d') # yyyy-mm-dd for use with cdsapi
    end     = end.strftime('%Y-%m-%d')
    
    # Download
    cs.download_era5_surface_level_data_to_netcdf(coords_era5, start, end, temp_folder/'ERA5'/'1_download'/f'ERA5_{yyyy_mm}_surface_variables.nc')
    cs.download_era5_pressure_level_data_to_netcdf(coords_era5, start, end, temp_folder/'ERA5'/'1_download'/f'ERA5_{yyyy_mm}_pressure_variables.nc')
    
    return print(f'Download {yyyy_mm} successful')

# Number of threads or processes to use
num_workers = 5  # Adjust as needed

with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
    results = list( executor.map(download_era5, start_list, end_list) )