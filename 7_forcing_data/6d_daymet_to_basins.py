# Uses DataTool to subset Daymet data to basins
import sys
import shutil
import subprocess
import math
import pandas as pd
from datetime import datetime
from pathlib import Path
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

# Hard-coded for now:
dt_folder = Path('/globalhome/wmk934/HPC/datatool/')

# --- General flag that can later be used to fix errors
overwrite_existing = True

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

# Temporary download path
temp_folder = Path( cs.read_from_config(config_file, 'temp_path') )
temp_daymet_folder = temp_folder / 'daymet' 

# Variables for the RDRS extraction
daymet_vars = cs.read_from_config(config_file, 'daymet_vars')

# --- Data loading
# CAMELS-spat metadata file
cs_meta_path = Path(data_path) / cs_meta_path
cs_meta = pd.read_csv(cs_meta_path / cs_meta_name)

# Open list of unusable stations; Enforce reading IDs as string to keep leading 0's
cs_unusable = pd.read_csv(cs_meta_path / cs_unusable_name, dtype={'Station_id': object})

# --- Processing
# Get the array index from the command-line argument
if len(sys.argv) != 3:
    print("Usage: python 7b_rdrs_to_basins.py <array_index> <temp_dir>")
    sys.exit(1)
else:
    ix = int(sys.argv[1])
    tmp_dir = sys.argv[2]

# Define cs_meta row
row = cs_meta.iloc[ix] # needs to be between 0  and 1697

# Get shapefile path to determine download coordinates, and forcing destination path
basin_id, shp_lump_path, shp_dist_path, _, _ = cs.prepare_delineation_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder)
raw_fold, _, _ = cs.prepare_forcing_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder) # Returns folders only, not file names
print('--- Now running basin {}. {}'.format(ix, basin_id))

# Check if we need to run for this station at all
missing = cs.flow_obs_unavailable(cs_unusable, row.Country, row.Station_id)
if 'iv' in missing and 'dv' in missing: 
    print(f'No flow observations for basin {basin_id}. Exiting.')
    sys.exit(0) # with next station, because we have no observations at all for this station. Error code 0: clean exit, no problems

# Specifically for Daymet, create a temporary 'year' folder where we can keep 
#  these files before merging them to monthly later
temp_basin_folder = raw_fold / 'daymet_year'

# RERUNS 2024-06-06
# Remove the existing folder
if temp_basin_folder.exists() and overwrite_existing:
    print(f'Removing existing folder {temp_basin_folder}')
    shutil.rmtree(temp_basin_folder)
# End RERUNS

temp_basin_folder.mkdir(parents=True, exist_ok=True)

# Update 2024-06-18: We need to use the bounds after all (before
# we simply gave the shapefile to datatool) because:
# 1. Datatool doesn't subset accurately with rotated lat/lon grids:
#    it sometimes misses the edges
# 2. We later process the gridded files into basin averages with 
#    easymore: with rotated grids easymore cannot estimate the 
#    edges of the outer grid cells, and thus doesn't use those. 
#    Therefore we need an extra buffer so that easymore can create
#    a forcing grid that actually covers the entire basin.
#
# From shapefile, get bounding coordinates
bounds = cs.find_shapefile_bounds(shp_lump_path)

# Apply an appropriate buffer. Daymet is at 1km resolution, which
# is approximately 0.008983 degrees latitude, and 
# 1/(111.32 * cos(lat)) longitude
lat_buffer = 0.008983
lat = math.radians(row['Station_lat'])
lon_buffer = 1/(111.32 * math.cos(lat))

# Applying a 2-cell buffer should theoretically be enough, 
# but let's go for 3 anyway
n_buffer = 3
bounds[0] = bounds[0] - n_buffer*lon_buffer
bounds[1] = bounds[1] - n_buffer*lat_buffer
bounds[2] = bounds[2] + n_buffer*lon_buffer
bounds[3] = bounds[3] + n_buffer*lat_buffer

# From meta-data, get flow obs period
times_flow = cs.find_flow_obs_times_from_metadata(row, missing)
times_era5 = cs.round_flow_obs_to_days(times_flow)
start_date = datetime.strptime(times_era5[0], '%Y-%m-%d')
final_date = datetime.strptime(times_era5[1], '%Y-%m-%d')
print(f'    Basin coordinates:         {bounds}')
print(f'    Flow obs unavailable:      {missing}')
print(f'    Flow obs times:            {times_era5}')

# Limit the dates to the Daymet period
# According to Kasra, datatool should be able to handle periods that are outside what the data has
# BUT this does not apply to Daymet (datatool: 2024-06-01)
if start_date.year < 1980:
    start_date = datetime(1980,1,1)
if final_date.year > 2023:
    final_date = datetime(2023,12,31)

# Fix a further issue where datatool doesn't like dates that are not on the start/end of the year
start_date = start_date.replace(day=1, month=1)
final_date = final_date.replace(day=31, month=12)

# Convert to strings
dt_start_date = start_date.strftime('%Y-%m-%d %H:%M:%S')
dt_final_date = final_date.strftime('%Y-%m-%d %H:%M:%S')

# -- DEV
#tmp_dir = Path('/scratch/gwf/gwf_cmt/wknoben/daymet_cache')
#tmp_dir.mkdir(parents=True, exist_ok=True)
# -- END DEV

# Create the datatool string
dt_string = f'{dt_folder}/extract-dataset.sh' \
            f' --dataset="daymet"' \
            f' --dataset-dir="{str(temp_daymet_folder)}"' \
            f' --output-dir="{str(temp_basin_folder)}"' \
            f' --start-date="{dt_start_date}"' \
            f' --end-date="{dt_final_date}"' \
            f' --variable="{daymet_vars}"'\
            f' --prefix="{basin_id}_subset_"' \
            f' --cache="{tmp_dir}"' \
            f' --lat-lims={bounds[1]},{bounds[3]}' \
            f' --lon-lims={bounds[0]},{bounds[2]}'
            #f' --shape-file="{str(shp_lump_path)}"'\
                      
# use subprocess to run the datatool command
print(f'Running datatool command: {dt_string}')
result = subprocess.run(dt_string, shell=True, 
                        capture_output=True, text=True)
print(result.stdout)
print(result.stderr)

# --- Check the results
out_files = list(temp_basin_folder.glob(f'{basin_id}_subset_*.nc'))
if len(out_files) == 0:
    print(f'\nNo files found in {temp_basin_folder} after datatool. Exiting.')
    sys.exit(1)

expected_vars = daymet_vars.split(',')
start_year = str(max(1980,start_date.year)) # limit to Daymet availability
end_year = str(min(2023,final_date.year))
for var in expected_vars:
    # Find the files with this variable in the name
    var_files = [f for f in out_files if var in f.name]
    var_files.sort()
    if len(var_files) == 0:
        print(f'Error: no files found for variable {var}')
        sys.exit(1)
    # Check that we have the right start and end files
    if start_year not in var_files[0].name:
        print(f'Error: expected {start_year} as first file but got {var_files[0].name} for variable {var}')
        sys.exit(1)
    if end_year not in var_files[-1].name:
        print(f'Error: expected {end_year} as last file but got {var_files[-1].name} for variable {var}')
        sys.exit(1)
    # Check that we have the right number of files
    expected_files = int(end_year) - int(start_year) + 1
    if len(var_files) != expected_files:
        print(f'Error: expected {expected_files} files but got {len(var_files)} for variable {var}')
        sys.exit(1)


'''
datatool extract example

./extract-dataset.sh  \
  --dataset="rdrs" \
  --dataset-dir="/project/rpp-kshook/Climate_Forcing_Data/meteorological-data/rdrsv2.1" \
  --output-dir="$HOME/scratch/rdrs_outputs/" \
  --start-date="2001-01-01 00:00:00" \
  --end-date="2001-12-31 23:00:00" \
  --lat-lims=49,51  \
  --lon-lims=-117,-115 \
  --variable="RDRS_v2.1_A_PR0_SFC,RDRS_v2.1_P_HU_09944" \
  --prefix="testing_"
'''