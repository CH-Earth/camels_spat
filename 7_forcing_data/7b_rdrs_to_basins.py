# Uses DataTool to subset RDRS data to basins
import glob
import os
import sys
import subprocess
import pandas as pd
from datetime import datetime
from pathlib import Path
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

# Hard-coded for now:
dt_folder = Path('/globalhome/wmk934/HPC/datatool/')

# --- General flag that can later be used to fix errors
overwrite_existing = False

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
temp_rdrs_folder = temp_folder / 'rdrs' # hard-coded like this in 7a.py; datatool expects the 1980,1981,.. subfolders

# Variables for the RDRS extraction
rdrs_vars = cs.read_from_config(config_file, 'rdrs_vars')

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

# Check if we need to run for this station at all
missing = cs.flow_obs_unavailable(cs_unusable, row.Country, row.Station_id)
if 'iv' in missing and 'dv' in missing: 
    print(f'No flow observations for basin {basin_id}. Exiting.')
    sys.exit(0) # with next station, because we have no observations at all for this station. Error code 0: clean exit, no problems

# Get shapefile path to determine download coordinates, and forcing destination path
basin_id, shp_lump_path, shp_dist_path, _, _ = cs.prepare_delineation_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder)
raw_fold, _, _ = cs.prepare_forcing_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder) # Returns folders only, not file names
print('--- Now running basin {}. {}'.format(ix, basin_id))

# Specifically for RDRS, create a temporary 'day' folder where we can keep 
#  these files before merging them to monthly later
temp_day_folder = raw_fold / 'rdrs_day'
temp_day_folder.mkdir(parents=True, exist_ok=True)

# From shapefile, get bounding coordinates
bounds = cs.find_shapefile_bounds(shp_lump_path)

# From meta-data, get flow obs period
times_flow = cs.find_flow_obs_times_from_metadata(row, missing)
times_era5 = cs.round_flow_obs_to_days(times_flow)
start_date = datetime.strptime(times_era5[0], '%Y-%m-%d')
final_date = datetime.strptime(times_era5[1], '%Y-%m-%d')
print(f'    Basin coordinates:         {bounds}')
print(f'    Flow obs unavailable:      {missing}')
print(f'    Flow obs times:            {times_era5}')

# According to Kasra, datatool should be able to handle periods that are outside what the data has
dt_start_date = start_date.strftime('%Y-%m-%d %H:%M:%S')
dt_final_date = final_date.strftime('%Y-%m-%d %H:%M:%S')

# For posterity: this can be used to finetune datatool input if turns out to be needed
# Set the RDRS data availability
# rdrs_start = datetime(1980,1,1,0,0)
# rdrs_final = datetime(2018,12,30,0,0)
# if start_date < rdrs_start:
#     dt_start_date = rdrs_start # datetime(1980,1,1,0,0) 
# else:
#     dt_start_date = start_date

# if final_date > rdrs_final:
#     dt_final_date = rdrs_final # datetime(2018,12,31,23,0)
# else:
#     dt_final_date = final_date

# dt_start_date = dt_start_date.strftime('%Y-%m-%d %H:%M:%S')
# dt_final_date = dt_final_date.strftime('%Y-%m-%d %H:%M:%S')

# Create the datatool string
dt_string = f'{dt_folder}/extract-dataset.sh' \
            f' --dataset="rdrs"' \
            f' --dataset-dir="{str(temp_rdrs_folder)}"' \
            f' --output-dir="{str(temp_day_folder)}"' \
            f' --start-date="{dt_start_date}"' \
            f' --end-date="{dt_final_date}"' \
            f' --variable="{rdrs_vars}"'\
            f' --prefix="{basin_id}_subset_"' \
            f' --shape-file="{str(shp_lump_path)}"'\
            f' --cache="{tmp_dir}"'
            #f' --lat-lims={bounds[1]},{bounds[3]}' \
            #f' --lon-lims={bounds[0]},{bounds[2]}'
                      
# use subprocess to run the datatool command
print(f'Running datatool command: {dt_string}')
result = subprocess.run(dt_string, shell=True, capture_output=True, text=True)
print(result)

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