# Subset merged ERA5 data to basins
import glob
import os
import sys
import pandas as pd
from datetime import datetime
from pathlib import Path
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

# --- General flag that can later be used to debug
overwrite_existing = False

# --- Reruns 2024-05-18
# These fix various small errors discovered during data use
rerun_file = Path('/globalhome/wmk934/HPC/camels_spat/7_forcing_data/forcing_check_logs/reruns_20240516.csv')
reruns = pd.read_csv(rerun_file)
overwrite_existing = True
# --- Reruns 2024-05-18

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

# --- Data loading
# CAMELS-spat metadata file
cs_meta_path = Path(data_path) / cs_meta_path
cs_meta = pd.read_csv(cs_meta_path / cs_meta_name)

# Open list of unusable stations; Enforce reading IDs as string to keep leading 0's
cs_unusable = pd.read_csv(cs_meta_path / cs_unusable_name, dtype={'Station_id': object})

# --- Processing
# Find the ERA5 files
era5_invariant = glob.glob(str(temp_folder/'*invariants.nc')) # list
era5_merged = sorted(glob.glob(str(temp_folder/'ERA5_[0-9][0-9][0-9][0-9]-[0-9][0-9].nc'))) # list
era5_files = era5_invariant + era5_merged

# Get the array index from the command-line argument
if len(sys.argv) != 2:
    print("Usage: python 2c_subset_EM_Earth_to_basins.py <array_index>")
    sys.exit(1)
else:
    ix = int(sys.argv[1])

# Define cs_meta row
row = cs_meta.iloc[ix] # needs to be between 0  and 1697

# --- Reruns 2024-05-18
this_basin = row['Country'] + '_' + row['Station_id']
if this_basin not in reruns['basin'].values:
    print(f'No reruns for basin {this_basin}. Exiting.')
    sys.exit(0) # with next station, because we have no reruns for this station. Error code 0: clean exit, no problems
else:
    print(f'Running reruns for basin {this_basin}.')
# --- Reruns 2024-05-18

# Run
debug_message = f'\n!!! CHECK DEBUGGING STATUS: \n-Full array job run'
print(debug_message)

# Get shapefile path to determine download coordinates, and forcing destination path
basin_id, shp_lump_path, shp_dist_path, _, _ = cs.prepare_delineation_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder)
raw_fold, _, _ = cs.prepare_forcing_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder) # Returns folders only, not file names
print('--- Now running basin {}. {}'.format(ix, basin_id))

# From shapefile, get bounding coordinates. Then determine download coordinates from those
bounds = cs.find_shapefile_bounds(shp_lump_path)
coords_era5, _, _ = cs.find_download_coords_from_bounds(bounds, target='ERA5')

# Check if we need to run downloads for this station at all
missing = cs.flow_obs_unavailable(cs_unusable, row.Country, row.Station_id)
if 'iv' in missing and 'dv' in missing: 
    print(f'No flow observations for basin {basin_id}. Exiting.')
    sys.exit(0) # with next station, because we have no observations at all for this station. Error code 0: clean exit, no problems

# From meta-data, get download period
times_flow = cs.find_flow_obs_times_from_metadata(row, missing)
times_era5 = cs.round_flow_obs_to_days(times_flow)
start_date = datetime.strptime(times_era5[0], '%Y-%m-%d')
final_date = datetime.strptime(times_era5[1], '%Y-%m-%d')

print(f'    Basin coordinates:         {bounds}')
print(f'    ERA5 download coordinates: [{coords_era5}]')
print(f'    Flow obs unavailable:      {missing}')
print(f'    Download times:            {times_era5}')

# Convert start and end dates into two lists of start and end dates, that we'll iterate over
date_list,_ = cs.convert_start_and_end_dates_to_era5_download_lists(start_date,final_date) # not the cleanest but this lets us reuse old code
subset_strings = [date_obj.strftime("%Y-%m") for date_obj in date_list] # convert datetime objects to yyyy-mm strings

# Subset the data files
infiles = [temp_folder/'ERA5_2023-01-01_invariants.nc'] + [temp_folder/f'ERA5_{yyyy_mm}.nc' for yyyy_mm in subset_strings]
for infile in infiles:
    if os.path.exists(infile):
        outfile = raw_fold/os.path.basename(infile)
        if (not os.path.exists(outfile)) or (overwrite_existing): # file doesn't exist yet or we don't care
            cs.extract_ERA5_subset(infile,outfile,coords_era5)
    else:
        print(f'    ERROR: source file {infile} not found.')

# Create a figure to check if we actually cover the right domain with this
fig_file = raw_fold.parent / f'{row.Country}_{row.Station_id}_era5_coverage.png'
cs.compare_forcing_data_and_shape_extents(fig_file, outfile, shp_lump_path, nc_var='t', nc_time=0)

fig_file = Path('/gpfs/tp/gwf/gwf_cmt/wknoben/camels_spat/camels-spat-data/TEMP_images_forcing_coverage_era5') / f'{row.Country}_{row.Station_id}_era5_coverage.png'
cs.compare_forcing_data_and_shape_extents(fig_file, outfile, shp_lump_path, nc_var='t', nc_time=0)

print(debug_message)

# Leaving the below for posterity - might need it at some point
'''
from concurrent.futures import ThreadPoolExecutor # new
import threading # new

# Initialize file locking for nc file opening with xarray inside cs.extract_ERA5_subset()
file_lock = threading.Lock()

# Functions for parallel processing
def era5_subsetting_workflow(ix,row):
    
    # Get shapefile path to determine download coordinates, and forcing destination path
    basin_id, shp_lump_path, shp_dist_path, _, _ = cs.prepare_delineation_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder)
    raw_fold, _, _ = cs.prepare_forcing_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder) # Returns folders only, not file names
    print('--- Now running basin {}. {}'.format(ix, basin_id))
    
    # From shapefile, get bounding coordinates. Then determine download coordinates from those
    bounds = cs.find_shapefile_bounds(shp_lump_path)
    coords_era5, _, _ = cs.find_download_coords_from_bounds(bounds, target='ERA5')
    
    # Check if we need to run downloads for this station at all
    missing = cs.flow_obs_unavailable(cs_unusable, row.Country, row.Station_id)
    if 'iv' in missing and 'dv' in missing: 
        print(f'Skipping basin {basin_id} because there are no discharge observations available.') # with next station, because we have no observations at all for this station
        return
    
    # From meta-data, get download period
    times_flow = cs.find_flow_obs_times_from_metadata(row, missing)
    times_era5 = cs.round_flow_obs_to_days(times_flow)
    start_date = datetime.strptime(times_era5[0], '%Y-%m-%d')
    final_date = datetime.strptime(times_era5[1], '%Y-%m-%d')
    
    print(f'    Basin coordinates:         {bounds}')
    print(f'    ERA5 download coordinates: [{coords_era5}]')
    print(f'    Flow obs unavailable:      {missing}')
    print(f'    Download times:            {times_era5}')
    
    # Convert start and end dates into two lists of start and end dates, that we'll iterate over
    date_list,_ = cs.convert_start_and_end_dates_to_era5_download_lists(start_date,final_date) # not the cleanest but this lets us reuse old code
    subset_strings = [date_obj.strftime("%Y-%m") for date_obj in date_list] # convert datetime objects to yyyy-mm strings
    
    # Subset the data files
    infiles = [temp_folder/'ERA5_2023-01-01_invariants.nc'] + [temp_folder/f'ERA5_{yyyy_mm}.nc' for yyyy_mm in subset_strings]
    for infile in infiles:
        if os.path.exists(infile):
            outfile = raw_fold/os.path.basename(infile)
            if not os.path.exists(outfile):
                with file_lock: # ensure only a single thread can access the file
                    cs.extract_ERA5_subset(infile,outfile,coords_era5)
        else:
            print(f'    ERROR: source file {infile} not found.')
    
    # Create a figure to check if we actually cover the right domain with this
    # Note: this doesn't work well with parallel processing - we'll move it to its own file to run in serial
    #fig_file = raw_fold.parent / f'{row.Country}_{row.Station_id}_era5_coverage.png'
    #cs.compare_forcing_data_and_shape_extents(fig_file, outfile, shp_lump_path, nc_var='t', nc_time=0)
    return

# Start processing in parallel
num_threads = 24

# Create a ThreadPoolExecutor with the desired number of threads
debug_message = f'\n!!! CHECK DEBUGGING STATUS: \n- Running all files \n- Running all basins'
print(debug_message)
with ThreadPoolExecutor(max_workers=num_threads) as executor:
    # Use executor.map to parallelize row processing
    # Pass the process_row function and the DataFrame rows as arguments
    executor.map(era5_subsetting_workflow, cs_meta.index, cs_meta.itertuples(index=True))
print(debug_message)
'''