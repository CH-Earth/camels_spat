# Check datatool outcomes
from datetime import datetime,timedelta 
import glob
import pandas as pd
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

# --- RDRS source data
# RDRS data availbility
rdrs_s = datetime(1980,1,1,12,0)
rdrs_e = datetime(2018,12,30,12,0)

# RDRS files
rdrs_files = sorted(glob.glob(str(temp_rdrs_folder / '*' / '*.nc'))) # E.g.: '/gpfs/tp/gwf/gwf_cmt/wknoben/camels_spat/temp/rdrs/1980/1980010112.nc'
rdrs_dates = [datetime.strptime(file.split('/')[-1].split('.')[0], '%Y%m%d%H') for file in rdrs_files]

# --- Helper function
def are_dates_consecutive(dates, interval=timedelta(days=1)):
    sorted_dates = sorted(dates)
    for i in range(1, len(sorted_dates)):
        if sorted_dates[i] - sorted_dates[i - 1] != interval:
            return False
    return True

# --- Checks
write_rerun_file = False
rerun_file = '7b_rdrs_to_basins_rerun_with_start_file_missing.txt'
if write_rerun_file: f = open(rerun_file, 'w')

for ix,row in cs_meta.iterrows():
    
    # Check if we need to run for this station at all
    missing = cs.flow_obs_unavailable(cs_unusable, row.Country, row.Station_id)
    if 'iv' in missing and 'dv' in missing: 
        print(f'No flow observations for basin {basin_id}. Skipping.')
        continue

    # Get shapefile path to determine download coordinates, and forcing destination path
    basin_id, shp_lump_path, shp_dist_path, _, _ = cs.prepare_delineation_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder)
    raw_fold, _, _ = cs.prepare_forcing_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder) # Returns folders only, not file names
    #print('--- Now running basin {:0004}. {}'.format(ix, basin_id))

    # Specifically for RDRS, create a temporary 'day' folder where we can keep 
    #  these files before merging them to monthly later
    temp_basin_folder = raw_fold / 'rdrs_day'

    # From meta-data, get flow obs period
    times_flow = cs.find_flow_obs_times_from_metadata(row, missing)
    times_era5 = cs.round_flow_obs_to_days(times_flow)
    start_date = datetime.strptime(times_era5[0], '%Y-%m-%d')
    final_date = datetime.strptime(times_era5[1], '%Y-%m-%d')

    # Limit this to the RDRS availability
    if start_date < rdrs_s:
        start_date = rdrs_s
        
    if final_date > rdrs_e:
        final_date = rdrs_e

    # Find the subsetted RDRS files for this basin
    basin_rdrs = sorted(glob.glob(str(temp_basin_folder / '*' / '*.nc')))
    basin_dates = [datetime.strptime(file.split('/')[-1].split('.')[0].split('_')[-1], '%Y%m%d%H') for file in basin_rdrs]

    # Check if start and end date match
    if basin_dates[0] > start_date:
        print(f'    First subset file for {basin_id} starts after expected date. Expected: {start_date}. Actual: {basin_dates[0]}')
        if write_rerun_file: f.write(f'{basin_id}\n')

    # End dates: a file 1984-11-01 12:00 has values up until 1984-11-02 12:00 (1st timestep is at 13:00 on day 1)
    if basin_dates[-1]+timedelta(hours=24) < final_date: # because every file extends 12h into the next day
        print(f'    Last subset file for {basin_id} ends before expected date. Expected: {final_date}. Actual: {basin_dates[-1]}')
        if write_rerun_file: f.write(f'{basin_id}\n')

    # Check if basin dates are consecutive
    if not are_dates_consecutive(basin_dates):
        print(f'    RDRS dates for {basin_id} are not consecutive.')
        if write_rerun_file: f.write(f'{basin_id}\n')

if write_rerun_file: f.close()
    

