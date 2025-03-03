# Checks a single basin for forcing file consistency
# Saves findings to a csv file
from datetime import datetime, timedelta
import glob
import numpy as np
import os
import pandas as pd
import sys
import time
import xarray as xr
from pathlib import Path

# Start timer
start = time.time()

# Define csv location
rerun_folder = '20240522_checks'
rerun_path = Path('/project/gwf/gwf_cmt/wknoben/camels_spat/camels-spat-data/forcing_reruns') / rerun_folder
rerun_path.mkdir(parents=True, exist_ok=True)

## Functions
def list_folders(directory):
    folders = []
    directory = str(directory)
    for item in os.listdir(directory):
        if os.path.isdir(os.path.join(directory, item)):
            folders.append(item)
    return sorted(folders)

def list_forcing_files(folder,data):
    if (data != 'EM_Earth') and (data != 'ERA5'):
        print(f'Error: Incorrect data: {data}. Must be "EM_Earth" or "ERA5"')
        return
    search_string = str(folder) + f'/{data}_*.nc'
    files = sorted(glob.glob(search_string))
    return files

def convert_filenames_to_yyyymm_strings(files):
    '''Extracts yyyy-mm values from file string and compares'''
    year_months = []
    for file in files:
        file_name = os.path.basename(file) # "<data>_yyyy-mm.nc" or "<data>_<spatial>_remapped_yyyy-mm-01-00-00-00.nc"
        file_parts = file_name.split('_') # final part should always be "yyyy-mm[..].nc"
        year_months.append(file_parts[-1][0:7])
    return year_months

def are_consecutive_dates(date_list, debug=False):
    for i in range(len(date_list) - 1):
        # Get current yyyy-mm
        this = datetime.strptime(date_list[i], "%Y-%m")
        
        # Define the expected next yyyy-mm value
        next_month = this.month + 1
        if next_month == 13: 
            next_month = 1
            next_year = this.year+1
        else:
            next_year = this.year
        next_string = f'{next_year}-{next_month:02}'

        # Compare expected next with real next
        real_next = datetime.strptime(date_list[i+1], "%Y-%m")
        want_next = datetime.strptime(next_string, "%Y-%m")

        if debug:
            print(f'Current:  {this}')
            print(f'Expected: {want_next}')
            print(f'Actual  : {real_next}\n')

        if real_next != want_next:
            print(f'Error: Non-consecutive dates: {want_next} missing until {real_next}')
            return False # Stop checks, because one missing month is one too many
    return True

def has_variables_in_files(files,data):
    '''Checks if the forcing files contain the expected variables'''
    
    # Variables assumed to be in files
    if data == 'ERA5':
        expected = ['mtpr', 'msdwswrf', 'msdwlwrf',
                    'msnswrf', 'msnlwrf', 'mper', 
                    't', 'u', 'v', 'q',
                    'e', 'rh', 'w', 'phi', 'time_bnds']
    elif data == 'EM_Earth':
        expected = ['prcp', 'tmean', 'time_bnds']

    # Initialize flags: we only need to find one missing variable to fail
    # Print statements will give further details in case of multiple files
    # missing variables
    has_vars = True
    has_nans = False

    # Run checks
    for file in files:
        data = xr.open_dataset(file)
        for var in expected:
            if var not in data:
                print(f'Error: Missing variable {var} in {file}')
                has_vars = False
            else: # This has to be here to prevent checking a non-existent variable for NaNs
                if np.isnan(data[var].values).any():
                    print(f'Error: NaNs in variable {var} in {file}')
                    has_nans = True
    return has_vars, has_nans

def check_if_list_contains_invariant(files):
    '''Checks if an ERA5 invariants file is present'''
    return any('invariant' in file for file in files)

def remove_invariants_from_list(files):
    '''Checks if an ERA5 invariants file is present and removes if so'''
    return [file for file in files if 'invariant' not in file]

def check_if_times_in_lst(files):
    '''Checks if the time variable is in local standard time, under the assumption that the first time entry in each file is not at midnight'''
    for file in files:
        data = xr.open_dataset(file)
        time = data['time'].values[0] # e.g. numpy.datetime64('1986-07-31T20:00:00.000000000')
        hour = time.astype('datetime64[h]').astype(int) % 24 # remainder after division, i.e. 0 if midnight
        if hour == 0:
            print(f'Error: Time series starts at midnight in {file}. This implies conversion to LST is not done.')
            return False
    return True

# Main functions
def check_forcing_file_status(folder, data):
    files = list_forcing_files(folder,data)
    files = remove_invariants_from_list(files)
    dates = convert_filenames_to_yyyymm_strings(files)
    is_consecutive = are_consecutive_dates(dates, debug=False)
    has_variables, has_nans = has_variables_in_files(files, data)
    in_lst = check_if_times_in_lst(files)
    return is_consecutive,has_variables, has_nans, in_lst

def check_invariant_status(folder, data):
    files = list_forcing_files(folder,data)
    has_invariant = check_if_list_contains_invariant(files)
    return has_invariant

## Main
reruns = []
reasons = []

# Workflow - data complete?
# 1. Define the data location
# 2. Loop over all basins
# 3. Loop over the three spatial discretization folders
# 4. Loop over the two different data sets
# 5. Check if we have an invariant file for the basin
# 6. Check if the forcing files are consecutive
# 7. Check if each file contains the expected variables without NaNs
main_folder = Path('/project/gwf/gwf_cmt/wknoben/camels_spat/camels-spat-data/basin_data')
basin_folders = list_folders(main_folder) # sorted listed so always same order

# Get command line input and store as ix
ix = int(sys.argv[1])
basin = basin_folders[ix] # e.g. 'CAN_09ED001'
print(f'Checking basin forcing data for {basin}')

for basin_folder in basin_folders:

    # Skip to next if we haven't selected the current basin
    if basin_folder != basin:
        continue

    for spatial_res in ['raw', 'lumped', 'distributed']:
        for dataset in ['ERA5','EM_Earth']:
            # Check file status
            data_folder = main_folder / basin_folder / 'forcing' / spatial_res
            is_consecutive,has_variables, has_nans, times_in_lst = check_forcing_file_status(data_folder, dataset)
            if not is_consecutive:
                reruns.append(basin_folder)
                reasons.append(f'Non-consecutive dates in {spatial_res} {dataset}')
            if not has_variables:
                reruns.append(basin_folder)
                reasons.append(f'Missing variables in {spatial_res} {dataset}')
            if has_nans:
                reruns.append(basin_folder)
                reasons.append(f'NaNs in variables in {spatial_res} {dataset}')
            if not times_in_lst:
                reruns.append(basin_folder)
                reasons.append(f'Times may not be converted to LST in {spatial_res} {dataset}')
            # Check invariant status
            if spatial_res == 'raw' and dataset == 'ERA5':
                has_invariant = check_invariant_status(data_folder, dataset)
                if not has_invariant:
                    reruns.append(basin_folder)
                    reasons.append('Missing invariant file')
    print(f'Finished checking {basin_folder}')
    print(f'Reruns: {reruns}')
    print(f'Reasons: {reasons}')

##### SPECIFIC CHECKS 2024-05-16 #####
# We can disable these now (2024-05-22), because we should have addressed the two issues below after out 2024-05-17 to 2024-05-22 reruns.
# ---------------------------------------------------------------------------------------------------------------------------------------      
# # Workflow - time zones is NST?
# # 1. Open the metadata file
# # 2. Check if the time zone is NST
# meta_path = Path('/project/gwf/gwf_cmt/wknoben/camels_spat/camels-spat-data/camels_spat_metadata.csv')
# meta = pd.read_csv(meta_path)
# mask = (meta['Country'] == basin[0:3]) & (meta['Station_id'] == basin[4:])
# if len(meta.loc[mask]) == 0:
#     print(f'Error: Basin {basin} not found in metadata')
# else:
#     if meta.loc[mask,'dv_flow_obs_timezone'].values[0] == 'NST':
#         reruns.append(basin)
#         reasons.append('Time zone is NST')

# # Workflow - check if we had to remove tiny HRUs from this basin
# # We updated these after manual checks in: ./5_forcing/5_fix_minor_errors.py.
# # We can ignore the 'double_comid' issue, because that relates to duplicated 
# # river segments, not HRUs. issue_dict coped directly from the file above, 
# # and commented the basins with ONLY 'double_comid' issues.
# issue_dict = {'CAN_02GG003': ['hru_area'], 
#               'CAN_04FC001': ['hru_area', 'double_comid'], 
#               'CAN_07AD002': ['hru_area'], 
#               'CAN_07HC002': ['hru_area'], 
#               'CAN_08NE077': ['hru_area'], 
#               'CAN_08NH007': ['hru_area'], 
#               'USA_07142300':['hru_area'], 
#               'USA_08198500':['hru_area'],
#               #'CAN_06AB001': ['double_comid'],
#               #'CAN_06BB005': ['double_comid'],
#               #'CAN_08MH076': ['double_comid']
#               }

# if basin in issue_dict.keys():
#     reruns.append(basin)
#     reasons.append('Basin had tiny HRUs removed')
# End commenting this out for 2024-05-22 check of re-runs.

# Create a dataframe for storage
out = pd.DataFrame(data={'basin':reruns, 'reason':reasons})
out.to_csv(rerun_path / f'{basin}_reruns.csv', index=False)

# End timer
end = time.time()
print(f'Elapsed time: {end-start} seconds')