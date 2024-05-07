from datetime import datetime, timedelta
import glob
import numpy as np
import os
from pathlib import Path

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
            return False # Stop checks, because one missing month is one too many
    return True

def remove_invariants_from_list(files):
    '''Checks if an ERA5 invariants file is present and removes if so'''
    return [file for file in files if 'invariant' not in file]

# main function
def check_forcing_file_names(folder, data):
    files = list_forcing_files(folder,data)
    files = remove_invariants_from_list(files)
    dates = convert_filenames_to_yyyymm_strings(files)
    is_consecutive = are_consecutive_dates(dates, debug=False)
    return is_consecutive

# processing
main_folder = Path('/project/gwf/gwf_cmt/wknoben/camels_spat/camels-spat-data/basin_data')
sub_folders = list_folders(main_folder)
for sub_folder in sub_folders:
    for dataset in ['ERA5','EM_Earth']:
        for spatial_res in ['raw', 'lumped', 'distributed']:
            data_folder = main_folder / sub_folder / 'forcing' / spatial_res
            is_consecutive = check_forcing_file_names(data_folder, dataset)   
            if not is_consecutive:
                print(f'Some monthly files are missing for {sub_folder}, at resolution {spatial_res:<11}, for dataset {dataset}')
