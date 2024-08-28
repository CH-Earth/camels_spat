# Loops over basins and checks if we have all the RDRS data we expect.
#
# Expectations:
# - Files are stored in 'raw', 'lumped', and 'distributed' sub-folders
# - 12 files per year, per sub-folder
# - Files have no pre-determined start and end year
# - Dates within consecutive files must be consecutive
# - Each file must contain variables: 
#    'RDRS_v2.1_A_PR0_SFC', 'RDRS_v2.1_P_FB_SFC', 'RDRS_v2.1_P_FI_SFC',
#    'pet', 'RDRS_v2.1_P_TT_1.5m', 'RDRS_v2.1_P_P0_SFC', 'RDRS_v2.1_P_HU_1.5m',
#    'RDRS_v2.1_P_HR_1.5m', 'e', 'RDRS_v2.1_P_UUC_10m', 'RDRS_v2.1_P_VVC_10m',
#    'RDRS_v2.1_P_UVC_10m', 'phi', 'time', 'time_bnds'
# - Variables in 'raw' may contain NaN data, but not in 'lumped' or 'distributed'
# - Only 1 HRU in the 'lumped' case

import glob
import numpy as np
import re
import sys
import pandas as pd
import xarray as xr
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

# --- Basin selection
# Get the array index from the command-line argument
if len(sys.argv) != 2:
    print("Usage: python 7j_rdrs_checks.py <array_index>")
    sys.exit(1)
else:
    ix = int(sys.argv[1])

row = cs_meta.iloc[ix]

# Check if we need to run downloads for this station at all
missing = cs.flow_obs_unavailable(cs_unusable, row.Country, row.Station_id)
if 'iv' in missing and 'dv' in missing: 
    sys.exit(0) # because we have no observations at all for this station

# --- Data checks
# Get forcing paths
basin_id = row['Country'] + '_' + row['Station_id']
raw_fold, lump_fold, dist_fold = cs.prepare_forcing_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder) # Returns folders only, not file names

# Prepare two checking functions
def generate_expected_dates(start_date, end_date):
    date_range = pd.date_range(start=start_date, end=end_date, freq='MS') # months
    dates = [f"{date.strftime('%Y-%m')}" for date in date_range]
    return dates

def extract_dates_from_filenames(files):
    return [file[-10:-3] for file in files]

# Loop over the subfolders
for subfolder in [raw_fold, lump_fold, dist_fold]:

    # Find the files
    rdrs_files = glob.glob(str(subfolder / 'RDRS_*.nc'))
    rdrs_files.sort() # critical for later checks
    
    # - Begin checks
    # Do we have data at all?
    assert len(rdrs_files) > 0, f"Missing {subfolder} files for {basin_id}"

    # Do we have all months between first and last?
    # raw file example: RDRS_1980-01.nc
    # lumped file example: RDRS_lumped_remapped_RDRS_1980-01.nc
    # distributed file example: RDRS_dist_remapped_RDRS_1980-01.nc
    start_date = rdrs_files[0][-10:-3]
    final_date = rdrs_files[-1][-10:-3]
    want_dates = generate_expected_dates(start_date, final_date)
    have_dates = extract_dates_from_filenames(rdrs_files)
    assert want_dates == have_dates, f"Missing month(s) in {subfolder} files for {basin_id}"

    # Loop over the files for internal checks
    times = []
    for i, file in enumerate(rdrs_files):
        data = xr.open_dataset(file)

        # Check for the required variables
        for var in ['RDRS_v2.1_A_PR0_SFC', 'RDRS_v2.1_P_FB_SFC', 'RDRS_v2.1_P_FI_SFC',
                    'pet', 'RDRS_v2.1_P_TT_1.5m', 'RDRS_v2.1_P_P0_SFC', 'RDRS_v2.1_P_HU_1.5m',
                    'RDRS_v2.1_P_HR_1.5m', 'e', 'RDRS_v2.1_P_UUC_10m', 'RDRS_v2.1_P_VVC_10m',
                    'RDRS_v2.1_P_UVC_10m', 'phi', 'time', 'time_bnds']:
            assert var in data, f"Missing {var} in {file} for {basin_id}"

        # Check for NaNs in the lumped and distributed cases
        if subfolder == lump_fold or subfolder == dist_fold:
            for var in ['RDRS_v2.1_A_PR0_SFC', 'RDRS_v2.1_P_FB_SFC', 'RDRS_v2.1_P_FI_SFC',
                    'pet', 'RDRS_v2.1_P_TT_1.5m', 'RDRS_v2.1_P_P0_SFC', 'RDRS_v2.1_P_HU_1.5m',
                    'RDRS_v2.1_P_HR_1.5m', 'e', 'RDRS_v2.1_P_UUC_10m', 'RDRS_v2.1_P_VVC_10m',
                    'RDRS_v2.1_P_UVC_10m', 'phi']:
                assert data[var].isnull().sum() == 0, f"NaNs in {var} in {file} for {basin_id}"

        # Check for only 1 HRU in the lumped case
        if subfolder == lump_fold:
            assert data['hru'].shape[0] == 1, f"More than 1 HRU in {file} for {basin_id}"

        # Check for consecutive dates by concatenating the dates from all files
        times.append(data['time'].values) # list of arrays of datetimes

        # Close the file
        data.close()
    
    # Check for consecutive dates
    # This implicitly also checks if we added the missing days in leap years correctly
    datetimes = np.concatenate(times) # array of datetimes
    expected_datetimes = pd.date_range(start=datetimes[0], end=datetimes[-1], freq='H')
    assert (datetimes == expected_datetimes).all(), f"Non-consecutive dates in {subfolder} files for {basin_id}"