# Loops over basins and checks if we have all the Daymet data we expect.
#
# Expectations:
# - Files are stored in 'raw', 'lumped', and 'distributed' sub-folders
# - One file per year, per sub-folder
# - Files have no pre-determined start and end year
# - Dates within consecutive files must be consecutive
# - Each file must contain variables: 'prcp', 'tmax', 'tmin', 'srad', 'vp', 'dayl', 'pet'
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
    print("Usage: python 6i_daymet_checks.py <array_index>")
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

# Loop over the subfolders
for subfolder in [raw_fold, lump_fold, dist_fold]:

    # Find the files
    daymet_files = glob.glob(str(subfolder / 'daymet_*.nc'))
    daymet_files.sort() # critical for later checks
    
    # - Begin checks
    # Do we have data at all?
    assert len(daymet_files) > 0, f"Missing {subfolder} files for {basin_id}"

    # Do we have all years between first and last?
    # raw file example: daymet_1980.nc
    # lumped file example: Daymet_lumped_remapped_1980-01-01-12-00-00.nc
    # distributed file example: Daymet_dist_remapped_1980-01-01-12-00-00.nc
    years = []
    for file in daymet_files:
        name = Path(file).name
        years.append(int(re.findall(r'\d+', name)[0]))
        assert (np.diff(np.array(years)) == 1).all(), f"Missing year(s) in {subfolder} files for {basin_id}"

    # Loop over the files for internal checks
    times = []
    for i, file in enumerate(daymet_files):
        data = xr.open_dataset(file)

        # Check for the required variables
        for var in ['prcp', 'tmax', 'tmin', 'srad', 'vp', 'dayl', 'pet', 'time', 'time_bnds']:
            assert var in data, f"Missing {var} in {file} for {basin_id}"

        # Check for NaNs in the lumped and distributed cases
        if subfolder == lump_fold or subfolder == dist_fold:
            for var in ['prcp', 'tmax', 'tmin', 'srad', 'vp', 'dayl', 'pet']:
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
    expected_datetimes = pd.date_range(start=datetimes[0], end=datetimes[-1], freq='D')
    assert (datetimes == expected_datetimes).all(), f"Non-consecutive dates in {subfolder} files for {basin_id}"

