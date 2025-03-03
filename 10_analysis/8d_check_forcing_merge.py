import glob
import numpy as np
import os
import re
import sys
import time
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

# Base folder
base_folder = Path('/scratch/gwf/gwf_cmt/wknoben/camels-spat-upload')

# 2. Meta-data subsetting to 1426 final basins
cs_meta_upload = cs_meta[~cs_meta.set_index(['Country', 'Station_id']).index.isin(cs_unusable.set_index(['Country', 'Station_id']).index)]

# -- ARGUMENTS
if len(sys.argv) != 2:
    print("Usage: python 8d_check_forcing_merge.py <array_index>")
    sys.exit(1)
else:
    ix = int(sys.argv[1])

# Get basic info
row = cs_meta_upload.iloc[ix]
basin_id = row['Country'] + '_' + row['Station_id']
category = row['subset_category']

# Generate the loops
forcing_data = ['daymet','em-earth','era5','rdrs']
spatial_sets = ['distributed','gridded','lumped']

# Check the files
print(f"{basin_id}")
for f_data in forcing_data:
    for spatial_set in spatial_sets:
        
        # Construct paths to data
        folder = base_folder / 'forcing' / category / f_data / f"{f_data}-{spatial_set}"
        filename = f"{basin_id}_{f_data.replace('-','_')}_{spatial_set}.nc"
        file = folder/filename
        
        # 1. Ensure file exists
        if not os.path.exists(file):
            print(f'file not found: {file}')
        else:
            data = xr.open_dataset(file)

        # 2. Run the variable checks
        if f_data == 'daymet':
            
            # Check for the required variables
            for var in ['prcp', 'tmax', 'tmin', 'srad', 'vp', 'dayl', 'pet', 'time', 'time_bnds']:
                assert var in data, f"Missing {var} in {file} for {basin_id}"
    
            # Check for NaNs in the lumped and distributed cases
            if spatial_set == 'lumped' or spatial_set == 'distributed':
                for var in ['prcp', 'tmax', 'tmin', 'srad', 'vp', 'dayl', 'pet']:
                    assert data[var].isnull().sum() == 0, f"NaNs in {var} in {file} for {basin_id}"
    
            # Check for only 1 HRU in the lumped case
            if spatial_set == 'lumped':
                assert data['hru'].shape[0] == 1, f"More than 1 HRU in {file} for {basin_id}"

        if f_data == 'em-earth':

            # Check for the required variables
            for var in ['prcp','tmean','time', 'time_bnds']:
                assert var in data, f"Missing {var} in {file} for {basin_id}"
    
            # Check for NaNs in the lumped and distributed cases       
            if spatial_set == 'lumped' or spatial_set == 'distributed':
                for var in ['prcp','tmean']:
                    # handle the special case missing prcp values in some cases
                    if var == 'prcp':
                        assert data['prcp'].isnull().sum() / data['prcp'].notnull().sum() < 0.0001, f"NaNs in {var} in {file} for {basin_id}" # skips the first 7 time steps
                    else:
                        assert data[var].isnull().sum() == 0, f"NaNs in {var} in {file} for {basin_id}"
    
            # Check for only 1 HRU in the lumped case
            if spatial_set == 'lumped':
                assert data['hru'].shape[0] == 1, f"More than 1 HRU in {file} for {basin_id}"

        if f_data == 'era5':

            # Check for the required variables
            for var in ['mtpr','msdwswrf','msdwlwrf','msnswrf','msnlwrf','mper',
                        't','sp','q','rh','e','u','v','w','phi',
                        'time', 'time_bnds']:
                assert var in data, f"Missing {var} in {file} for {basin_id}"
    
            # Check for NaNs in the lumped and distributed cases
            if spatial_set == 'lumped' or spatial_set == 'distributed':
                for var in ['mtpr','msdwswrf','msdwlwrf','msnswrf','msnlwrf','mper',
                            't','sp','q','rh','e','u','v','w','phi']:
                    assert data[var].isnull().sum() == 0, f"NaNs in {var} in {file} for {basin_id}"
    
            # Check for only 1 HRU in the lumped case
            if spatial_set == 'lumped':
                assert data['hru'].shape[0] == 1, f"More than 1 HRU in {file} for {basin_id}"

        if f_data == 'rdrs':

            # Check for the required variables
            for var in ['RDRS_v2.1_A_PR0_SFC', 'RDRS_v2.1_P_FB_SFC', 'RDRS_v2.1_P_FI_SFC',
                        'pet', 'RDRS_v2.1_P_TT_1.5m', 'RDRS_v2.1_P_P0_SFC', 'RDRS_v2.1_P_HU_1.5m',
                        'RDRS_v2.1_P_HR_1.5m', 'e', 'RDRS_v2.1_P_UUC_10m', 'RDRS_v2.1_P_VVC_10m',
                        'RDRS_v2.1_P_UVC_10m', 'phi', 'time', 'time_bnds']:
                assert var in data, f"Missing {var} in {file} for {basin_id}"
    
            # Check for NaNs in the lumped and distributed cases
            if spatial_set == 'lumped' or spatial_set == 'distributed':
                for var in ['RDRS_v2.1_A_PR0_SFC', 'RDRS_v2.1_P_FB_SFC', 'RDRS_v2.1_P_FI_SFC',
                        'pet', 'RDRS_v2.1_P_TT_1.5m', 'RDRS_v2.1_P_P0_SFC', 'RDRS_v2.1_P_HU_1.5m',
                        'RDRS_v2.1_P_HR_1.5m', 'e', 'RDRS_v2.1_P_UUC_10m', 'RDRS_v2.1_P_VVC_10m',
                        'RDRS_v2.1_P_UVC_10m', 'phi']:
                    assert data[var].isnull().sum() == 0, f"NaNs in {var} in {file} for {basin_id}"
    
            # Check for only 1 HRU in the lumped case
            if spatial_set == 'lumped':
                assert data['hru'].shape[0] == 1, f"More than 1 HRU in {file} for {basin_id}"

        # 3. Check the time period
        years = len(pd.unique(data['time'].dt.year))
        print(f"{f_data.ljust(8)} {spatial_set.ljust(12)}: {years} years")
        
        data.close()