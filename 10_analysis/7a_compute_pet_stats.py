# Calculate basin attributes
# Accepts basin ID as command line input

import glob
import xarray as xr
import numpy as np
import os
import pandas as pd
import shutil
import sys
from pathlib import Path
sys.path.append(str(Path().absolute().parent))
from python_cs_functions import config as cs, attributes as csa
from python_cs_functions.delineate import prepare_delineation_outputs

# --- Config handling
# Specify where the config file can be found
config_file = '../0_config/config.txt'

# Get the required info from the config file
data_path            = cs.read_from_config(config_file,'data_path')

# CAMELS-spat metadata
cs_meta_path = cs.read_from_config(config_file,'cs_basin_path')
cs_meta_name = cs.read_from_config(config_file,'cs_meta_name')
cs_unusable_name = cs.read_from_config(config_file,'cs_unusable_name')

# Basin folder
cs_basin_folder = cs.read_from_config(config_file, 'cs_basin_path')
basins_path = Path(data_path) / cs_basin_folder

# Get the attribute folder
att_folder = cs.read_from_config(config_file, 'att_path')
att_path = basins_path / att_folder / 'PET_analysis'
att_path.mkdir(parents=True, exist_ok=True)

# --- Data loading
# CAMELS-spat metadata file
cs_meta_path = Path(data_path) / cs_meta_path
cs_meta = pd.read_csv(cs_meta_path / cs_meta_name)

# Open list of unusable stations; Enforce reading IDs as string to keep leading 0's
cs_unusable = pd.read_csv(cs_meta_path / cs_unusable_name, dtype={'Station_id': object})

# --- Command line arguments
# Get the array index from the command-line argument
if len(sys.argv) != 2:
    print("Usage: python 1_calculate_attributes.py <array_index>")
    sys.exit(1)
else:
    ix = int(sys.argv[1])

# Get metadata
row = cs_meta.iloc[ix]

# Skip the basins for which we don't have daily flow obs
if (row['Station_id'] == cs_unusable['Station_id']).any():
    ixs = np.where(row['Station_id'] == cs_unusable['Station_id'])
    unusable = cs_unusable.iloc[ixs]    
    for i,r in unusable.iterrows():
        if r['Missing'] == 'dv':
            print('no daily streamflow observations available for this basin.')
            sys.exit(0)

# --- Processing
# Define various conversion constants
water_density = 1000 # kg m-3
mm_per_m = 1000 # mm m-1
seconds_per_hour = 60*60 # s h-1
seconds_per_day = seconds_per_hour*24 # s d-1
days_per_month = np.array([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]).reshape(-1, 1).flatten()
days_per_year = days_per_month.sum()

# Storage lists
l_rdrs_longterm_mean_p = []
l_rdrs_longterm_mean_pet = []
l_rdrs_annual_mean_p = []
l_rdrs_annual_mean_pet = []
l_daym_longterm_mean_p = []
l_daym_longterm_mean_pet = []
l_daym_annual_mean_p = []
l_daym_annual_mean_pet = []
l_era5_longterm_mean_p = []
l_era5_longterm_mean_pet = []
l_era5_annual_mean_p = []
l_era5_annual_mean_pet = []

# Get the paths
basin_id, _, _, _, _ = prepare_delineation_outputs(cs_meta, ix, basins_path)
met_folder = basins_path / 'basin_data' / basin_id / 'forcing'

# Get the files
era5_files = glob.glob(str(met_folder / 'lumped' / 'ERA5_*.nc'))
rdrs_files = glob.glob(str(met_folder / 'lumped' / 'RDRS_*.nc'))
daymet_files = glob.glob(str(met_folder / 'lumped' / 'daymet_*.nc'))

# Open the files with xarray
ds_era5 = xr.open_mfdataset(era5_files, combine="by_coords")
ds_rdrs = xr.open_mfdataset(rdrs_files, combine="by_coords") # Don't use 'decode_cf=False' > this somehow loses most of the timesteps
ds_daymet = xr.open_mfdataset(daymet_files, combine="by_coords")

# Calculate long-term stats RDRS
l_rdrs_longterm_mean_p.append((  ds_rdrs['RDRS_v2.1_A_PR0_SFC'].mean().compute() * seconds_per_day * mm_per_m / water_density).values) # mm day-1
l_rdrs_longterm_mean_pet.append((ds_rdrs['pet'].mean().compute() * seconds_per_day * mm_per_m / water_density).values)
l_rdrs_annual_mean_p.append(    (ds_rdrs['RDRS_v2.1_A_PR0_SFC'].resample(time='1Y').mean().compute() * seconds_per_day * mm_per_m / water_density).mean().values)
l_rdrs_annual_mean_pet.append(  (ds_rdrs['pet'].resample(time='1Y').mean().compute() * seconds_per_day * mm_per_m / water_density).mean().values)

# Calculate Daymet-equivalent of the RDRS mean-mean values
l_daym_longterm_mean_p.append(  ds_daymet['prcp'].mean().compute().values) # mm d-1
l_daym_longterm_mean_pet.append(ds_daymet['pet'].mean().compute().values)
l_daym_annual_mean_p.append((   ds_daymet['prcp'].resample(time='1Y').mean().compute()).mean().values)
l_daym_annual_mean_pet.append(( ds_daymet['pet'].resample(time='1Y').mean().compute()).mean().values)

# ERA5
l_era5_longterm_mean_p.append(  (ds_era5['mtpr'].mean().compute() * seconds_per_day * mm_per_m / water_density).values)
l_era5_longterm_mean_pet.append((ds_era5['mper'].mean().compute() * seconds_per_day * mm_per_m / water_density).values)
l_era5_annual_mean_p.append(    (ds_era5['mtpr'].resample(time='1Y').mean().compute() * seconds_per_day * mm_per_m / water_density).mean().values)
l_era5_annual_mean_pet.append(  (ds_era5['mper'].resample(time='1Y').mean().compute() * seconds_per_day * mm_per_m / water_density).mean().values)


# --- Make the dataframe
df = pd.DataFrame(data={'basin_id': basin_id,
                   'rdrs_long_mean_p':     l_rdrs_longterm_mean_p,
                   'rdrs_long_mean_e':     l_rdrs_longterm_mean_pet,
                   'rdrs_annual_mean_p':   l_rdrs_annual_mean_p,
                   'rdrs_annual_mean_e':   l_rdrs_annual_mean_pet,
                   'daymet_long_mean_p':   l_daym_longterm_mean_p,
                   'daymet_long_mean_e':   l_daym_longterm_mean_pet,
                   'daymet_annual_mean_p': l_daym_annual_mean_p,
                   'daymet_annual_mean_e': l_daym_annual_mean_pet,
                   'era5_long_mean_p':     l_era5_longterm_mean_p,
                   'era5_long_mean_e':     l_era5_longterm_mean_pet,
                   'era5_annual_mean_p':   l_era5_annual_mean_p,
                   'era5_annual_mean_e':   l_era5_annual_mean_pet})

# Save to file
att_file = f'pet_stats_{basin_id}.csv'
df.to_csv(att_path/att_file, index=False)
        