# Loop over the Daymet files in raw, as well as lumped 
# and dsitributed, and copy the time_bnds varaible from 
# raw over to the other two.
import glob
import sys
import pandas as pd
import re
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

# --- Select basin
# Get the array index from the command-line argument
if len(sys.argv) != 2:
    print("Usage: python 6i_copy_time_bnds_to_lump_and_dist.py <array_index>")
    sys.exit(1)
else:
    ix = int(sys.argv[1])
 
# Define cs_meta row
row = cs_meta.iloc[ix] # needs to be between 0  and 1697

# Check if we need to run downloads for this station at all
missing = cs.flow_obs_unavailable(cs_unusable, row.Country, row.Station_id)
if 'iv' in missing and 'dv' in missing: 
    sys.exit(0) # exit gracefully, because we have no observations at all for this station

# --- Processing
# Get shapefile path to determine download coordinates, and forcing destination path
basin_id = row['Country'] + '_' + row['Station_id']
raw_fold, lump_fold, dist_fold = cs.prepare_forcing_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder) # Returns folders only, not file names
print('--- Now running basin {}. {}'.format(ix, basin_id))

# Find the daymet files
raw_files = sorted(glob.glob(str(raw_fold/'daymet_*.nc')))
lum_files = sorted(glob.glob(str(lump_fold/'daymet_*.nc')))
dis_files = sorted(glob.glob(str(dist_fold/'daymet_*.nc')))

# Loop over the files
for raw_f, lum_f, dis_f in zip(raw_files, lum_files, dis_files):
    # Check that we're looking at the same years
    raw_year = int(re.findall(r'\d+', Path(raw_f).name)[0])
    lum_year = int(re.findall(r'\d+', Path(lum_f).name)[0])
    dis_year = int(re.findall(r'\d+', Path(dis_f).name)[0])
    assert raw_year == lum_year == dis_year, 'Years do not match'
    
    # Open the datasets
    raw_ds = xr.open_dataset(raw_f, mode='a')
    lum_ds = xr.open_dataset(lum_f, mode='a')
    dis_ds = xr.open_dataset(dis_f, mode='a')

    # Copy the time_bnds variable
    lum_ds['time_bnds'] = raw_ds['time_bnds']
    dis_ds['time_bnds'] = raw_ds['time_bnds']

    # Drop the empty License attribute while we're here
    del lum_ds.attrs['License']
    del dis_ds.attrs['License']

    # To file
    lum_ds.to_netcdf(lum_f, mode='a')
    dis_ds.to_netcdf(dis_f, mode='a')

    # Close the datasets
    raw_ds.close()
    lum_ds.close()
    dis_ds.close()
