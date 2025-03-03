import os
import pandas as pd
from pathlib import Path
import shutil
import sys
import xarray as xr
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

# --- Command line arguments
if len(sys.argv) != 2:
    print("Usage: python 6e_merge_daymet_vars_to_single_year <array_index>")
    sys.exit(1)
else:
    ix = int(sys.argv[1]) # Basin index in cs_meta

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

# --- Processing
# Define cs_meta row
row = cs_meta.iloc[ix] # needs to be between 0  and 1697
basin_id = row['Country'] + '_' + row['Station_id']

# Check if we need to run for this station at all
missing = cs.flow_obs_unavailable(cs_unusable, row.Country, row.Station_id)
if 'iv' in missing and 'dv' in missing: 
    print(f'No flow observations for basin {basin_id}. Exiting.')
    sys.exit(0) # with next station, because we have no observations at all for this station. Error code 0: clean exit, no problems
print('--- Now running basin {}. {}'.format(ix, basin_id))

# Define where the Daymet folder is
raw_fold, _, _ = cs.prepare_forcing_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder) # Returns folders only, not file names
temp_daymet_folder = raw_fold / 'daymet_year'

# Find the files in the Daymet folder
daymet_files = list(temp_daymet_folder.glob('*_subset_daymet_v4_daily_na_*.nc')) 
daymet_files.sort() # Doesn't matter but feels cleaner

# Loop over the PET files and ensure lat/lon are coordinates
for file in daymet_files:
    if 'pet' in file.name:
        tmp_file = str(file).replace('.nc', '.in.nc')
        shutil.copy(file, tmp_file)
        ds = xr.open_dataset(tmp_file)
        ds = ds.set_coords(('lat','lon'))
        ds.to_netcdf(file)
        os.remove(tmp_file)

# Ensure we don't have any '*.in.nc' files left
in_files = list(temp_daymet_folder.glob('*.in.nc'))
for file in in_files:
    os.remove(file)

# Extract the years stored in the daymet file paths
years = []
for file in daymet_files:
    year = int(file.name.split('_')[-1].split('.')[0])
    years.append(year)

# Get the unique years
years = list(set(years))

# Loop over years and create a single file for each year
for year in years:
    # Subset the file list by year
    ds_list = [xr.open_dataset(f) for f in daymet_files if f'_{year}.nc' in str(f)]
    # For some reason xr.concat() creates a time dimension of length 2555 (365*7 variables)
    # This happens with xr.concat(ds_list, dim='time', coords='minimal', compat='override')
    # and xr.concat(ds_list, dim='time'). xr.merge has the same problem so we need a workaround.
    data = None
    for ds in ds_list:
        if data is None:
            data = ds
        else:
            for var in ds.data_vars:
                if var in ['dayl', 'prcp', 'srad', 'tmax', 'tmin', 'vp', 'pet']:
                    #print(f'Adding {var}')
                    data[var] = ds[var] 
    
    data.to_netcdf(temp_daymet_folder / f'daymet_{year}_stillLeapYearGaps.nc')


