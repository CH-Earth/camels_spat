from collections import defaultdict
from datetime import datetime
import os
import pandas as pd
from pathlib import Path
import sys
import xarray as xr
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

# --- Command line arguments
if len(sys.argv) != 2:
    print("Usage: python 7c_merge_rdrs_to_month.py <array_index>")
    sys.exit(1)
else:
    ix = int(sys.argv[1]) # Basin index in cs_meta

# -- Set the overwrite flag
overWrite = True

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

# Define where the RDRS folder with daily data is
raw_fold, _, _ = cs.prepare_forcing_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder) # Returns folders only, not file names
temp_day_folder = raw_fold / 'rdrs_day'
temp_month_folder = raw_fold / 'rdrs_month'
temp_month_folder.mkdir(exist_ok=True, parents=True)

# Find the yearly folders in the RDRS temp directory
year_folders = [f.name for f in temp_day_folder.iterdir() if f.is_dir()]
year_folders.sort() # Doesn't matter but feels cleaner

# main function
def merge_daily_to_monthly(input_dir, output_dir, variables_to_keep=None):
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Compression settings
    comp = dict(zlib=True, complevel=4)  # Use zlib compression with level 4

    # Gather all NetCDF files in the input directory
    files = [f for f in os.listdir(input_dir) if f.endswith('12.nc')]
    files.sort()

    # Group files by month
    monthly_files = defaultdict(list)
    for file in files:
        # Extract date from filename
        date_str = file.split('_')[-1].split('.')[0][:8]  # Assumes the format is country_basin_subset_yyyymmdd12.nc
        date = datetime.strptime(date_str, '%Y%m%d')
        month_key = date.strftime('%Y-%m')
        monthly_files[month_key].append(file)

    # Loop over the months and concatenate
    for yyyymm, files in monthly_files.items():
        output_file = os.path.join(output_dir, f'RDRS_{yyyymm}.nc')
        print(f"Processing month: {yyyymm} into {output_file}")

        # Skip if the file already exists
        if os.path.exists(output_file) and not overWrite:
            print(f"File already exists: {output_file}")
            return # out of function

        # Process each file in the month
        datasets = []
        for file in files:
            filepath = os.path.join(input_dir, file)
            tmp_ds = xr.open_dataset(filepath)
            if variables_to_keep is not None:
                variables_to_drop = [var for var in tmp_ds.variables if (var not in variables_to_keep and (var not in tmp_ds.coords)) and (var != 'rotated_pole')]
                tmp_ds = tmp_ds.drop_vars(variables_to_drop)
            datasets.append(tmp_ds)

        # Merge the datasets
        ds = xr.concat(datasets, dim='time')
        encoding = {var: comp for var in ds.data_vars}
        ds.to_netcdf(output_file, mode='w', encoding=encoding)
        ds.close()
        print(f"Saved monthly file: {output_file}")

        # Close the datasets   
        for ds in datasets:
            ds.close()

    return

# Loop over the years, and the months within each year, then merge directly into raw_folder
for year in year_folders:
    input_dir = temp_day_folder / year 
    merge_daily_to_monthly(input_dir, temp_month_folder)        
        