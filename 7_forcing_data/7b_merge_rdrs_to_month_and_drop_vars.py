import os
import sys
import xarray as xr
from pathlib import Path
from datetime import datetime
from collections import defaultdict
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

#### Command line arguments
if len(sys.argv) > 1:
    ix = int(sys.argv[1]) # 0, 1, 2, 3, 4. This will be the index of the year we process

#### Config handling
# Specify where the config file can be found
config_file = '../0_config/config.txt'

# Temporary download path
temp_folder = Path( cs.read_from_config(config_file, 'temp_path') )
temp_folder = temp_folder / 'rdrs' # Overwrite with actual location, so we can run a few things side-by-side
dest_folder = Path('/scratch/gwf/gwf_cmt/wknoben/') / 'rdrs' / 'monthly'
temp_folder.mkdir(parents=True, exist_ok=True)
dest_folder.mkdir(parents=True, exist_ok=True)

# Download URL
rdrs_vars = cs.read_from_config(config_file, 'rdrs_vars')
rdrs_vars_to_keep = rdrs_vars.split(',')

# main function
def merge_daily_to_monthly(input_dir, output_dir, variables_to_keep, month_index=None):
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Compression settings
    comp = dict(zlib=True, complevel=4)  # Use zlib compression with level 4

    # Gather all NetCDF files in the input directory
    files = [f for f in os.listdir(input_dir) if f.endswith('.nc')]
    files.sort()

    # Group files by month
    monthly_files = defaultdict(list)
    for file in files:
        # Extract date from filename
        date_str = file[:8]  # Assumes the format is yyyymmdd12.nc
        date = datetime.strptime(date_str, '%Y%m%d')
        month_key = date.strftime('%Y%m')
        monthly_files[month_key].append(file)
    
    # if year_index is provided, only process the months of that year
    if month_index is not None:
        months = list(monthly_files.keys())
        monthly_files = {months[month_index]: monthly_files[months[month_index]]}
        print(f"WARNING: Running reduced set. Processing year: {months[month_index]}")

    # Continue with regular processing, though monthly_files may contain only 1 year now
    for month, files in monthly_files.items():
        print(f"Processing month: {month}")
        output_file = os.path.join(output_dir, f"{month}.nc")

        # Skip if the file already exists
        if os.path.exists(output_file):
            print(f"File already exists: {output_file}")
            return # out of function

        # Process each file in the month
        datasets = []
        for file in files:
            filepath = os.path.join(input_dir, file)
            tmp_ds = xr.open_dataset(filepath)
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

        # Having merged the file, remove the individual daily files
        #for file in files:
        #    remove_this = os.path.join(input_dir, file)
        #    print(f'Removing {remove_this}')
        #    #os.remove(remove_this)
    return

#### Processing
merge_daily_to_monthly(temp_folder, dest_folder, rdrs_vars_to_keep, month_index=ix)        
        
