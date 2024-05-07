# Subset EM-Earth data to basins
import glob
import os
import sys
import pandas as pd
from datetime import datetime
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

# Temporary download path
temp_folder = Path( cs.read_from_config(config_file, 'temp_path') )
temp_folder = Path( '/gpfs/tp/gwf/gwf_cmt/wknoben/camels_spat/temp_em_earth' ) # Overwrite with actual location, so we can run a few things side-by-side

# --- Data loading
# CAMELS-spat metadata file
cs_meta_path = Path(data_path) / cs_meta_path
cs_meta = pd.read_csv(cs_meta_path / cs_meta_name)

# Open list of unusable stations; Enforce reading IDs as string to keep leading 0's
cs_unusable = pd.read_csv(cs_meta_path / cs_unusable_name, dtype={'Station_id': object})

# --- Processing
# Find the ERA5 files
em_earth_fold = temp_folder / 'EM_Earth_v1' / 'deterministic_hourly' / 'merged'
em_earth_files = sorted(glob.glob(str(em_earth_fold/'*.nc'))) # list

# Get the array index from the command-line argument
if len(sys.argv) != 2:
    print("Usage: python 2c_subset_EM_Earth_to_basins.py <array_index>")
    sys.exit(1)
else:
    ix = int(sys.argv[1])

# Define cs_meta row
row = cs_meta.iloc[ix] # needs to be between 0  and 1697

# Run
debug_message = f'\n!!! CHECK DEBUGGING STATUS: \n-Full array job run'
print(debug_message)

# Get shapefile path to determine download coordinates, and forcing destination path
basin_id, shp_lump_path, shp_dist_path, _, _ = cs.prepare_delineation_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder)
raw_fold, _, _ = cs.prepare_forcing_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder) # Returns folders only, not file names
print('--- Now running basin {}. {}'.format(ix, basin_id))

# From shapefile, get bounding coordinates. Then determine download coordinates from those
bounds = cs.find_shapefile_bounds(shp_lump_path)
coords_eme, _, _ = cs.find_download_coords_from_bounds(bounds, target='EM-Earth')

# Check if we need to run downloads for this station at all
missing = cs.flow_obs_unavailable(cs_unusable, row.Country, row.Station_id)
if 'iv' in missing and 'dv' in missing: 
    print(f'No flow observations for basin {basin_id}. Exiting.')
    sys.exit(0) # with next station, because we have no observations at all for this station. Error code 0: clean exit, no problems

# From meta-data, get download period
times_flow = cs.find_flow_obs_times_from_metadata(row, missing)
times_era5 = cs.round_flow_obs_to_days(times_flow)
start_date = datetime.strptime(times_era5[0], '%Y-%m-%d')
final_date = datetime.strptime(times_era5[1], '%Y-%m-%d')

print(f'    Basin coordinates:            {bounds}')
print(f'    EM-Earth subset coordinates: [{coords_eme}]')
print(f'    Flow obs unavailable:         {missing}')
print(f'    Download times:               {times_era5}')

# Convert start and end dates into two lists of start and end dates, that we'll iterate over
date_list,_ = cs.convert_start_and_end_dates_to_era5_download_lists(start_date,final_date) # not the cleanest but this lets us reuse old code
subset_strings = [date_obj.strftime("%Y-%m") for date_obj in date_list] # convert datetime objects to yyyy-mm strings

# Subset the data files
infiles = [file for file in em_earth_files if any(subset_string in file for subset_string in subset_strings)]

for infile in infiles:
    if os.path.exists(infile):
        file_name = os.path.basename(infile).replace('deterministic_hourly_NorthAmerica_','') # Make the name more similar to ERA5_YYYY-MM.nc
        outfile = raw_fold/file_name
        if not os.path.exists(outfile):
            cs.extract_ERA5_subset(infile,outfile,coords_eme)
    else:
        print(f'    ERROR: source file {infile} not found.')

# Create a figure to check if we actually cover the right domain with this
fig_file = raw_fold.parent / f'{row.Country}_{row.Station_id}_em_earth_coverage.png'
cs.compare_forcing_data_and_shape_extents(fig_file, outfile, shp_lump_path, nc_var='tmean', nc_time=0)

fig_file = Path('/gpfs/tp/gwf/gwf_cmt/wknoben/camels_spat/camels-spat-data/TEMP_images_forcing_coverage_em_earth') / f'{row.Country}_{row.Station_id}_em_earth_coverage.png'
cs.compare_forcing_data_and_shape_extents(fig_file, outfile, shp_lump_path, nc_var='tmean', nc_time=0)

print(debug_message)