# Add wind direction
import glob
import sys
import netCDF4 as nc4
import pandas as pd
from pathlib import Path
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

## Indices for which we need to rerun
ixs = [963,986] # Only 2 basins that were incomplete
last_completes =['/gpfs/tp/gwf/gwf_cmt/wknoben/camels_spat/camels-spat-data/basin_data/CAN_09ED001/forcing/distributed/ERA5_dist_remapped_1986-10-01-00-00-00.nc','/gpfs/tp/gwf/gwf_cmt/wknoben/camels_spat/camels-spat-data/basin_data/CAN_10ED002/forcing/distributed/ERA5_dist_remapped_1977-08-01-00-00-00.nc']

## Config handling
# Specify where the config file can be found
config_file = '../0_config/config.txt'

# Get the required info from the config file
data_path = cs.read_from_config(config_file,'data_path')

# CAMELS-spat metadata
cs_meta_path = cs.read_from_config(config_file,'cs_basin_path')
cs_meta_name = cs.read_from_config(config_file,'cs_meta_name')

# Basin folder
cs_basin_folder = cs.read_from_config(config_file, 'cs_basin_path')
basins_path = Path(data_path) / cs_basin_folder

## Data loading
# CAMELS-spat metadata file
cs_meta_path = Path(data_path) / cs_meta_path
cs_meta = pd.read_csv(cs_meta_path / cs_meta_name)

## Build the list of files we need to process still
file_list = []
for ix,last_complete in zip(ixs, last_completes):
    row = cs_meta.iloc[ix] # needs to be between 0  and 1697
    _, _, dist_fold = cs.prepare_forcing_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder)
    era5_dist_files = sorted(glob.glob(str(dist_fold/'ERA5_dist_*.nc'))) # list
    [file_list.append(file) for file in era5_dist_files if file > last_complete]

# This gives us 919 files, meaning we need array 0-918 for indexing this list

# Select the file based on SLURM_ARRAY_ID that is the input to this script
ix = int(sys.argv[1])
file = file_list[ix]

print(f'Processing {file}')
with nc4.Dataset(file, 'r+') as f:
    # Add wind direction, function argument 'dims' toggles a switch away from
    #  default behavior in make_nc_variable() which assumes lat & lon dims exist
    if 'hru' in f.dimensions:
        f = cs.derive_wind_direction(f,dims='hru')
    else:
        f = cs.derive_wind_direction(f)
