# Fix RDRS PET units
# We did the PET computations in `kg m-2 s-1`, 
# but stored the `units` attribute as `mm hr-1`. 
# Here we loop over all the RDRS files and 
# correct this oversight.

import glob
import sys
import pandas as pd
from pathlib import Path
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

from itertools import chain
from netCDF4 import Dataset

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

# --- Data loading
# CAMELS-spat metadata file
cs_meta_path = Path(data_path) / cs_meta_path
cs_meta = pd.read_csv(cs_meta_path / cs_meta_name)

# Open list of unusable stations; Enforce reading IDs as string to keep leading 0's
cs_unusable = pd.read_csv(cs_meta_path / cs_unusable_name, dtype={'Station_id': object})

# --- Handle command line inputs
if len(sys.argv) != 2:
    print("Usage: python 7e_fix_rdrs_pet_units.py <array_index>")
    sys.exit(1)
else:
    ix = int(sys.argv[1])

# Define cs_meta row
row = cs_meta.iloc[ix] # needs to be between 0  and 1697
basin_id = row.Country + '_' + row.Station_id

# Check if we need to run downloads for this station at all
missing = cs.flow_obs_unavailable(cs_unusable, row.Country, row.Station_id)
if 'iv' in missing and 'dv' in missing: 
    sys.exit(0) # gracefully exit run, because we have no observations at all for this station

# --- Run the actual code
def netcdf_pet_units_mm_hr_to_kg_m2_s1(file, var_name='pet'):
    with Dataset(file, 'a') as nc_file:
        variable = nc_file.variables[var_name]
        if variable.units == 'mm hr**-1':
            variable.units = 'kg m**-2 s**-1'

# Get the paths
met_folder = basins_path / 'basin_data' / basin_id / 'forcing'

# Find the RDRS files
sub_folders = ['raw', 'lumped', 'distributed']
rdrs_files = []
for sub_folder in sub_folders:
    rdrs_files.append(glob.glob(str(met_folder / sub_folder / 'RDRS_*.nc')))

# Flatten the list
rdrs_files = list(chain.from_iterable(rdrs_files))
rdrs_files.sort()

# Fix the units
print(f"updated PET units in:")
for file in rdrs_files:
    netcdf_pet_units_mm_hr_to_kg_m2_s1(file)
    print(f' - {file}')