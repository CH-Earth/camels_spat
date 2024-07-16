# Average Daymet forcing data into HRUs
import shutil
import sys
import pandas as pd
from pathlib import Path
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

# Figure out which easymore version we're dealing with and load accordingly
import pkg_resources

def get_package_version(package_name):
    try:
        version = pkg_resources.get_distribution(package_name).version
        return version
    except pkg_resources.DistributionNotFound:
        return f"Package '{package_name}' not found."

package_name = 'easymore'
esmr_version = get_package_version(package_name)
print(f"The version of {package_name} is {esmr_version}")

#if esmr_version == '2.0.0':
#    from easymore import Easymore
#else:
#    from easymore import easymore

# --- Warnings
import warnings

# Save the current warning filter settings
original_filter = warnings.filters[:]

# Temporarily filter out future warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

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
# Get the array index from the command-line argument
if len(sys.argv) != 2:
    print("Usage: python 6g_daymet_to_subbasins.py <array_index>")
    sys.exit(1)
else:
    ix = int(sys.argv[1])

# Temporary storage for easymore files
tmp_dir_main = Path('/scratch/gwf/gwf_cmt/wknoben/easymore_tmp')
 
# Define cs_meta row
row = cs_meta.iloc[ix] # needs to be between 0  and 1697

# Check if we need to run downloads for this station at all
missing = cs.flow_obs_unavailable(cs_unusable, row.Country, row.Station_id)
if 'iv' in missing and 'dv' in missing: 
    sys.exit(0) # exit gracefully, because we have no observations at all for this station

# Get shapefile path to determine download coordinates, and forcing destination path
basin_id, shp_lump_path, shp_dist_path, _, _ = cs.prepare_delineation_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder)
raw_fold, lump_fold, dist_fold = cs.prepare_forcing_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder) # Returns folders only, not file names
shp_dist_path = Path( str(shp_dist_path).format('basin') )
print('--- Now running basin {}. {}'.format(ix, basin_id))
print(' DEBUG WARNING: Dev run')

# Determine where the forcing grid needs to go
forcing_grid_folder = Path(data_path) / 'camels-spat-data' / 'basin_data' / basin_id / 'shapefiles' / 'forcing_grids'

# Define the temporary folder location
tmp_dir = tmp_dir_main / basin_id
tmp_dir.mkdir(parents=True, exist_ok=True)

# Do the remapping for the lumped and distributed cases
for case, basin_shp, out_dir in zip(['lumped', 'dist'],
                                    [shp_lump_path, shp_dist_path], 
                                    [lump_fold, dist_fold]):

    # Define the EASYMORE inputs
    esmr,_ = cs.get_easymore_settings('Daymet', case, None, basin_shp, tmp_dir, out_dir,
                                      in_files=str(raw_fold / 'daymet_*.nc'), version=esmr_version)
    
    # Set the approximate grid resolution for cases where we only have a single grid cell
    # We'll overestimate this a little, so we don't accidentally cut part of the basin of
    # degrees 1km / 111km = 0.0090 degress latitude, longitude would be smaller but we can 
    # only specify one value in EASYMORE v2.
    esmr.source_nc_resolution = 0.1 

    esmr.nc_remapper()

    # If we're in the lumped case, retain EASYMORE's forcing grid 
    if case == 'lumped':
        tmp_grid_exts = ['cpg', 'dbf', 'prj', 'shp', 'shx']
        tmp_grid_file = esmr.case_name + '_source_shapefile_corrected_frame.'
        des_grid_file = f'Daymet_grid_{basin_id}.'
        for ext in tmp_grid_exts:
            src = tmp_dir / (tmp_grid_file + ext)
            des = forcing_grid_folder / (des_grid_file + ext)
            shutil.copy(src, des)

# Remove the temporary folder (disabled for now, until we know everything works)
#shutil.rmtree(tmp_dir)

# --- Warnings
# Restore the original warning filter settings
warnings.filters = original_filter


# # --- DEV --- Define the EASYMORE inputs
# esmr = easymore() # initialize object
# 
# esmr.author_name = 'CAMELS-spat workflow'
# esmr.case_name = 'daymet_lumped'
# esmr.save_csv  = False # Flag that we do not want the data stored in .csv in addition to .nc
# esmr.remap_csv = '' # Initialize with 'no remap file available' - we'll update this later
# esmr.sort_ID = False # Enforce that we want our HRUs returned in the order we put them in
# 
# esmr.temp_dir = str(tmp_dir) + '/' # easymore demands this extra /
# 
# esmr.target_shp = str(shp_lump_path)
# esmr.target_shp_ID = 'FID'
# 
# esmr.source_nc = str(raw_fold / 'daymet_*.nc')
# esmr.var_names = ['prcp', 'tmax', 'tmin', 'srad', 'vp', 'dayl', 'pet']
# esmr.var_lon = 'lon'
# esmr.var_lat = 'lat'
# esmr.var_time = 'time'
# 
# esmr.output_dir = str(lump_fold)  + '/' # easymore demands this extra /
# 
# esmr.format_list = ['f4'] # everything as float32, like it is in the source file
# esmr.fill_value_list = [-9999.0] # both expand automatically to the length of var_names
# 
# remap_csv = esmr.case_name + '_remapping.csv'
# 
# esmr.nc_remapper()
# 
# # Copy EASYMORE's forcing grid shapefile to the basin's shapefile folder
# tmp_grid_exts = ['cpg', 'dbf', 'prj', 'shp', 'shx']
# tmp_grid_file = esmr.case_name + '_source_shapefile_corrected_frame.'
# des_grid_file = f'Daymet_grid_{basin_id}.'
# for ext in tmp_grid_exts:
#     src = tmp_dir / (tmp_grid_file + ext)
#     des = forcing_grid_folder / (des_grid_file + ext)
#     shutil.copy(src, des)
# # --- END DEV ---