# Average forcing grids into HRUs
import glob
import shutil
import sys
import pandas as pd
from pathlib import Path
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

# --- Reruns 2024-05-20
# These fix various small errors discovered during data use
rerun_file = Path('/globalhome/wmk934/HPC/camels_spat/7_forcing_data/forcing_check_logs/reruns_20240516.csv')
reruns = pd.read_csv(rerun_file)
# --- Reruns 2024-05-20

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
    print("Usage: python 3a_average_focring_to_hrus.py <array_index>")
    sys.exit(1)
else:
    ix = int(sys.argv[1])

# Define cs_meta row
row = cs_meta.iloc[ix] # needs to be between 0  and 1697

# --- Reruns 2024-05-19
this_basin = row['Country'] + '_' + row['Station_id']
if this_basin not in reruns['basin'].values:
    print(f'No reruns for basin {this_basin}. Exiting.')
    sys.exit(0) # with next station, because we have no reruns for this station. Error code 0: clean exit, no problems
else:
    print(f'Running reruns for basin {this_basin}.')
# --- Reruns 2024-05-19

# Set the spacing
era_spacing = 0.25
eme_spacing = 0.10

# Check if we need to run downloads for this station at all
missing = cs.flow_obs_unavailable(cs_unusable, row.Country, row.Station_id)
if 'iv' in missing and 'dv' in missing: 
    sys.exit(0) # exit gracefully, because we have no observations at all for this station

# Get shapefile path to determine download coordinates, and forcing destination path
basin_id, shp_lump_path, shp_dist_path, _, _ = cs.prepare_delineation_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder)
raw_fold, lump_fold, dist_fold = cs.prepare_forcing_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder) # Returns folders only, not file names
shp_dist_path = Path( str(shp_dist_path).format('basin') )
print('--- Now running basin {}. {}'.format(ix, basin_id))
print(' DEBUG WARNING: running EM-EARTH only')

# Ensure we have the CRS set in these shapes, because EASYMORE needs this
for shp in [shp_lump_path, shp_dist_path]:
    cs.add_crs_to_shapefile(shp)

# Get the forcing files
eme_merged_files = sorted(glob.glob(str(raw_fold/'EM_Earth_[0-9][0-9][0-9][0-9]-[0-9][0-9].nc'))) # list
era_merged_files = sorted(glob.glob(str(raw_fold/'ERA5_[0-9][0-9][0-9][0-9]-[0-9][0-9].nc'))) # list
era_invariant = glob.glob(str(raw_fold/'ERA5_*_invariants.nc'))

# Make forcing shapefiles
era_grid_shp, eme_grid_shp = cs.prepare_forcing_grid_shapefiles(row.Country, row.Station_id, Path(data_path)/cs_basin_folder)
for infile, outfile in zip([era_merged_files[0],eme_merged_files[0]], [era_grid_shp,eme_grid_shp]):
    cs.make_forcing_grid_shapefile(infile,outfile)

# Add geopotential to ERA5 forcing grid shapefile
cs.add_geopotential_to_era5_grid(era_invariant[0], era_grid_shp)

# Prepare for remapping
esmr_temp = cs.prepare_easymore_temp_folder(row.Country, row.Station_id, Path(data_path)/cs_basin_folder)

# Check if can do the remapping with EASYMORE with files as-is, and act accordingly
# reruns 2024-05-20: disabling ERA5 because the error we're fixing occurred with the EM-Earth forcing
era5_can_remap = cs.check_can_remap_as_is(era_merged_files[0]) # we can assume that if this applies to one file, it applies to all
#if era5_can_remap:
#    print('Remapping ERA5')
#    era_lump_esmr = cs.easymore_workflow('ERA5', 'lumped', esmr_temp, era_grid_shp, shp_lump_path, lump_fold, era_merged_files)
#    era_dist_esmr = cs.easymore_workflow('ERA5', 'dist',   esmr_temp, era_grid_shp, shp_dist_path, dist_fold, era_merged_files)
#else: 
#    # Files are 1x1 (lat x lon), use the workflow that adds padding of empty cells around this so we can keep using EASYMORE
#    era_lump_esmr = cs.easymore_workflow_with_cell_padding('ERA5', 'lumped', esmr_temp, era_grid_shp, shp_lump_path, lump_fold, 
#                                                        era_merged_files, grid_spacing=era_spacing)
#    era_dist_esmr = cs.easymore_workflow_with_cell_padding('ERA5', 'dist',   esmr_temp, era_grid_shp, shp_dist_path, dist_fold, 
#                                                        era_merged_files, grid_spacing=era_spacing)
#
# Temporary lines for reruns - we need these easymore objects for later plotting
era_lump_esmr, _ = cs.get_easymore_settings('ERA5', 'lumped', era_grid_shp, shp_lump_path, esmr_temp, lump_fold)
era_lump_esmr.source_nc = era_merged_files[-1] 
era_dist_esmr, _ = cs.get_easymore_settings('ERA5', 'dist', era_grid_shp, shp_dist_path, esmr_temp, dist_fold)
era_dist_esmr.source_nc = era_merged_files[-1] 
# reruns 2024-05-20

# Repeat for EM-Earth: because EM-Earth has a smaller spacing than ERA5, it is possible that we can remap one but not the other, hence separate
eme_can_remap = cs.check_can_remap_as_is(eme_merged_files[0]) # we can assume that if this applies to one file, it applies to all
if eme_can_remap:
    print('Remapping EM-Earth')
    eme_lump_esmr = cs.easymore_workflow('EM-Earth', 'lumped', esmr_temp, eme_grid_shp, shp_lump_path, lump_fold, eme_merged_files)
    eme_dist_esmr = cs.easymore_workflow('EM-Earth', 'dist',   esmr_temp, eme_grid_shp, shp_dist_path, dist_fold, eme_merged_files)
else:
    eme_lump_esmr = cs.easymore_workflow_with_cell_padding('EM-Earth', 'lumped', esmr_temp, eme_grid_shp, shp_lump_path, lump_fold, 
                                                        eme_merged_files, grid_spacing=eme_spacing)
    eme_dist_esmr = cs.easymore_workflow_with_cell_padding('EM-Earth', 'dist',   esmr_temp, eme_grid_shp, shp_dist_path, dist_fold, 
                                                        eme_merged_files, grid_spacing=eme_spacing)

# Create a graphical check of what we just did
fig_file = esmr_temp.parent / f'{row.Country}_{row.Station_id}_spatial_averaging.png'
cs.era5_eme_easymore_plotting_loop( [era_lump_esmr, era_dist_esmr, eme_lump_esmr, eme_dist_esmr], esmr_temp, fig_file )

fig_file = Path('/gpfs/tp/gwf/gwf_cmt/wknoben/camels_spat/camels-spat-data/TEMP_images_forcing_averaging') / f'{row.Country}_{row.Station_id}_spatial_averaging.png'
cs.era5_eme_easymore_plotting_loop( [era_lump_esmr, era_dist_esmr, eme_lump_esmr, eme_dist_esmr], esmr_temp, fig_file )

# Clean up
shutil.rmtree(esmr_temp)

# --- Warnings
# Restore the original warning filter settings
warnings.filters = original_filter








