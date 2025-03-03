# Restructure data
# This takes the work-in-progress files and moves them into the final data structure we want to upload.

# Almost all of the data is catchment-based so we can put this into a parallel run. Things to remember:
# - Do not redistribute the raw WorldClim data - this is not allowed.
# - Remember to put the main attribute file into the resulting attributes folder.

# imports
import geopandas as gpd
import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import shutil
import sys
import warnings
import xarray as xr

from datetime import datetime
from pathlib import Path

sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs


# - CONFIG HANDLING
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

# Attributes folder
cs_att_folder = cs.read_from_config(config_file, 'att_path')
att_path  = basins_path / 'camels_spat_attributes.csv'

# Destination folder
final_fold = cs.read_from_config(config_file, 'final_path')


# -- DATA LOADING
# CAMELS-spat metadata file
cs_meta_path = Path(data_path) / cs_meta_path
cs_meta = pd.read_csv(cs_meta_path / cs_meta_name)

# Open list of unusable stations; Enforce reading IDs as string to keep leading 0's
cs_unusable = pd.read_csv(cs_meta_path / cs_unusable_name, dtype={'Station_id': object})


# -- ARGUMENTS
if len(sys.argv) != 2:
    print("Usage: python 8b_restructure_basin.py <array_index>")
    sys.exit(1)
else:
    ix = int(sys.argv[1])


# -- GENERAL
# Set the top level path
dest_root = Path(final_fold) / 'camels-spat-upload'

# Subset Meta-data to the 1426 basins we want
cs_meta_upload = cs_meta[~cs_meta.set_index(['Country', 'Station_id']).index.isin(cs_unusable.set_index(['Country', 'Station_id']).index)]

# Ensure we have the expected number of basins (1426)
meta_basins = (cs_meta_upload['Country'] + '_' + cs_meta_upload['Station_id']).values
assert len(meta_basins) == 1426, "Number of basins not 1426"


# -- PROCESSING
# Get the basin's info
row = cs_meta_upload.iloc[ix]

# Set the spatial category
category = row['subset_category']

# Get the location of the source files
basin_id   = row['Country'] + '_' + row['Station_id']
geo_folder = basins_path / 'basin_data' / basin_id / 'geospatial'
met_folder = basins_path / 'basin_data' / basin_id / 'forcing'
hyd_folder = basins_path / 'basin_data' / basin_id / 'observations'
shp_folder = basins_path / 'basin_data' / basin_id / 'shapefiles'
att_folder = basins_path / 'attributes'

# Track the main files we expect
copied_files = []
expected_files = [
    # shapefiles excluding reference files because we know we don't have them for all cases
    'shapefiles-delineation-plot',
    'shapefiles-dist-basin',
    'shapefiles-dist-river',
    'shapefiles-forc-daymet',
    'shapefiles-forc-emearth',
    'shapefiles-forc-era5',
    'shapefiles-forc-rdrs',
    'shapefiles-lump-basin',
    # observations
    'observations-daily', # we know these exist because if not the shutil.copy() will fail
    'observations-hourly',# must exist or copy error
    # forcing
    'forcing-daym-grid',
    'forcing-daym-lump', 
    'forcing-daym-dist',
    'forcing-emea-grid',
    'forcing-emea-lump',
    'forcing-emea-dist',
    'forcing-era5-grid',
    'forcing-era5-lump',
    'forcing-era5-dist',
    'forcing-era5-inva',
    'forcing-rdrs-grid',
    'forcing-rdrs-lump',
    'forcing-rdrs-dist',
    # attributes
    'attributes-dist', # must exist or copy error
    # geospatial
    'geospatial-forest-height-2000', # must exist or copy error
    'geospatial-forest-height-2020', # must exist or copy error
    'geospatial-glclu2019-map', # must exist or copy error
    'geospatial-glclu2019-strata', # must exist or copy error
    'geospatial-glhymps',
    'geospatial-hydrolakes',
    'geospatial-lai',
    'geospatial-lgrip30', # must exist or copy error
    'geospatial-merit-aspect', # must exist or copy error
    'geospatial-merit-dem', # must exist or copy error
    'geospatial-merit-slope', # must exist or copy error
    'geospatial-modis-land-mode', # must exist or copy error
]


# - SHAPEFILES

# 1. Shapefiles
# Source:      /basin_data/[basin_id]/shapesfiles/distributed
#                                                /forcing_grids
#                                                /lumped
#                                                /reference
#                                                /[basin_id]_delineation_results.png
#
# Destination: /camels-spat-upload/shapefiles/[category]/delineation-outcomes
#                                                       /shapes-distributed
#                                                       /shapes-forcing
#                                                       /shapes-lumped
#                                                       /shapes-reference
#
# Too keep things clean we'll need to create a dedicated basin folder inside each 
# of the 'shapes' subfolders in destination.

# 1a. Delineation plot
# This will error automatically if the source file does not exist as specified
file = f"{basin_id}_delineation_results.png"
src = shp_folder / file
dst = dest_root / 'shapefiles' / category / 'delineation-outcomes' / file
shutil.copy(src, dst)
copied_files.append('shapefiles-delineation-plot')

# 1b. Distributed shapefiles
dst = dest_root / 'shapefiles' / category / 'shapes-distributed' / basin_id
dst.mkdir(exist_ok=True)

src_files = glob.glob(str(shp_folder / 'distributed' / f"{basin_id}_distributed_*"))
for src_file in src_files:
    file = os.path.basename(src_file)
    dst_file = dst / file
    shutil.copy(src_file, dst_file)
    if file == f"{basin_id}_distributed_basin.shp": copied_files.append('shapefiles-dist-basin')
    if file == f"{basin_id}_distributed_river.shp": copied_files.append('shapefiles-dist-river')

# Run a quick fix on the distributed river shapefile, to drop the existing 'lengthkm' column
river_file = dst / f"{basin_id}_distributed_river.shp"
gdf = gpd.read_file(river_file)
gdf = gdf.drop(columns=['lengthkm'])
if len(gdf) == 0: # empty shapefile, and this will interfere with saving
    empty_row = pd.DataFrame([{col: np.nan if col != 'geometry' else None for col in gdf.columns}]) # add Nan/None (this may also help users)
    gdf = pd.concat([gdf, empty_row], ignore_index=True)
gdf.to_file(river_file)

# 1c. Forcing shapefiles
dst = dest_root / 'shapefiles' / category / 'shapes-forcing' / basin_id
dst.mkdir(exist_ok=True)

src_files = glob.glob(str(shp_folder / 'forcing_grids' / f"*_grid_{basin_id}.*"))
for src_file in src_files:
    file = os.path.basename(src_file)
    file = file.replace(f"_{basin_id}","")     # swap the basin id to the front
    dst_file = f"{basin_id}_{file}"
    dst_file = dst / dst_file
    shutil.copy(src_file, dst_file)
    if os.path.basename(dst_file) == f"{basin_id}_Daymet_grid.shp": copied_files.append('shapefiles-forc-daymet')
    if os.path.basename(dst_file) == f"{basin_id}_EM_Earth_grid.shp": copied_files.append('shapefiles-forc-emearth')
    if os.path.basename(dst_file) == f"{basin_id}_ERA5_grid.shp": copied_files.append('shapefiles-forc-era5')
    if os.path.basename(dst_file) == f"{basin_id}_RDRS_grid.shp": copied_files.append('shapefiles-forc-rdrs')

# 1d. Lumped shapefiles
dst = dest_root / 'shapefiles' / category / 'shapes-lumped' / basin_id
dst.mkdir(exist_ok=True)
src_files = glob.glob(str(shp_folder / 'lumped' / f"{basin_id}_lumped.*"))
for src_file in src_files:
    file = os.path.basename(src_file)
    dst_file = dst / file
    shutil.copy(src_file, dst_file)
    if file == f"{basin_id}_lumped.shp": copied_files.append('shapefiles-lump-basin')

# 1e. Reference shapefiles
dst = dest_root / 'shapefiles' / category / 'shapes-reference' / basin_id
dst.mkdir(exist_ok=True)
src_files = glob.glob(str(shp_folder / 'reference' / f"{basin_id}_reference.*"))
if len(src_files) > 0:
    for src_file in src_files:
        file = os.path.basename(src_file)
        dst_file = dst / file
        shutil.copy(src_file, dst_file)
else:
    with open(str(dst/f"{basin_id}.txt"), 'w') as f:
        f.write(f"No reference file available for basin {basin_id}")


# -- OBSERVATIONS

# 2. Observations
# Source:      /basin_data/[basin_id]/observations/[basin_id]_daily_flow_observations.nc
#                                                 /[basin_id]_hourly_flow_observations.nc
#
# Destination: /camels-spat-upload/observations/[category]/obs-daily/[basin_id]_daily_flow_observations.nc
#                                                         /obs-hourly/[basin_id]_hourly_flow_observations.nc

# 2a. Daily
file = f"{basin_id}_daily_flow_observations.nc"
src = hyd_folder / file
dst = dest_root / 'observations' / category / 'obs-daily' / file
shutil.copy(src, dst)
copied_files.append('observations-daily')

# 2b. Hourly
file = f"{basin_id}_hourly_flow_observations.nc"
src = hyd_folder / file
dst = dest_root / 'observations' / category / 'obs-hourly' / file
shutil.copy(src, dst)
copied_files.append('observations-hourly')


# -- FORCING

def check_netcdf_compatibility(file1, file2, file3, expected_vars, expected_time_steps):
    """
    Opens three NetCDF files and checks if they:
    - Cover the same time period.
    - Contain the expected variables (also ensures variables have a 'time' dimension).
    - Have the expected number of time steps.

    Parameters:
    - file1, file2, file3 (str): File paths to the NetCDF files.
    - expected_vars (list): List of expected variables that should have a "time" dimension.
    - expected_time_steps (int): Expected number of time steps in each file.

    Returns:
    - dict: A dictionary with results for time period matching, variable existence, and time step count.
    """

    # Open the NetCDF files
    ds1 = xr.open_dataset(file1)
    ds2 = xr.open_dataset(file2)
    ds3 = xr.open_dataset(file3)

    # Extract time ranges (if time exists in dataset)
    time_range1 = (ds1["time"].min().values, ds1["time"].max().values) if "time" in ds1 else None
    time_range2 = (ds2["time"].min().values, ds2["time"].max().values) if "time" in ds2 else None
    time_range3 = (ds3["time"].min().values, ds3["time"].max().values) if "time" in ds3 else None

    # Check if time periods match
    time_match = (time_range1 == time_range2 == time_range3)

    # Extract the number of time steps in each dataset
    time_steps1 = len(ds1["time"]) if "time" in ds1 else None
    time_steps2 = len(ds2["time"]) if "time" in ds2 else None
    time_steps3 = len(ds3["time"]) if "time" in ds3 else None

    # Check if time step counts match the expected value
    time_step_match = (time_steps1 == time_steps2 == time_steps3 == expected_time_steps)

    # Check if expected variables exist in each dataset
    var_existence = {
        var: {
            "file1": var in ds1.data_vars and "time" in ds1[var].dims,
            "file2": var in ds2.data_vars and "time" in ds2[var].dims,
            "file3": var in ds3.data_vars and "time" in ds3[var].dims,
        }
        for var in expected_vars
    }

    # Find any missing variables
    missing_vars = {
        var: [file for file in ["file1", "file2", "file3"] if not var_existence[var][file]]
        for var in expected_vars
        if not all(var_existence[var].values())  # Only keep variables that are missing in at least one file
    }

    # Find any files with incorrect time step count
    incorrect_time_steps = {
        file: actual_steps
        for file, actual_steps in {
            "file1": time_steps1,
            "file2": time_steps2,
            "file3": time_steps3,
        }.items()
        if actual_steps != expected_time_steps
    }

    # Raise assertion error if any expected variable is missing in any file
    assert not missing_vars, f"Missing expected variables in files: {missing_vars}"

    # Raise assertion error if any file has an incorrect number of time steps
    assert not incorrect_time_steps, f"Incorrect time steps in files: {incorrect_time_steps}. Expected {expected_time_steps}."

    # Close datasets to free memory
    ds1.close()
    ds2.close()
    ds3.close()

    return

def generate_encoding(nc, dim_order, time_chunk=100):
    """
    Generates an encoding dictionary for NetCDF compression and chunking.

    Parameters:
    - nc (xr.Dataset): Xarray dataset containing the variables and dimensions.
    - dim_order (tuple): Tuple specifying the expected dimension order (e.g., ("time", "rlat", "rlon")).
    - time_chunk (int): Chunk size for the "time" dimension.

    Returns:
    - dict: Encoding dictionary for use with `to_netcdf()`, ensuring 'source' and 'coordinates' are removed.
    """

    # Create chunk_sizes where "time" has "time_chunk" size, and other dimensions keep their full length
    chunk_sizes = tuple(time_chunk if dim == "time" else nc.dims[dim] for dim in dim_order)

    # Find variables that match the specified dimension order
    matching_vars = [var for var in nc.data_vars if nc[var].dims == dim_order]

    # Construct encoding dictionary, preserving existing settings
    encoding = {}
    for var in matching_vars:
        # Copy existing encoding or start with an empty dict
        existing_encoding = nc[var].encoding.copy()

        # Update only the necessary keys
        updated_encoding = {
            **existing_encoding,  # Retain existing encoding
            "zlib": True,
            "complevel": 4,
            "chunksizes": chunk_sizes,
        }

        # Remove unwanted encoding keys
        for key in ["source", "coordinates"]:
            updated_encoding.pop(key, None)  # Remove if exists

        # Assign updated encoding back to dictionary
        encoding[var] = updated_encoding

    return encoding

def fix_nc_encoding(ds, variables=None, missing_values=-999.0):
    """
    Replaces specified missing_values in an xarray Dataset with NaN and updates encoding.

    Parameters:
    - ds (xr.Dataset): The xarray Dataset to modify.
    - variables (list, optional): List of variable names to process. If None, all variables are processed.
    - missing_values (float, int, or list): Single missing value or a list of values matching the variables.

    Returns:
    - xr.Dataset: The updated Dataset with NaNs and modified history.
    """

    # If no variable list is provided, process all data variables
    if variables is None:
        variables = list(ds.data_vars)

    # Ensure missing_values is a list matching the variables list
    if not isinstance(missing_values, list):
        missing_values = [missing_values] * len(variables)  # Expand single value to list

    if len(missing_values) != len(variables):
        raise ValueError("Length of missing_values must match the number of variables.")

    modified_vars = []  # Track which variables were updated

    # Loop over specified variables and their corresponding missing values
    for var, missing_value in zip(variables, missing_values):
        if var in ds:
            # Preserve existing encoding while removing/updating _FillValue and missing_value
            updated_encoding = ds[var].encoding.copy()  # Copy the existing encoding
                        
            # Replace missing_value with NaN
            ds[var] = ds[var].where(ds[var] != missing_value, np.nan)
            
            # Update encoding
            updated_encoding["_FillValue"] = np.nan
            updated_encoding["missing_value"] = np.nan

            # Apply the updated encoding back
            ds[var].encoding = updated_encoding

            # Track changes
            modified_vars.append(f"{var} (Replaced {missing_value} with NaN)")

    # Update global history attribute (handling case insensitivity)
    if modified_vars:
        timestamp = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC")
        change_log = f"{timestamp}: Replaced FillValues with NaN in: {', '.join(modified_vars)}."

        # Find existing history key (case-insensitive search)
        history_key = next((key for key in ds.attrs if key.lower() == "history"), "History")

        # Append to existing history
        existing_history = ds.attrs.get(history_key, "")
        ds.attrs[history_key] = f"{existing_history} {change_log}".strip()

    return ds

# 3. Forcing
# Source:      /basin_data/[basin_id]/forcing/distributed/[lots of individual files]
#                                            /lumped/[lots of individual files]
#                                            /raw/[lots of individual files]
#
# Destination: /camels-spat-upload/forcing/[category]/daymet/daymet-distributed
#                                                    /daymet/daymet-gridded
#                                                    /daymet/daymet-lumped
#                                                    /em-earth
#                                                    /era5
#                                                    /rdrs
#
# We need to merge the individual files for each basin into a single one

# 3a. DAYMET - 
# Special treatment: we'll load these without decoding times (decode_times=False) to avoid:
# UserWarning: Variable 'time' has datetime type and a bounds variable but time.encoding 
# does not have units specified. The units encodings for 'time' and 'time_bnds' will be 
# determined independently and may not be equal, counter to CF-conventions. If this is a 
# concern, specify a units encoding for 'time' before writing to a file.
with warnings.catch_warnings():
    warnings.simplefilter("ignore", UserWarning)

    # gridded
    dst1 = dest_root / 'forcing' / category / 'daymet' / 'daymet-gridded' / f"{basin_id}_daymet_gridded.nc"
    src_files = glob.glob(str(met_folder / 'raw' / f"*aymet_*.nc"))
    nc = xr.open_mfdataset(src_files, combine="by_coords")#, decode_times=False)
    expected_length = len(nc['time'])
    encoding = generate_encoding(nc, ("time","x","y"))
    nc.to_netcdf(dst1, encoding=encoding)
    nc.close()
    
    # lumped
    dst2 = dest_root / 'forcing' / category / 'daymet' / 'daymet-lumped' / f"{basin_id}_daymet_lumped.nc"
    src_files = glob.glob(str(met_folder / 'lumped' / f"*aymet_*.nc"))
    nc = xr.open_mfdataset(src_files, combine="by_coords")#, decode_times=False)
    encoding = generate_encoding(nc, ("time","hru"))
    nc.to_netcdf(dst2, encoding=encoding)
    nc.close()
    
    # distributed
    dst3 = dest_root / 'forcing' / category / 'daymet' / 'daymet-distributed' / f"{basin_id}_daymet_distributed.nc"
    src_files = glob.glob(str(met_folder / 'distributed' / f"*aymet_*.nc"))
    nc = xr.open_mfdataset(src_files, combine="by_coords")#, decode_times=False)
    encoding = generate_encoding(nc, ("time","hru"))
    nc.to_netcdf(dst3, encoding=encoding)
    nc.close()

    # Check
    expected_vars = ['dayl','time_bnds','pet','prcp','srad','tmax','tmin','vp']
    check_netcdf_compatibility(dst1, dst2, dst3, expected_vars, expected_length)

    copied_files.append('forcing-daym-grid')
    copied_files.append('forcing-daym-lump')
    copied_files.append('forcing-daym-dist')


# 3b. EM-Earth - 
# Special treatment: we can't use decode_times=False, and combine="by_coords" at the same time
# because every monthly EM-Earth file has a different time encoding: they are all "hours since
# start of month" and thus, when decoded, are all 0..744 (or something, depending on month).
# When combining, this simply overwrites everything. It is probably better to use 
# combine="by_coords" to ensure we retain the correct temporal order.
with warnings.catch_warnings():
    warnings.simplefilter("ignore", UserWarning)

    # gridded
    dst1 = dest_root / 'forcing' / category / 'em-earth' / 'em-earth-gridded' / f"{basin_id}_em_earth_gridded.nc"
    src_files = glob.glob(str(met_folder / 'raw' / f"EM*.nc"))
    nc = xr.open_mfdataset(src_files, combine="by_coords")
    expected_length = len(nc['time']) # used later to check
    nc = fix_nc_encoding(nc, variables=['tmean','prcp'],missing_values=[
            nc['tmean'].encoding['missing_value'],nc['prcp'].encoding['missing_value']])
    nc = fix_nc_encoding(nc, variables=['tmean','prcp'],missing_values=-9999.0) # just to be sure
    encoding = generate_encoding(nc, ("time","latitude","longitude"))
    nc.to_netcdf(dst1, encoding=encoding)
    nc.close()
    
    # lumped
    dst2 = dest_root / 'forcing' / category / 'em-earth' / 'em-earth-lumped' / f"{basin_id}_em_earth_lumped.nc"
    src_files = glob.glob(str(met_folder / 'lumped' / f"EM*.nc"))
    nc = xr.open_mfdataset(src_files, combine="by_coords")
    encoding = generate_encoding(nc, ("time","hru"))
    nc.to_netcdf(dst2, encoding=encoding)
    nc.close()
    
    # distributed
    dst3 = dest_root / 'forcing' / category / 'em-earth' / 'em-earth-distributed' / f"{basin_id}_em_earth_distributed.nc"
    src_files = glob.glob(str(met_folder / 'distributed' / f"EM*.nc"))
    nc = xr.open_mfdataset(src_files, combine="by_coords")
    encoding = generate_encoding(nc, ("time","hru"))
    nc.to_netcdf(dst3, encoding=encoding)
    nc.close()

    # Check
    expected_vars = ['tmean','time_bnds','prcp']
    check_netcdf_compatibility(dst1, dst2, dst3, expected_vars, expected_length)

    copied_files.append('forcing-emea-grid')
    copied_files.append('forcing-emea-lump')
    copied_files.append('forcing-emea-dist')


# 3c. ERA5
with warnings.catch_warnings():
    warnings.simplefilter("ignore", UserWarning)

    # gridded
    dst1 = dest_root / 'forcing' / category / 'era5' / 'era5-gridded' / f"{basin_id}_era5_gridded.nc"
    src_files = glob.glob(str(met_folder / 'raw' / f"ERA5*.nc"))
    src_files = [file for file in src_files if not 'invariant' in file]
    nc = xr.open_mfdataset(src_files, combine="by_coords")
    expected_length = len(nc['time']) # used later to check
    nc = fix_nc_encoding(nc, variables=[
        'mper','msdwlwrf','msdwswrf','msnlwrf','msnswrf','mtpr','q','sp','t','u','v'], 
                    missing_values=-999.0)
    encoding = generate_encoding(nc, ("time","latitude","longitude"))
    nc.to_netcdf(dst1, encoding=encoding)
    nc.close()
    
    # lumped
    dst2 = dest_root / 'forcing' / category / 'era5' / 'era5-lumped' / f"{basin_id}_era5_lumped.nc"
    src_files = glob.glob(str(met_folder / 'lumped' / f"ERA5*.nc"))
    nc = xr.open_mfdataset(src_files, combine="by_coords")
    encoding = generate_encoding(nc, ("time","hru"))
    nc.to_netcdf(dst2, encoding=encoding)
    nc.close()
    
    # distributed
    dst3 = dest_root / 'forcing' / category / 'era5' / 'era5-distributed' / f"{basin_id}_era5_distributed.nc"
    src_files = glob.glob(str(met_folder / 'distributed' / f"ERA5*.nc"))
    nc = xr.open_mfdataset(src_files, combine="by_coords")
    encoding = generate_encoding(nc, ("time","hru"))
    nc.to_netcdf(dst3, encoding=encoding)
    nc.close()
    
    # invariant file
    dst_file = dest_root / 'forcing' / category / 'era5' / 'era5-invariants' / f"{basin_id}_era5_invariants.nc"
    src_file = glob.glob(str(met_folder / 'raw' / f"ERA5*_invariants.nc"))
    assert len(src_file) == 1
    shutil.copy(src_file[0], dst_file)

    # Check
    expected_vars = ['e','mper','msdwlwrf','msdwswrf','msnlwrf','msnswrf','mtpr','q','rh',
                     'sp','t','u','v','w','phi','time_bnds']
    check_netcdf_compatibility(dst1, dst2, dst3, expected_vars, expected_length)
    
    copied_files.append('forcing-era5-grid')
    copied_files.append('forcing-era5-lump')
    copied_files.append('forcing-era5-dist')
    copied_files.append('forcing-era5-inva')


# 3d. RDRS
with warnings.catch_warnings():
    warnings.simplefilter("ignore", UserWarning)

    # gridded
    dst1 = dest_root / 'forcing' / category / 'rdrs' / 'rdrs-gridded' / f"{basin_id}_rdrs_gridded.nc"
    src_files = glob.glob(str(met_folder / 'raw' / f"RDRS*.nc"))
    nc = xr.open_mfdataset(src_files, combine="by_coords")
    expected_length = len(nc['time']) # used later to check
    encoding = generate_encoding(nc,("time", "rlat", "rlon"))
    nc.to_netcdf(dst1, encoding=encoding)
    nc.close()
    
    # lumped
    dst2 = dest_root / 'forcing' / category / 'rdrs' / 'rdrs-lumped' / f"{basin_id}_rdrs_lumped.nc"
    src_files = glob.glob(str(met_folder / 'lumped' / f"RDRS*.nc"))
    nc = xr.open_mfdataset(src_files, combine="by_coords")
    encoding = generate_encoding(nc,("time", "hru"))
    nc.to_netcdf(dst2, encoding=encoding)
    nc.close()
    
    # distributed
    dst3 = dest_root / 'forcing' / category / 'rdrs' / 'rdrs-distributed' / f"{basin_id}_rdrs_distributed.nc"
    src_files = glob.glob(str(met_folder / 'distributed' / f"RDRS*.nc"))
    nc = xr.open_mfdataset(src_files, combine="by_coords")
    encoding = generate_encoding(nc,("time", "hru"))
    nc.to_netcdf(dst3, encoding=encoding)
    nc.close()

    # Check
    expected_vars = ['RDRS_v2.1_P_P0_SFC','RDRS_v2.1_P_TT_1.5m','RDRS_v2.1_P_HU_1.5m',
                     'RDRS_v2.1_P_HR_1.5m','RDRS_v2.1_P_UUC_10m','RDRS_v2.1_P_VVC_10m',
                     'RDRS_v2.1_P_UVC_10m','RDRS_v2.1_P_FI_SFC','RDRS_v2.1_P_FB_SFC',
                     'RDRS_v2.1_P_GZ_SFC','RDRS_v2.1_A_PR0_SFC','e','pet','phi','time_bnds']
    check_netcdf_compatibility(dst1, dst2, dst3, expected_vars, expected_length)
    
    copied_files.append('forcing-rdrs-grid')
    copied_files.append('forcing-rdrs-lump')
    copied_files.append('forcing-rdrs-dist')


# -- ATTRIBUTES
# 4. Attributes
# Source:      /attributes/distributed/attributes_[basin_id].csv
#
# Destination: /camels-spat-upload/attributes/[category]/[basin_id]_attributes.csv

file_src = f"attributes_{basin_id}.csv"
file_dst = f"{basin_id}_attributes.csv"
src = att_folder / 'distributed' / file_src
dst = dest_root / 'attributes' / category / file_dst
shutil.copy(src, dst)
copied_files.append('attributes-dist')


# -- GEOSPATIAL

# 5. Geospatial
# Source:       /basin_data/[basin_id]/geospatial/forest_height/raw
#                                                /glclu2019/raw
#                                                /glhymps/raw
#                                                /hydrolakes/raw
#                                                /lai/raw
#                                                /lgrip30/raw
#                                                /merit/aspect
#                                                /merit/raw
#                                                /merit/slope
#                                                /modis_land/raw
#                                                /pelletier/raw
#                                                /soilgrids/raw/bdod
#                                                /soilgrids/raw/cfvo
#                                                /soilgrids/raw/clay
#                                                /soilgrids/raw/conductivity
#                                                /soilgrids/raw/porosity
#                                                /soilgrids/raw/sand
#                                                /soilgrids/raw/silt
#                                                /soilgrids/raw/soc
#                                                /worldclim/raw/annual
#                                                /worldclim/raw/aridity2
#                                                /worldclim/raw/fracsnow2
#                                                /worldclim/raw/pet
#                                                /worldclim/raw/snow2
#
# Destination: /camels-spat-upload/attributes/[category]/forest_height/[basin_id]

# Special treatment: we'll need to make sure we put the basin_id in every single filename.
# This will avoid confusion for users.


# Forest height

# General
dst_fold = dest_root / 'geospatial' / category / 'forest-height' / basin_id
dst_fold.mkdir(exist_ok=True)

# 2000
file_src = "forest_height_2000.tif"
file_dst = f"{basin_id}_forest_height_2000.tif"
src = geo_folder / 'forest_height' / 'raw' / file_src
dst = dst_fold / file_dst
shutil.copy(src, dst)
copied_files.append('geospatial-forest-height-2000')

# 2020
file_src = "forest_height_2020.tif"
file_dst = f"{basin_id}_forest_height_2020.tif"
src = geo_folder / 'forest_height' / 'raw' / file_src
dst = dst_fold / file_dst
shutil.copy(src, dst)
copied_files.append('geospatial-forest-height-2020')


# GLCLU2019

# General
dst_fold = dest_root / 'geospatial' / category / 'glclu2019' / basin_id
dst_fold.mkdir(exist_ok=True)

# Map
file_src = "glclu2019_map.tif" # NOTE that we switch strata and map on purpose here, so that 
file_dst = f"{basin_id}_glclu2019_strata.tif" # the names match the legend in the excel file
src = geo_folder / 'glclu2019' / 'raw' / file_src
dst = dst_fold / file_dst
shutil.copy(src, dst)
copied_files.append('geospatial-glclu2019-strata')

# Strata
file_src = "glclu2019_strata.tif"
file_dst = f"{basin_id}_glclu2019_map.tif"
src = geo_folder / 'glclu2019' / 'raw' / file_src
dst = dst_fold / file_dst
shutil.copy(src, dst)
copied_files.append('geospatial-glclu2019-map')


# GLHYMPS

# General
dst_fold = dest_root / 'geospatial' / category / 'glhymps' / basin_id
dst_fold.mkdir(exist_ok=True)

# Find files
src_files = glob.glob(str(geo_folder / 'glhymps' / 'raw' / "glhymps.*"))
for src_file in src_files:
    file = os.path.basename(src_file)
    dst_file = dst_fold / f"{basin_id}_{file}"
    shutil.copy(src_file, dst_file)
    if os.path.basename(dst_file) == f"{basin_id}_glhymps.shp": copied_files.append('geospatial-glhymps')


# HYDROLAKES

# General
dst_fold = dest_root / 'geospatial' / category / 'hydrolakes' / basin_id
dst_fold.mkdir(exist_ok=True)

# Find files
src_files = glob.glob(str(geo_folder / 'hydrolakes' / 'raw' / "HydroLAKES*"))

# Check if we have a .shp file: this suggests we have an actual lake shapefile
if any('.shp' in src_file for src_file in src_files):
    for src_file in src_files:
        file = os.path.basename(src_file)
        dst_file = dst_fold / f"{basin_id}_{file}"
        shutil.copy(src_file, dst_file)
        if os.path.basename(dst_file) == f"{basin_id}_HydroLAKES_polys_v10_NorthAmerica.shp": copied_files.append('geospatial-hydrolakes')
# If we don't have a .shp file, confirm we have a .txt file: 
# This would ahve been generated if the hydrolakes data has no lake in this basin polygon
# In this case, add an empty shapefile. This is likely easier for users (all basins have same file
# just some are empty) than mixing file types (some basins have .shp, others have .txt)
elif any('.txt' in src_file for src_file in src_files):
    gdf = gpd.GeoDataFrame(
        [{'no_lakes': np.nan, 'geometry': None}],
        geometry='geometry',
        crs="EPSG:4326"
    )
    dst_file = dst_fold / f"{basin_id}_HydroLAKES_polys_v10_NorthAmerica.shp"
    gdf.to_file(dst_file)
    copied_files.append('geospatial-hydrolakes')


# LAI

# General
dst_fold = dest_root / 'geospatial' / category / 'lai' / basin_id
dst_fold.mkdir(exist_ok=True)

# Find files
src_files = glob.glob(str(geo_folder / 'lai' / 'raw' / "*_Lai_500m.tif"))
for src_file in src_files:
    file = os.path.basename(src_file)
    dst_file = dst_fold / f"{basin_id}_{file}"
    shutil.copy(src_file, dst_file)
    
if len(src_files) > 0:
    copied_files.append('geospatial-lai')


# LGRIP30

# General
dst_fold = dest_root / 'geospatial' / category / 'lgrip30' / basin_id
dst_fold.mkdir(exist_ok=True)

# Main
file_src = "lgrip30_agriculture.tif"
file_dst = f"{basin_id}_{file_src}"
src = geo_folder / 'lgrip30' / 'raw' / file_src
dst = dst_fold / file_dst
shutil.copy(src, dst)
copied_files.append('geospatial-lgrip30')


# MERIT

# General
dst_fold = dest_root / 'geospatial' / category / 'merit' / basin_id
dst_fold.mkdir(exist_ok=True)

# Aspect
file_src = "merit_hydro_aspect.tif"
file_dst = f"{basin_id}_{file_src}"
src = geo_folder / 'merit' / 'aspect' / file_src
dst = dst_fold / file_dst
shutil.copy(src, dst)
copied_files.append('geospatial-merit-aspect')

# DEM
file_src = "merit_hydro_elv.tif"
file_dst = f"{basin_id}_{file_src}"
src = geo_folder / 'merit' / 'raw' / file_src
dst = dst_fold / file_dst
shutil.copy(src, dst)
copied_files.append('geospatial-merit-dem')

# Slope
file_src = "merit_hydro_slope.tif"
file_dst = f"{basin_id}_{file_src}"
src = geo_folder / 'merit' / 'slope' / file_src
dst = dst_fold / file_dst
shutil.copy(src, dst)
copied_files.append('geospatial-merit-slope')


# MODIS land

# General
dst_fold = dest_root / 'geospatial' / category / 'modis-land' / basin_id
dst_fold.mkdir(exist_ok=True)

# Mode file
file_src = "2001_2022_mode_MCD12Q1_LC_Type1.tif"
file_dst = f"{basin_id}_{file_src}"
src = geo_folder / 'modis_land' / 'raw' / file_src
dst = dst_fold / file_dst
shutil.copy(src, dst)
copied_files.append('geospatial-modis-land-mode')

# Annual files
for year in range(2001,2023):
    file_src = f"{year}0101_MCD12Q1_LC_Type1.tif"
    file_dst = f"{basin_id}_{file_src}"
    src = geo_folder / 'modis_land' / 'raw' / file_src
    dst = dst_fold / file_dst
    shutil.copy(src, dst)


# Pelletier

# General
dst_fold = dest_root / 'geospatial' / category / 'pelletier' / basin_id
dst_fold.mkdir(exist_ok=True)

# Files
src_files = ['average_soil_and_sedimentary-deposit_thickness.tif',
             'hill-slope_valley-bottom.tif',
             'land_cover_mask.tif',
             'upland_hill-slope_regolith_thickness.tif',
             'upland_hill-slope_soil_thickness.tif',
             'upland_valley-bottom_and_lowland_sedimentary_deposit_thickness.tif']
for file_src in src_files:
    file_dst = f"{basin_id}_{file_src}"
    src = geo_folder / 'pelletier' / 'raw' / file_src
    dst = dst_fold / file_dst
    shutil.copy(src, dst)


# soilgrids

# General
dst_fold = dest_root / 'geospatial' / category / 'soilgrids' / basin_id
dst_fold.mkdir(exist_ok=True)

sub_folders = ['bdod','cfvo','clay','conductivity','porosity','sand','silt','soc']
for sub_folder in sub_folders:
    src_fold = geo_folder / 'soilgrids' / 'raw' / sub_folder
    files_src = glob.glob(str(src_fold/'*.tif'))
    for file_src in files_src:
        file_name = os.path.basename(file_src)
        file_dst = f"{basin_id}_{file_name}"
        dst = dst_fold / file_dst
        shutil.copy(file_src, dst)


# WorldClim

# General
dst_fold = dest_root / 'geospatial' / category / 'worldclim-derived' / basin_id
dst_fold.mkdir(exist_ok=True)

# annual files
src_fold = geo_folder / 'worldclim' / 'raw' / 'annual'
files_src = ['aridity.tif','fracsnow.tif','pet_sum.tif','prec_sum.tif',
             'seasonality.tif','snow_sum.tif','t_avg.tif']
files_dst = [f"{basin_id}_wc2.1_30s_derived_annual_aridity.tif",
             f"{basin_id}_wc2.1_30s_derived_annual_snow_fraction.tif",
             f"{basin_id}_wc2.1_30s_derived_annual_pet_total.tif",
             f"{basin_id}_wc2.1_30s_derived_annual_prec_total.tif",
             f"{basin_id}_wc2.1_30s_derived_annual_seasonality.tif",
             f"{basin_id}_wc2.1_30s_derived_annual_snow_total.tif",
             f"{basin_id}_wc2.1_30s_derived_annual_mean_temp.tif"]

for file_src,file_dst in zip(files_src,files_dst):
    src = src_fold / file_src
    dst = dst_fold / file_dst
    shutil.copy(src, dst)

# Monthly files
sub_folders = ['aridity2','fracsnow2','pet','snow2']
new_names = ['aridity','snow_fraction','pet_total','snow_total']
for sub_folder,new_name in zip(sub_folders,new_names):
    src_fold = geo_folder / 'worldclim' / 'raw' / sub_folder
    for mm in range(1,13):
        file_src = f"wc2.1_30s_{sub_folder}_{mm:02}.tif"
        file_dst = f"{basin_id}_wc2.1_30s_derived_{new_name}_month_{mm:02}.tif"
        src = src_fold / file_src
        dst = dst_fold / file_dst
        shutil.copy(src, dst)


# -- CHECK
# Check if we missed anything
missing_elements = [item for item in expected_files if item not in copied_files]

# Assertion with error message
assert not missing_elements, f"{basin_id}: {missing_elements} not successfully copied"