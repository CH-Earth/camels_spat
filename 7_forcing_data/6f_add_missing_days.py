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
    print("Usage: python 6f_add_missing_days.py <array_index>")
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

# Get the variables we need
daymet_vars = cs.read_from_config(config_file, 'daymet_vars').split(',')

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
daymet_files = list(temp_daymet_folder.glob('daymet_*_stillLeapYearGaps.nc')) 
daymet_files.sort() # Doesn't matter but feels cleaner

# --- Processing
# Copy over all the files into the the main directory
# We'll overwrite the leap years in the next step
for file in daymet_files:
    # clean-up: remove the incorrectly named earlier copies from raw_fold (2024-06-11), if the file exists at all
    if (raw_fold / file.name).exists():
        os.remove(raw_fold / file.name)
    # drop the _stillLeapYearGaps part of the name
    new_name = file.name.replace('_stillLeapYearGaps', '')
    # copy the file with the old name to the new location with the new name
    shutil.copy(file, raw_fold / new_name)
    

# Open the first file and extract the encoding for each variable
data = xr.open_dataset(daymet_files[0])
encoding = {var: data[var].encoding for var in data.data_vars}
encoding['time'] = data['time'].encoding
data.close()

# Set fill values for yearday, lambert_conformal_conic and apply the same
# fillValue to pet as used in other meteo variables
encoding['yearday'] = {'_FillValue': -9999, 'missing_value': -9999}
encoding['lambert_conformal_conic'] = {'_FillValue': -9999.0, 'missing_value': -9999.0}
encoding['pet'] = {'_FillValue': -9999.0, 'missing_value': -9999.0}

# Now loop over the files pair-wise, and update the leap years
# We do this pairwise to avoid loading all the data at once
# We don't use open_mfdataset because ...
for file1,file2 in zip(daymet_files[:-1], daymet_files[1:]):
    
    # check if file1 is a leap year
    year = int(file1.name.split('_')[1]) # daymet_1980_stillLeapYearGaps.nc
    leap = (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)
    if not leap: continue # nothing to process here
    
    # Load the data
    data = xr.open_mfdataset([file1,file2])

    # Update the time index: add the missing day
    complete_time_index = pd.date_range(start=data['time'].min().values, 
                                        end=data['time'].max().values, freq='D')
    data = data.reindex(time=complete_time_index) # This will add NaNs for missing days
    data = data.chunk(dict(time=-1)) # Originally we had time chunks  of size 1, and this doesn't work with interpolate. Now we have size [long]

    # Interpolate the data variables
    for var in daymet_vars:
        data[var] = data[var].interpolate_na(dim='time', method='linear')

    # Update the time_bnds variable
    assert pd.isna(data['time_bnds'].isel(time=365)).all(), 'time_bnds not empty in expected location'
    data['time_bnds'].loc[f'{year}-12-31T12:00:00'] = data['time_bnds'].isel(time=364) + pd.Timedelta(days=1)

    # Update the yearday variable
    assert pd.isna(data['yearday'].isel(time=365)), 'yearday not empty in expected location'
    data['yearday'].loc[f'{year}-12-31T12:00:00'] = 366

    # Remove all the incorrectly added lambert_conformal_conic values
    llc_var = data['lambert_conformal_conic'].isel(time=0)
    llc_var = llc_var.drop_vars('time').squeeze() # Get rid of all the time stuff
    data['lambert_conformal_conic'] = llc_var

    # Apply the encoding where needed
    for var,enc in encoding.items():
        if data[var].encoding == {}:
            data[var].encoding = enc
        # Specify a fillValue for these variables specifically,
        # because it's not there in existing files and leads
        # to warnings during saving.
        if var == 'yearday' and data[var].encoding != {}:
            data[var].encoding['_FillValue'] = -9999
        if var == 'lambert_conformal_conic' and data[var].encoding != {}: 
            data[var].encoding['_FillValue'] = -9999.0

    # Save the first year to file
    new_name = file1.name.replace('_stillLeapYearGaps', '')
    data_y1 = data.sel(time=f'{year}')
    data_y1.to_netcdf(raw_fold / new_name)

'''
DEV NOTES

# Attempt 2: open_mfdataset
file1 = daymet_files[0]
file2 = daymet_files[1]

data = xr.open_mfdataset([file1,file2])
complete_time_index = pd.date_range(start=data['time'].min().values, 
                                        end=data['time'].max().values, freq='D')
data = data.reindex(time=complete_time_index) # This will add NaNs for missing days
data = data.chunk(dict(time=-1)) # Originally we had time chunks  of size 1, and this doesn't work with interpolate. Nopw we have size [long]

# Interpolate the data variables
for var in daymet_vars:
    data[var] = data[var].interpolate_na(dim='time', method='linear')

# Update the time_bnds variable
assert pd.isna(data['time_bnds'].isel(time=365)).all(), 'time_bnds not empty in expected location'
data['time_bnds'].loc[f'{year}-12-31T12:00:00'] = data['time_bnds'].isel(time=364) + pd.Timedelta(days=1)

# Update the yearday variable
assert pd.isna(data['yearday'].isel(time=365)), 'yearday not empty in expected location'
data['yearday'].loc[f'{year}-12-31T12:00:00'] = 366

# Remove all the incorrectly added lambert_conformal_conic values
llc_var = data['lambert_conformal_conic'].isel(time=0)
llc_var = llc_var.drop_vars('time').squeeze() # Get rid of all the time stuff
data['lambert_conformal_conic'] = llc_var

# Apply the encoding where needed
for var,enc in encoding.items():
    if data[var].encoding == {}:
        data[var].encoding = enc

# Save the first year to file
data_y1 = data.sel(time=f'{year}')
data_y1.to_netcdf(raw_fold / file1.name)
'''


'''
Attempt 1: xr.open_dataset()
# pro: works
# con: needs more memory, I think

# Load the data
data = xr.open_dataset(file1)
data = data.merge(xr.open_dataset(file2))

'''


'''
OLD STUFF: initial attempts to do everything in one file
This doesn't seem too feasible given that we need to run this
through EASYMORE later. More memory efficient to do with smaller files.

# Merge all of these into a single xarray dataset
# We need these arguments (combine, compat, coords) to avoid a small issue in some files
#  where the latitude values are slightly different for god knows what reason
data = xr.open_mfdataset(daymet_files, combine='by_coords', compat='override', coords='minimal')
#ds_list = [xr.open_dataset(f) for f in daymet_files] 
#data = xr.concat(ds_list, dim='time', coords='minimal', compat='override')

# Update the time index
complete_time_index = pd.date_range(start=data['time'].min().values, 
                                    end=data['time'].max().values, freq='D')
data = data.reindex(time=complete_time_index) # This will add NaNs for missing days
data = data.chunk(dict(time=-1)) # Originally we had time chunks  of size 1, and this doesn't work with interpolate. Nopw we have size [long]
data = data.interpolate_na(dim='time', method='linear') # This is the default but clearer this way
data = data.chunk(dict(time=1)) # back to original chunk size

# Ensure the time_bnds variable has encoding specified
if data['time_bnds'].encoding == {}:
    # Based on; https://github.com/pydata/xarray/pull/2965#issuecomment-493518414
    #data['time_bnds'].encoding = {'dtype': data['time'].encoding['dtype'], 
    #                              'units': data['time'].encoding['units'], 
    #                              'calendar': data['time'].encoding['calendar']}
    # This is what the encoding is in a regularly opened file: open_mfdataset strips the encoding
    data['time_bnds'].encoding = {'dtype': 'float32', 
                         'zlib': True, 'shuffle': False, 'complevel': 4, 
                         'fletcher32': False, 'contiguous': False, 
                         'chunksizes': (1, 2), 'preferred_chunks': {'time': 1, 'nv': 2}, 
                         'source': '/gpfs/tp/gwf/gwf_cmt/wknoben/camels_spat/camels-spat-data/basin_data/CAN_01AD002/forcing/raw/daymet_year/CAN_01AD002_subset_daymet_v4_daily_na_dayl_1980.nc', 
                         'original_shape': (365, 2), 
                         'units': 'days since 1950-01-01 00:00:00', 'calendar': 'standard'}

# Save individual years to file
years_grouped = data.groupby('time.year')
for year, year_ds in years_grouped:
    output_path = raw_fold / f'daymet_{year}.nc'
    year_ds.to_netcdf(output_path)
    print(f'Saved year {year} to {output_path}')

# Alternative:
years, datasets = zip(*data.groupby("time.year"))
output_paths = [raw_fold / f'daymet_{y}.nc' for y in years]
xr.save_mfdataset(datasets, paths)


data = xr.open_mfdataset(daymet_files[39:41])
data = xr.open_mfdataset(daymet_files[40:42], combine='by_coords', compat='override', coords='minimal')


d1 = xr.open_dataset(daymet_files[0]) # dayl
d2 = xr.open_dataset(daymet_files[42]) # pet

ds_list = [xr.open_dataset(f) for f in daymet_files] 
data = xr.concat(ds_list, dim='time', coords='minimal', compat='override')

# Subset the file list by year
year = 1980
ds_list = [xr.open_dataset(f) for f in daymet_files if f'{year}' in str(f) or f'{year+1}' in str(f)]
data = xr.concat(ds_list, dim='time', coords='minimal', compat='override')

data_y1 = []
for ds in [xr.open_dataset(f) for f in daymet_files if f'{year}' in str(f)]:
    if data_y1 == []:
        data_y1 = ds
    else:
        for var in ds.data_vars:
            if var in ['dayl', 'prcp', 'srad', 'tmax', 'tmin', 'vp', 'pet']:
                print(f'Adding {var}')
                data_y1[var] = ds[var]

data_y2 = []
for ds in [xr.open_dataset(f) for f in daymet_files if f'{year+1}' in str(f)]:
    if data_y2 == []:
        data_y2 = ds
    else:
        for var in ds.data_vars:
            if var in ['dayl', 'prcp', 'srad', 'tmax', 'tmin', 'vp', 'pet']:
                data_y2[var] = ds[var]
'''