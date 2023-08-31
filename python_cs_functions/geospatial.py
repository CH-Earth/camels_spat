# Code to download geospatial data (forcing files, parameter fields)
import cdsapi
from datetime import datetime, timedelta
import geopandas as gpd
import netCDF4 as nc4
import numpy as np
import math
import os
import time
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# --- Folder setup ---
def prepare_forcing_outputs(df,i,data_path):
    
    '''Prepares output folders for lumped and distributed forcing downloads outcomes'''
    
    from pathlib import Path
    
    # Get identifiers
    country = df.iloc[i].Country
    basin_id = df.iloc[i].Station_id
    full_id = country + '_' + basin_id
    
    # Construct the paths
    main_folder = Path(data_path) / 'basin_data' / (country + '_' + basin_id) / 'forcing' 
    lump_folder = main_folder / 'lumped'
    dist_folder = main_folder / 'distributed'
    raw_folder  = main_folder / 'raw'
    
    # Make the paths
    lump_folder.mkdir(parents=True, exist_ok=True)
    dist_folder.mkdir(parents=True, exist_ok=True)
    raw_folder.mkdir(parents=True, exist_ok=True)
   
    return raw_folder, lump_folder, dist_folder

# --- Handling of time periods for which we have flow records ---
def convert_to_date_only(time_string, original_format="%Y-%m-%d %H:%M:%S", target_format="%Y-%m-%d"):
    try:
        datetime_obj = datetime.strptime(time_string, original_format)
        date_only_string = datetime_obj.strftime(target_format)
        return date_only_string
    except ValueError:
        return None

def round_flow_obs_to_days(times):
    
    '''Takes two times ([time1,time2]) in 'YYYY-MM-DD hh:mm:ss' and rounds to days '''
    
    return [convert_to_date_only(times[0]),
            convert_to_date_only(times[1])]

def flow_obs_unavailable(df,country,station):
    
    '''Checks in the "unusable" dataframe if iv, dv or both are unavailable for a station'''
    
    missing = []
    for ix,row in df.iterrows():
        if row.Country != country:
            continue
        if row.Station_id == station:
            missing.append(row.Missing)
    
    return missing

def find_flow_obs_times_from_metadata(row,missing):
    
    '''Finds required data start and end times from flow observations in meta-data file'''
    
    # Start and end dates come as 'YYYY-MM-DD hh:mm:ss' strings so we can directly compare them
    # Source: https://stackoverflow.com/a/54987418
    
    # Shorthands
    iv_s = row.iv_flow_obs_availability_start
    iv_e = row.iv_flow_obs_availability_end
    dv_s = row.dv_flow_obs_availability_start
    dv_e = row.dv_flow_obs_availability_end
    
    # Missing data cases
    if ('iv' in missing) and ('dv' in missing):
        return []
    
    elif 'iv' in missing:
        times = [dv_s,dv_e]
        for time in times: 
            assert is_valid_date_format(time), f'{time} not in expected format'
        return [dv_s,dv_e]
    
    elif 'dv' in missing:
        times = [iv_s,iv_e]
        for time in times: 
            assert is_valid_date_format(time), f'{time} not in expected format'
        return [iv_s,iv_e]
    
    else:
        times = [iv_s,iv_e,dv_s,dv_e]
        for time in times: 
            assert is_valid_date_format(time), f'{time} not in expected format'
    return [min(iv_s,dv_s),max(iv_e,dv_e)]

def is_valid_date_format(date_string, date_format='%Y-%m-%d %H:%M:%S'):
    try:
        datetime.strptime(date_string, date_format)
        return True
    except ValueError:
        return False

# --- Define spatial extent of download domain ---
def find_shapefile_bounds(path):
    
    # Modified from: https://github.com/CH-Earth/CWARHM/blob/main/0_tools/ERA5_find_download_coordinates_from_shapefile.ipynb
    shp = gpd.read_file(path)
    
    return shp.total_bounds

def find_download_coords_from_bounds(coords, target='ERA5'):
    
    '''
    Determines download coordinates from shapefile bounds for a given data set.
    Assumes coodinates are an array: [lon_min, lat_min, lon_max, lat_max] (bottom-left, top-right).
    Returns separate lat and lon vectors.
    '''

    # Source: https://github.com/CH-Earth/CWARHM/blob/main/3a_forcing/1a_download_forcing/download_ERA5_pressureLevel_annual.ipynb   
    
    # Extract values
    lon = [coords[0],coords[2]]
    lat = [coords[1],coords[3]]
    
    if target == 'ERA5':
        
        # Round to ERA5 0.25 degree resolution
        rounded_lon = [math.floor(lon[0]*4)/4, math.ceil(lon[1]*4)/4]
        rounded_lat = [math.floor(lat[0]*4)/4, math.ceil(lat[1]*4)/4]

        # Find if we are still in the representative area of a different ERA5 grid cell
        if lat[0] > rounded_lat[0]+0.125:
            rounded_lat[0] += 0.25
        if lon[0] > rounded_lon[0]+0.125:
            rounded_lon[0] += 0.25
        if lat[1] < rounded_lat[1]-0.125:
            rounded_lat[1] -= 0.25
        if lon[1] < rounded_lon[1]-0.125:
            rounded_lon[1] -= 0.25
    
        # Make a download string ready for ERA5 (cdsapi) format
        dl_string = '{}/{}/{}/{}'.format(rounded_lat[1],rounded_lon[0],rounded_lat[0],rounded_lon[1])
    
    return dl_string, rounded_lat, rounded_lon

# --- ERA5 downloads ---
def download_era5_surface_level_data_to_netcdf(coordinates,date_s,date_e,file,retries_max=10):
    
    '''Downloads specified ERA5 surface level parameters for a specified bounding box'''
    # Pressure level variables: 

    if not os.path.isfile(file):
        # Make sure the connection is re-tried if it fails
        retries_cur = 1
        while retries_cur <= retries_max:
            try:
    
                # connect to Copernicus (requires .cdsapirc file in $HOME)
                c = cdsapi.Client()
            
                # specify and retrieve data
                c.retrieve('reanalysis-era5-single-levels', { # do not change this!
                           'product_type': 'reanalysis',
                           'format'      : 'netcdf',
                           'variable'    : [
                               'mean_surface_downward_long_wave_radiation_flux',
                               'mean_surface_net_long_wave_radiation_flux',
                               'mean_surface_downward_short_wave_radiation_flux',
                               'mean_surface_net_short_wave_radiation_flux',
                               'mean_total_precipitation_rate',
                               'surface_pressure',
                               'mean_potential_evaporation_rate',
                           ],
                           'date': f'{date_s}/{date_e}', # 'yyyy-mm-dd'
                           'time': '00/to/23/by/1',
                           'area': coordinates, # expected as [lat_max, lon_min, lat_min, lon_max], e.g. [51.75/-116.5/51.0/-115.5]
                           'grid': '0.25/0.25', # Latitude/longitude grid: east-west (longitude) and north-south resolution (latitude).
                    },
                    file) # file path and name

                # track progress
                print('Successfully downloaded ' + str(file))

            except:
                print('Error downloading ' + str(file) + ' on try ' + str(retries_cur))
                retries_cur += 1
                continue
            else:
                break
    return

def download_era5_pressure_level_data_to_netcdf(coordinates,date_s,date_e,file,retries_max=10):
    
    '''Downloads specified ERA5 pressure level parameters for a specified bounding box'''
    # Pressure level variables: https://confluence.ecmwf.int/pages/viewpage.action?pageId=82870405#ERA5:datadocumentation-Table9
    # MARS requests: https://confluence.ecmwf.int/display/UDOC/HRES%3A+Atmospheric+%28oper%29%2C+Model+level+%28ml%29%2C+Forecast+%28fc%29%3A+Guidelines+to+write+efficient+MARS+requests
    
    if not os.path.isfile(file):
        # Make sure the connection is re-tried if it fails
        retries_cur = 1
        while retries_cur <= retries_max:
            try:
    
                # connect to Copernicus (requires .cdsapirc file in $HOME)
                c = cdsapi.Client()
            
                # specify and retrieve data
                c.retrieve('reanalysis-era5-complete', {    # do not change this!
                           'class'   : 'ea',
                           'expver'  : '1',
                           'stream'  : 'oper',
                           'type'    : 'an',
                           'levtype' : 'ml',
                           'levelist': '137',
                           'param'   : '130/131/132/133', # i.e., Temperature, U and V wind, specific humidity
                           'date'    : f'{date_s}/to/{date_e}', # 'yyyy-mm-dd'
                           'time'    : '00/to/23/by/1', 
                           'area'    : coordinates, # expected as [lat_max, lon_min, lat_min, lon_max], e.g. [51.75/-116.5/51.0/-115.5]
                           'grid'    : '0.25/0.25', # Latitude/longitude grid: east-west (longitude) and north-south resolution (latitude).
                           'format'  : 'netcdf',
                    }, file)

                # track progress
                print('Successfully downloaded ' + str(file))

            except:
                print('Error downloading ' + str(file) + ' on try ' + str(retries_cur))
                retries_cur += 1
                continue
            else:
                break
    return

def download_era5_time_invariant_data_to_netcdf(coordinates,file,retries_max=10):
    
    '''Downloads all ERA5 time-invariant parameters for a specified bounding box'''
    # Time-invariants: https://confluence.ecmwf.int/pages/viewpage.action?pageId=82870405#ERA5:datadocumentation-Table1

    if not os.path.isfile(file):
        # Make sure the connection is re-tried if it fails
        retries_cur = 1
        while retries_cur <= retries_max:
            try:
    
                # connect to Copernicus (requires .cdsapirc file in $HOME)
                c = cdsapi.Client()
            
                # specify and retrieve data
                c.retrieve('reanalysis-era5-complete', {    # do not change this!
                           'stream' : 'oper',
                           'levtype': 'sf',
                           'param'  : '26/228007/27/28/29/30/43/74/129/160/161/162/163/172', # i.e., all time-invariant values
                           'date'   : '2023-01-01', # arbitrary date
                           'time'   : '00', # time-invariant data; no need to get more than a single time step
                           'area'   : coordinates, # expected as [lat_max, lon_min, lat_min, lon_max], e.g. [51.75/-116.5/51.0/-115.5]
                           'grid'   : '0.25/0.25', # Latitude/longitude grid: east-west (longitude) and north-south resolution (latitude).
                           'format' : 'netcdf',
                    }, file)

                # track progress
                print('Successfully downloaded ' + str(file))

            except:
                print('Error downloading ' + str(file) + ' on try ' + str(retries_cur))
                retries_cur += 1
                continue
            else:
                break
    return

def convert_start_and_end_dates_to_era5_download_lists(start,end):

    '''Takes two datetime.datetime(y,m,d,h,min) objects and returns two lists with start and end dates for ERA5 downloads at monthly intervals'''

    # Initiate the date we're current working with
    cur = start 
    
    # Initiate the lists
    start_l = []
    end_l = []
    
    # Loop over the dates, until we have the end date
    while cur < end:
        
        # Add to start list
        start_l.append(cur)
        
        # Figure out the index of the next month and if the year changes:
        tmp = cur + timedelta(days=31) # Add 31 days to current date to ensure we're in the next month, might also switch the year
        next_month = tmp.month #  Extract 'month' from this object
        next_year = tmp.year # Extract the year too. If we ticked over into a new year we need to track this, otherwise we never increment the year

        # Create the end-of-month date
        cur = cur.replace(year=next_year, month=next_month) - timedelta(days=1)
        
        # Ensure this does not step over our end date
        if cur >= end:
            cur = end

        # Add to end list
        end_l.append(cur)

        # Add 1 day to create the new start-of-month date
        cur = cur + timedelta(days=1)
        
    return start_l,end_l

# --- ERA5 processing ---
def extract_ERA5_subset(infile, outfile, coords):
    
    '''Subsets an existing ERA5 forcing file by coordinates "latmax / lonmin / latmin / lonmax"'''

    # Source: https://github.com/CH-Earth/CWARHM/blob/main/0_tools/ERA5_subset_forcing_file_by_lat_lon.py

    # Notes:
    # 1. Works for North American continent
    # 2. Works for two specific ERA5 file layouts

    # Assumptions
    # 1. Latitude = [-90,90]
    # 2. Longitude = [-180,180]
    
    # Split coordinates
    coords = coords.split('/') # split string
    coords = [float(value) for value in coords] # string to array
    latmax = coords[0]
    lonmin = coords[1]
    latmin = coords[2]
    lonmax = coords[3]
    
    with xr.open_dataset(infile) as ds:

        # Handle specific cases
        if (ds['longitude'] > 180).any(): # convert ds longitude form 0/360 to -180/180
            lon = ds['longitude'].values
            lon[lon > 180] = lon[lon > 180] - 360
            ds['longitude'] = lon
        
        if (ds['latitude'] > 90).any(): # convert ds longitude form 0/180 to -90/90
            lat = ds['latitude'].values
            lat[lat > 90] = lat[lat > 90] - 180
            ds['latitue'] = lat

        # Subset
        ds_sub = ds.sel(latitude = slice(latmax, latmin), longitude = slice(lonmin, lonmax))
        ds_sub.to_netcdf(outfile)
        ds_sub.close()

def compare_forcing_data_and_shape_extents(save_path, nc_file, shp_file, nc_var, nc_time=0):

    '''Plots forcing grid and shapefile on top of each other'''

    catchment = gpd.read_file(shp_file)
    with xr.open_dataset(nc_file) as forcing:
        ax = plt.subplot()   
        if (len(forcing['latitude']) == 1) or (len(forcing['longitude']) == 1):
            # Case 1: in at least one direction, ERA5 file covers only a single grid cell. Needs special plotting attention
            lon = forcing['longitude'].values
            lat = forcing['latitude'].values
            grid = Rectangle( (lon.min()-0.125,lat.min()-0.125), lon.max()-lon.min()+0.25, lat.max()-lat.min()+0.25 )
            ax.add_patch(grid)
        else:
            # Case 2: more then a single ERA5 cell in either direction
            forcing[nc_var].isel(time=nc_time).plot(ax=ax,cbar_kwargs={'shrink': 0.95})
        catchment.plot(color='None', edgecolor='k', ax=ax);
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()

def compare_era5_netcdf_dimensions(surf_file, pres_file):

    '''Opens two netCDF files containing ERA5 data, reads any variables in _vars_ and checks if contents all match'''

    # Extract space and time info for checks
    with nc4.Dataset(pres_file) as src1, nc4.Dataset(surf_file) as src2: # implicitly closes files
        pres_lat = src1.variables['latitude'][:]
        pres_lon = src1.variables['longitude'][:]
        pres_time = src1.variables['time'][:]
        surf_lat = src2.variables['latitude'][:]
        surf_lon = src2.variables['longitude'][:]
        surf_time = src2.variables['time'][:]

     # Update the pressure level coordinates
    pres_lat[pres_lat > 90] = pres_lat[pres_lat > 90] - 180
    pres_lon[pres_lon > 180] = pres_lon[pres_lon > 180] - 360

    # Perform the check
    if not all( [all(pres_lat == surf_lat), all(pres_lon == surf_lon), all(pres_time == surf_time)] ):
        err_txt = (f'ERROR: Dimension mismatch while merging {pressure_file} and {surface_file}. '
                    'Check latitude, longitude and time dimensions in both files. Continuing with next files.')
        return False,err_txt
    else:
        return True,'Variables match'

def merge_era5_surface_and_pressure_files(surf_file, pres_file, dest_file):

    '''Merges two ERA5 data files into dest_file'''

    # Transfer definitions
    dimensions_surf_transfer = ['longitude','latitude','time']
    variables_surf_transfer = ['msdwlwrf','msnlwrf','msdwswrf','msnswrf','mtpr','sp','mper']
    variables_pres_transfer = ['t','q','u','v']
    attr_names_expected = ['scale_factor','add_offset','_FillValue','missing_value','units','long_name','standard_name'] # these are the attributes we think each .nc variable has             
    loop_attr_copy_these = ['units','long_name','standard_name'] # we will define new values for _FillValue and missing_value when writing the .nc variables' attributes

    with nc4.Dataset(pres_file) as src1, nc4.Dataset(surf_file) as src2, nc4.Dataset(dest_file, 'w') as dest: # implicitly closes files

        # Set general attributes
        dest.setncattr('History','Created ' + time.ctime(time.time()))
        dest.setncattr('Reason','Merging separate ERA5 download files')
        dest.setncattr('Source','github.com/cH-Earth/camels_spat')

        # Copy meta data from the two source files
        if src1.getncattr('Conventions') == 'CF-1.6' and src2.getncattr('Conventions') == 'CF-1.6':
            dest.setncattr('Conventions','CF-1.6')
        dest.setncattr('History_source_file_1',src1.getncattr('history'))
        dest.setncattr('History_source_file_2',src2.getncattr('history'))

        # --- Dimensions: latitude, longitude, time
        # NOTE: we can use the lat/lon from the surface file (src2), because those are already in proper units. 
        # If there is a mismatch between surface and pressure we shouldn't have reached this point at all due to the earlier check
        for name, dimension in src2.dimensions.items():
            if dimension.isunlimited():
                dest.createDimension( name, None)
            else:
                dest.createDimension( name, len(dimension))

        # --- Get the surface level dimension variables (lat, lon, time)
        for name, variable in src2.variables.items():
    
            # Transfer lat, long and time variables because these don't have scaling factors
            if name in dimensions_surf_transfer:
                dest.createVariable(name, variable.datatype, variable.dimensions, fill_value = -999)
                dest[name].setncatts(src1[name].__dict__)
                dest.variables[name][:] = src2.variables[name][:]

        # ---  Transfer the surface level data first, for no particular reason
        # This should contain surface pressure (sp), downward longwave (msdwlwrf), downward shortwave (msdwswrf) and precipitation (mtpr)
        for name, variable in src2.variables.items():

            # Check that we are only using the names we expect, and thus the names for which we have the required code ready
            if name in variables_surf_transfer:
        
                # 0. Reset the dictionary that we keep attribute values in
                loop_attr_source_values = {name: 'n/a' for name in attr_names_expected}
        
                # 1a. Get the values of this variable from the source (this automatically applies scaling and offset)
                loop_val = variable[:]
        
                # 1b. Get the attributes for this variable from source
                for attrname in variable.ncattrs():
                    loop_attr_source_values[attrname] = variable.getncattr(attrname)
               
                # 2. Create the .nc variable 
                # Inputs: variable name as needed by SUMMA; data type: 'float'; dimensions; no need for fill value, because thevariable gets populated in this same script
                dest.createVariable(name, 'f4', ('time','latitude','longitude'), fill_value = False, zlib=True, shuffle=True)
        
                # 3a. Select the attributes we want to copy for this variable, based on the dictionary defined before the loop starts
                loop_attr_copy_values = {use_this: loop_attr_source_values[use_this] for use_this in loop_attr_copy_these}
        
                # 3b. Copy the attributes FIRST, so we don't run into any scaling/offset issues
                dest[name].setncattr('missing_value',-999)
                dest[name].setncatts(loop_attr_copy_values)
        
                # 3c. Copy the data SECOND
                dest[name][:] = loop_val

        # --- Transfer the pressure level variables next, using the same procedure as above
        for name, variable in src1.variables.items():
            if name in variables_pres_transfer:
        
                # 0. Reset the dictionary that we keep attribute values in
                loop_attr_source_values = {name: 'n/a' for name in attr_names_expected}
        
                # 1a. Get the values of this variable from the source (this automatically applies scaling and offset)
                loop_val = variable[:] 
        
                # 1b. Get the attributes for this variable from source
                for attrname in variable.ncattrs():
                    loop_attr_source_values[attrname] = variable.getncattr(attrname)
               
                # 2a. Create the .nc variable with the proper SUMMA name
                # Inputs: variable name as needed by SUMMA; data type: 'float'; dimensions; no need for fill value, because thevariable gets populated in this same script
                dest.createVariable(name, 'f4', ('time','latitude','longitude'), fill_value = False, zlib=True, shuffle=True)
        
                # 3a. Select the attributes we want to copy for this variable, based on the dictionary defined before the loop starts
                loop_attr_copy_values = {use_this: loop_attr_source_values[use_this] for use_this in loop_attr_copy_these}
        
                # 3b. Copy the attributes FIRST, so we don't run into any scaling/offset issues
                dest[name].setncattr('missing_value',-999)
                dest[name].setncatts(loop_attr_copy_values)
        
                # 3c. Copy the data SECOND
                dest[name][:] = loop_val

# --- ERA5 derived variables
def add_derived_variables(file):
    
    '''Adds various derived meteorological variables to a netCDF file with existing ERA5 forcing data'''

    with nc4.Dataset(file, 'r+') as f:
        f = derive_reflected_shortwave_radiation(f)
        f = derive_net_radiation(f)
        f = derive_vapor_pressure(f)
        f = derive_relative_humidity(f)
        f = derive_mean_wind_speed(f)

def make_nc_variable(file, name, units, values, long_name='', standard_name='', history=''):

    '''Makes a netcdf4 variable with dimensions (time,latitude,longitude) in file'''

    # Assumptions
    # - Dimension in input file are 'time', 'latitude', 'longitude'
    # - No fill value specified
    # - Compression options zlib=True, shuffle=True are active

    # 1. Create the .nc variable 
    file.createVariable(name,'f4',('time','latitude','longitude'), fill_value = False, zlib=True, shuffle=True)

    # 2. Set the attributes FIRST, so we don't run into any scaling/offset issues
    file[name].setncattr('missing_value',-999)
    file[name].setncattr('units',units)
    if long_name: file[name].setncattr('long_name',long_name)
    if standard_name: file[name].setncattr('standard_name', standard_name)

    # 3. Copy the data SECOND
    file[name][:] = values

    # 4. Update history
    if history:
        current_history = file.getncattr('History')
        new_history = f'{current_history}{history}'
        file.setncattr('History', new_history)
    
    return file

def derive_reflected_shortwave_radiation(file, incoming_shortwave='msdwswrf', net_shortwave='msnswrf', new_name='reflected_sw'):

    '''Adds reflected shortwave radiation to a netCDF4 file that cantains incoming and net shortwave radiation'''

    # 1. Get the values of this variable from the source (this automatically applies scaling and offset)
    ds = file.variables[incoming_shortwave][:]
    ns = file.variables[net_shortwave][:]

    # 2. Compute the values 
    values_new = ds - ns
    
    # 3. Define other inputs
    unit_new = 'W m**-2'
    long_name = 'Mean surface reflected shortwave radiation flux, computed from ERA5 mean surface downward short-wave radiation flux and mean surface net short-wave radiation flux'
    new_history = f' On {time.ctime(time.time())}: derive_reflected_shortwave_radiation().'

    # 4. Make the variable
    file = make_nc_variable(file, new_name, unit_new, values_new, long_name=long_name, history=new_history)
    
    return file

def derive_net_radiation(file, net_shortwave='msnswrf', net_longwave='msnlwrf', new_name='net_radiation'):

    '''Adds net radiation to a netCDF4 file that contains net shortwave and net longwave radiation'''

    # 1. Get the values of this variable from the source (this automatically applies scaling and offset)
    ns = file.variables[net_shortwave][:]
    nl = file.variables[net_longwave][:]

    # 2. Compute the values 
    values_new = ns + nl
    
    # 3. Define other inputs
    unit_new = 'W m**-2'
    long_name = 'Mean surface net radiation flux, computed from ERA5 mean surface net short-wave radiation flux and mean surface net long-wave radiation flux'
    new_history = f' On {time.ctime(time.time())}: derive_net_radiation().'

    # 4. Make the variable
    file = make_nc_variable(file, new_name, unit_new, values_new, long_name=long_name, history=new_history)
    
    return file

def derive_relative_humidity(file, air_temperature='t', vapor_pressure='e', new_name='rh'):

    '''Adds relative humidity to a netCDF4 file containing vapor pressure and air temperature'''

    # 0. Constants
    Rv = 461 # [J K-1 kg-1]
    T0 = 273.15 # [K]
    e0 = 0.6113 # [kPa]
    Lv = 2.5*10**6 # [J kg-1], latent heat of vaporization for liquid water
    Ld = 2.83*10**6 # [J kg-1], latent heat of deposition for ice

    # 1. Get the values of this variable from the source (this automatically applies scaling and offset)
    e = file.variables[vapor_pressure][:]
    T = file.variables[air_temperature][:]

    # 2. Compute the values
    is_ice = T < 273.15
    is_liq = T >= 273.15    
    es = T*0 # initialize the variable
    es[is_ice] = e0*np.exp(Ld/Rv * (1/T0 - 1/T[is_ice]))
    es[is_liq] = e0*np.exp(Lv/Rv * (1/T0 - 1/T[is_liq]))   
    values_new = e / es
    
    # 3. Define other inputs
    unit_new = 'kPa kPa**-1'
    long_name = 'relative humidity, computed from ERA5 temperature and derived vapor pressure'
    new_history = f' On {time.ctime(time.time())}: derive_relative_humidity().'

    # 4. Make the variable
    file = make_nc_variable(file, new_name, unit_new, values_new, long_name=long_name, history=new_history)
    return file

def derive_vapor_pressure(file, air_pressure='sp', specific_humidity='q', Pa_to_kPa=True, new_name='e', approximate=False):

    '''Adds vapor pressure to a netCDF4 file containing air pressure and specific humidity'''

    # Assumptions:
    # - New variable's attribute values are hard-coded
    
    # 0. Constants
    Rd = 2.871*10**(-4) # [kPa K-1 m3 kg-1]
    Rv = 4.61*10**(-4) # [kPa K-1 m3 kg-1]
    epsilon = Rd/Rv # Should be ~ 0.622
    
    # 1. Get the values of this variable from the source (this automatically applies scaling and offset)
    P = file.variables[air_pressure][:]
    q = file.variables[specific_humidity][:]

    # 2. Convert air pressure units if flagged
    if Pa_to_kPa: P = P / 1000

    # 3. Compute the values
    values_new = -1*(q*P)/(q*epsilon - q - epsilon)
    if approximate: values_new = q*P/(epsilon)

    # 4. Define other inputs
    unit_new = 'kPa'
    long_name = 'vapor pressure, computed from ERA5 specific humidity and surface pressure'
    new_history = f' On {time.ctime(time.time())}: derive_vapor_pressure().'

    # 5. Make the variable
    file = make_nc_variable(file, new_name, unit_new, values_new, long_name=long_name, history=new_history)
    return file

def derive_mean_wind_speed(file, var_1='u', var_2='v', new_name='w'):
    
    '''Adds mean wind speed to a netCDF4 file containing U and V wind direction data'''

    # Assumptions:
    # - New variable's attribute values are hard-coded
    
    # 1. Get the values of this variable from the source (this automatically applies scaling and offset)
    values_1 = file.variables[var_1][:]
    values_2 = file.variables[var_2][:]

    # 2. Create the variable attribute 'units' from the source data. This lets us check if the source units match (they should match)
    unit_1 = file.variables[var_1].getncattr('units')
    unit_2 = file.variables[var_2].getncattr('units')
    assert unit_1 == unit_2, f'WARNING: units of source variables do not match: variable {var1} in {unit_1}, varialbe {var_2} in {unit_2}'
    unit_new = '(({})**2 + ({})**2)**0.5'.format(unit_1,unit_2)

    # 3. Compute the values
    values_new = ((values_1**2)+(values_2**2))**0.5
    
    # 4. Define other inputs
    long_name = 'mean wind speed at the measurement height, computed from ERA5 U and V component of wind'
    standard_name = 'wind_speed'
    new_history = f'. On {time.ctime(time.time())}: derive_mean_wind_speed().'

    # 5. Make the variable
    file = make_nc_variable(file, new_name, unit_new, values_new, long_name=long_name, standard_name=standard_name, history=new_history)
    return file