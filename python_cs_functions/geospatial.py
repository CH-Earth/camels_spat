# Code to download geospatial data (forcing files, parameter fields)
import cdsapi
from datetime import datetime, timedelta
import geopandas as gpd
import math
import os

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

