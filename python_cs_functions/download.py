'''Contains functions to perform downloads.'''

import time
import shutil
import requests
import pandas as pd
import urllib.request
from pathlib import Path

def download_url_into_folder(url,folder, retries_max=10, requests_kwargs={}, overwrite=False):
    
    # Extract the filename from the URL
    file_name = url.split('/')[-1].strip() # Get the last part of the url, strip whitespace and characters
    
    # Early exit
    if Path(folder/file_name).is_file():
        if not overwrite:
            print(f'File {folder/file_name} exists and download_url_into_folder() argument overwrite is False. Skipping file.')
            return
    
    # Make sure the connection is re-tried if it fails
    retries_cur = 1
    while retries_cur <= retries_max:
        try: 

            # Send a HTTP request to the server and save the HTTP response in a response object called response
            # kwargs:
            # - stream: (True/False) ensures that only response headers are downloaded initially 
            # - header: (dict)       specifies details about the request, such as User-Agent
            # - auth  : (tuple)      user name, password
            with requests.get(url.strip(), stream=True, **requests_kwargs) as response:

                # Decode the response
                response.raw.decode_content = True
                content = response.raw

                # Write to file
                with open(folder / file_name, 'wb') as data:
                    shutil.copyfileobj(content, data)
                    
                # print a completion message
                print('Successfully downloaded ' + url)
                    
        except Exception as e:
            print('Error downloading ' + file_url + ' on try ' + str(retries_cur) + ' with error: ' + str(e))
            retries_cur += 1
            continue
        else:
            break           
    
    return

# File paths and names
# ------------------------------------------------------------------------------------------------------------------------
def prepare_flow_download_outputs(df,i,data_path,time='hourly'):
    
    '''Prepares output folder and file paths for flow observation downloads'''
      
    # Get identifiers
    country = df.iloc[i].Country
    basin_id = df.iloc[i].Station_id
    full_id = country + '_' + basin_id
    
    # Construct the paths
    main_folder = Path(data_path) / 'basin_data' / (country + '_' + basin_id) / 'observations'
    
    # Make the paths
    main_folder.mkdir(parents=True, exist_ok=True)
    
    # Make the output file paths
    raw_file  = main_folder / (full_id + '_{}_raw.txt'.format(time)) # To contain server response
    data_file = main_folder / (full_id + '_{}_raw_flow_observations.csv'.format(time)) # To contain data only
    head_file = main_folder / (full_id + '_{}_raw_header.txt'.format(time)) # To contain metadata/header info only
    nc_file = main_folder / (full_id + '_{}_flow_observations.nc'.format(time)) # Final file, to contain hourly data and metadata
    
    return basin_id, full_id, raw_file, data_file, head_file, nc_file

# USGS flow downloads
# ------------------------------------------------------------
def download_usgs_values(main_url, site, time_s, time_e, var, raw_path, dnf):
    
    '''Downloads Instananeous or Daily Values data from USGS'''
    
    # Construct the download URL
    url = f'{main_url}?format=rdb&sites={site}&startDT={time_s}&endDT={time_e}&parameterCd={var}&siteStatus=all'
        
    # Download the URL to a temporary location
    urllib.request.urlretrieve(url, raw_path)      
        
    # Checks
    df = pd.read_csv(raw_path, delimiter='\t', comment='#', low_memory=False) # skip comments (#); low_mem prevents mixed datatype warning 
    if len(df) < 3: # Sites with no data still have a 2-line df. 1st: data format. 2nd: NaNs
        print(f'No data downloaded for {site}')
        dnf.append(site)
    else:
        print(f'Completed {main_url} for {site}')

    return dnf

# WSC flow downloads
# ------------------------------------------------------------
def download_wsc_values(main_url, site, time_s, time_e, var, raw_path, dnf):
    
    '''Downloads real-time values from WSC'''
    
    # Construct the download URL
    url = f'{main_url}inline?stations[]={site}&parameters[]={var}&start_date={time_s}%2000:00:00&end_date={time_e}%2023:59:59'
    
    # Read the URL into a dataframe
    tmp = pd.read_csv(url, index_col=[1], parse_dates=True)
    tmp.to_csv(raw_path)
    
    # Checks
    if len(tmp) == 0: # Sites with no data only have a header, and that doesn't count towards dataframe length
        print(f'No data downloaded for {site}')
        dnf.append(site)
    else:
        print(f'Completed {site}')
    
    return dnf