'''Contains functions to perform downloads.'''

def download_url_into_folder(url,folder, retries_max=10, requests_kwargs={}):
    
    import shutil
    import requests
          
    # Extract the filename from the URL
    file_name = url.split('/')[-1].strip() # Get the last part of the url, strip whitespace and characters
    
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
def prepare_flow_download_outputs(df,i,data_path):
    
    '''Prepares output folder and file paths for flow observation downloads'''
    
    from pathlib import Path
    
    # Get identifiers
    country = df.iloc[i].Country
    basin_id = df.iloc[i].Station_id
    full_id = country + '_' + basin_id
    
    # Construct the paths
    main_folder = Path(data_path) / 'basin_data' / (country + '_' + basin_id) / 'observations'
    
    # Make the paths
    main_folder.mkdir(parents=True, exist_ok=True)
    
    # Make the output file paths
    raw_file  = main_folder / (full_id + '_raw.txt') # To contain server response
    data_file = main_folder / (full_id + '_flow_observations_raw.csv') # To contain data only
    meta_file = main_folder / (full_id + '_header.txt') # To contain metadata/header info only
    hour_file = main_folder / (full_id + '_flow_observations_hourly.nc') # Final file, to contain hourly data and metadata
    
    return basin_id, full_id, raw_file, data_file, meta_file, hour_file