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