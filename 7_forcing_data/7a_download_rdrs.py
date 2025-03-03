import os
import sys
import concurrent.futures
from pathlib import Path
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

#### Config handling
# Specify where the config file can be found
config_file = '../0_config/config.txt'

# Temporary download path
temp_folder = Path( cs.read_from_config(config_file, 'temp_path') )
temp_folder = temp_folder / 'rdrs' # Overwrite with actual location, so we can run a few things side-by-side
temp_folder.mkdir(parents=True, exist_ok=True)

# Download URL
rdrs_url = cs.read_from_config(config_file, 'rdrs_url')
rdrs_vars = cs.read_from_config(config_file, 'rdrs_vars')
rdrs_vars_to_keep = rdrs_vars.split(',')

#### Processing
# We know we want RDRS covers only 1980-2018 and we want all of that, 
# so we can simply find all files on the webpage and get those
folder_urls = cs.find_folders_on_webpage(rdrs_url,product='rdrs')
file_urls = []
for folder_url in folder_urls:
    file_urls.append(cs.find_file_urls_in_webpage_folders([folder_url], extension='.nc'))
file_urls = [item for sublist in file_urls for item in sublist] # flatten the list

# Define the processing function
# RDRS containa more variables than we want, so we remove the superfluous ones immediately
# after download to save some space.

# Variable list
# Precipitation (model)  [m]       RDRS_v2.1_P_PR0_SFC
# Precipitation (CaPA)   [m]       RDRS_v2.1_A_PR0_SFC
# Air temperature        [C]       RDRS_v2.1_P_TT_1.5m
# Downward shortwave rad [W m-2]   RDRS_v2.1_P_FB_SFC
# Downward longwave rad  [W m-2]   RDRS_v2.1_P_FI_SFC
# Relative humidity      [-]       RDRS_v2.1_P_HR_1.5m
# Specific humidty       [kg kg-1] RDRS_v2.1_P_HU_1.5m
# Surface pressure       [mb]      RDRS_v2.1_P_P0_SFC
# U wind component       [kts]     RDRS_v2.1_P_UUC_10m
# V wind component       [kts]     RDRS_v2.1_P_VVC_10m
# Wind speed             [kts]     RDRS_v2.1_P_UVC_10m
# Geopotential height    [dam]     RDRS_v2.1_P_GZ_10000

# Dev
#file_url = file_urls[0] # only download the first file for now
# -- end dev

def download_rdrs(file_url):
    cs.download_url_into_folder(file_url, temp_folder)
    file_name = file_url.split('/')[-1].strip()
    #print(f'Downloaded {file_name}')
    return

# Number of threads or processes to use
num_workers = 10
with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
    results = list( executor.map(download_rdrs, file_urls) )


# Subset variables
# We'll do this after using datatool to subset to basins
#def process_rdrs(file_url):
#    file_name = file_url.split('/')[-1].strip()
#    file_in = str(temp_folder / file_name) # str() for next replace and use with nc4
#    file_out = file_in.replace('12.nc','.nc') # we know all source files have the '12'[h] in the name
#    cs.remove_vars_from_rdrs_download(file_in,file_out,rdrs_vars_to_keep)
#    os.remove(file_in) # remove the original file
#    print(f'Processed {file_name}')
#    return

# DEV check - all True in test
#import xarray as xr
#src = xr.open_dataset(file_in)
#des = xr.open_dataset(file_out)
#for var in rdrs_vars_to_keep:
#    flag = (src[var] == des[var]).all().values
#    print(f'{var} match: {flag}')

# This takes too long.
    #for file_url in sub_list:
    #    process_rdrs(file_url)