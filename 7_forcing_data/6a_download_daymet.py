import sys
from pathlib import Path
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

# --- Config handling
# Specify where the config file can be found
config_file = '../0_config/config.txt'

# Temporary download path
temp_folder = Path( cs.read_from_config(config_file, 'temp_path') )

# --- Daymet download location
temp_daymet_folder = temp_folder / 'daymet'
temp_daymet_folder.mkdir(parents=True, exist_ok=True)

# --- Loop over the years and variables, and download the files
# https://daac.ornl.gov/daacdata/daymet/Daymet_Daily_V4R1/data/daymet_v4_daily_na_dayl_1980.nc
src_url = 'https://daac.ornl.gov/daacdata/daymet/Daymet_Daily_V4R1/data/'

# Single threaded download
#for year in range(1980, 2024):
#    for var in ['prcp', 'tmax', 'tmin', 'srad', 'vp', 'dayl']:
#        file_name = f'daymet_v4_daily_na_{var}_{year}.nc'
#        file_url = src_url + f'{file_name}'
#        cs.download_url_into_folder(file_url, temp_daymet_folder, retries_max=10, requests_kwargs={}, overwrite=False)

import concurrent.futures

def download_file(year, var):
    file_name = f'daymet_v4_daily_na_{var}_{year}.nc'
    file_url = src_url + f'{file_name}'
    cs.download_url_into_folder(file_url, 
                                temp_daymet_folder, 
                                retries_max=1,
                                overwrite=False)

num_cores = 10
with concurrent.futures.ThreadPoolExecutor(max_workers=num_cores) as executor:
    futures = []
    for year in range(1980, 2024):
        for var in ['dayl', 'prcp', 'tmax', 'tmin', 'srad', 'vp']:
            futures.append(executor.submit(download_file, year, var))
    
    for future in concurrent.futures.as_completed(futures):
        try:
            future.result()
        except Exception as e:
            print(f"An error occurred: {e}")
