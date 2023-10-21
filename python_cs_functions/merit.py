# Code adapted from CWARHM toolbox (Knoben et al., 2022a; Knoben et al., 2022b)
# 
# Knoben, W. J. M., Clark, M. P., Bales, J., Bennett, A., Gharari, S., Marsh, C. B., 
#  Nijssen, B., Pietroniro, A., Spiteri, R. J., Tang, G., Tarboton, D. G., & Wood, A. W. 
#  (2022a). Community Workflows to Advance Reproducibility in Hydrologic Modeling: 
#  Separating model-agnostic and model-specific configuration steps in applications 
#  of large-domain hydrologic models [Preprint]. Hydrology. 
#  https://doi.org/10.1002/essoar.10509195.2
#    
# Knoben, W. J. M., Marsh, C. B., & Tang, G. (2022b). CH-Earth/CWARHM: Initial 
#   release (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.6968609

import numpy as np
import os
from osgeo import gdal
import sys
from pathlib import Path
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

def read_merit_credentials(file_path = 'default'):
    
    '''Reads MERIT Hydro login credentials from [file_path].
       Unless specified, file is assumed to be ~/.merit'''
       
    if file_path == 'default':
        file_path = os.path.expanduser("~/.merit")
      
    merit_login = {}
    with open(file_path) as file:
        for line in file:
            (key, val) = line.split(':')
            merit_login[key] = val.strip() # remove whitespace, newlines
    
    # Get the authentication details
    usr = merit_login['name']
    pwd = merit_login['pass']
    
    return usr,pwd


def read_hydroshare_credentials(file_path = 'default'):
    
    '''Reads HydroShare login credentials from [file_path].
       Unless specified, file is assumed to be ~/.hydroshare'''

    if file_path == 'default':
        file_path = os.path.expanduser("~/.hydroshare")
      
    hs_login = {}
    with open(file_path) as file:
        for line in file:
            (key, val) = line.split(':')
            hs_login[key] = val.strip() # remove whitespace, newlines
    
    # Get the authentication details
    usr = hs_login['name']
    pwd = hs_login['pass']
    
    return usr,pwd


def convert_coordinates_to_merit_download_lists(coords):
    
    '''Converts [coords] as (lon_min,lon_max,lat_min,lat_max) to lists that 
       can be used to download various MERIT Hydro files for that area.'''
    
    # Convert area string into list
    coords = coords.split(',')

    # Store coordinates as floats in individual variables
    domain_min_lon = np.array(float(coords[0]))
    domain_max_lon = np.array(float(coords[1]))
    domain_min_lat = np.array(float(coords[2]))
    domain_max_lat = np.array(float(coords[3]))
    
    # Define the edges of the download areas
    lon_right_edge  = np.array([-150,-120, -90,-60,-30,  0,30,60,90,120,150,180])
    lon_left_edge   = np.array([-180,-150,-120,-90,-60,-30, 0,30,60, 90,120,150])
    lat_bottom_edge = np.array([-60,-30,0, 30,60]) # NOTE: latitudes -90 to -60 are NOT part of the MERIT domain
    lat_top_edge    = np.array([-30,  0,30,60,90]) 
    
    # Define the download variables
    dl_lon_all = np.array(['w180','w150','w120','w090','w060','w030','e000','e030','e060','e090','e120','e150'])
    dl_lat_all = np.array(['s60','s30','n00','n30','n60'])
   
    # Find the lower-left corners of each download square
    dl_lons = dl_lon_all[(domain_min_lon < lon_right_edge) & (domain_max_lon >= lon_left_edge)]
    dl_lats = dl_lat_all[(domain_min_lat < lat_top_edge) & (domain_max_lat >= lat_bottom_edge)]

    return dl_lons,dl_lats


def download_merit_hydro_grid(url,usr,pwd,dest_folder,
                              retries_max=10):
    
    '''Downloads MERIT Hydro [URL] into [dest_folder],
       using [usr] and [pwd] as login credentials.'''
    
    # Extract the filename from the URL
    file_name = url.split('/')[-1].strip() # Get the last part of the url, strip whitespace and characters
    
    # Check if file already exists in destination
    if os.path.isfile(dest_folder / file_name):
        print('WARNING: download_merit_hydro_grid: file {} already exists. Aborting download.'.format(dest_folder/file_name))
        return
        
    # Check if there is data for this specific location
    if ('n00' in file_name and 'w150' in file_name) or \
       ('s60' in file_name and 'w150' in file_name) or \
       ('s60' in file_name and 'w120' in file_name):
        print('WARNING: download_merit_hydro_grid: MERIT Hydro data does not contain data for {}. Aborting download.'.format(file_name))
        return
        
    # Make the keywords for request.get() and download the file
    kwargs={'auth': (usr,pwd)}
    cs.download_url_into_folder(url,dest_folder,requests_kwargs=kwargs)
    
    return


def merge_merit_downloads_into_area_of_interest(input_files,output_file,subset_window,no_data_value=''):
    
    '''Merges individual downloaded Merit Hydro GeoTIFF files into a single GeoTIFF file.
       Performs a subsetting of the individual files to only include area of interest in single file.
       Assumption: subsetting window is covered by the individual files - no checks to confirm this'''
    
    # Create a virtual dataset (VRT) of all individual GeoTIFF files that's not written to disk (argument '')
    vrt_options = gdal.BuildVRTOptions(resolution='highest')
    vrt = gdal.BuildVRT('', input_files, options=vrt_options)
    
    # Convert the VRT into a GeoTIFF file on disk
    tif_options = gdal.TranslateOptions(format='GTiff', 
                                        projWin=subset_window, 
                                        creationOptions=['COMPRESS=DEFLATE','BIGTIFF=YES']) 
    if no_data_value:
        tif_options = gdal.TranslateOptions(format='GTiff', 
                                            projWin=subset_window,
                                            noData=no_data_value,
                                            creationOptions=['COMPRESS=DEFLATE','BIGTIFF=YES']) 
    
    tif = gdal.Translate(output_file, vrt, options=tif_options)
    
    # Ensure merged GeoTIFF is actually written to disk
    # See: https://gdal.org/api/python_gotchas.html#saving-and-closing-datasets-datasources
    vrt = None
    tif = None
    
    return