'''Functions for processing of geospatial datasets'''

from bs4 import BeautifulSoup
import glob
import numpy as np
import os
from osgeo import gdal
import requests
from urllib.parse import urljoin

import sys
from pathlib import Path
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

# --- General
def geospatial_coordinates_to_download_coordinates(coords, product):

    '''Converts general download coodinates (lon_min, lon_max,lat_min,lat_max) to the data-specific ones'''

    # Store coordinates as floats in individual variables
    coords = coords.split(',')
    domain_min_lon = np.array(float(coords[0]))
    domain_max_lon = np.array(float(coords[1]))
    domain_min_lat = np.array(float(coords[2]))
    domain_max_lat = np.array(float(coords[3]))

    # Round, if necessary
    if product.lower() == 'merit':
        
        # Download edge values
        lon_left_edge   = np.array([-180,-150,-120,-90,-60,-30, 0,30,60, 90,120,150])
        lat_bottom_edge = np.array([-60,-30,0, 30,60]) # NOTE: latitudes -90 to -60 are NOT part of the MERIT domain

        # Indices if closest lowest
        lon_min_i = np.where(lon_left_edge <= domain_min_lon)[0]
        lon_max_i = np.where(lon_left_edge <= domain_max_lon)[0]
        lat_min_i = np.where(lat_bottom_edge <= domain_min_lat)[0]
        lat_max_i = np.where(lat_bottom_edge <= domain_max_lat)[0]

        # Convert to coordinate output (string)
        out = f'{lon_left_edge[lon_min_i[-1]]},{lon_left_edge[lon_max_i[-1]]},{lat_bottom_edge[lat_min_i[-1]]},{lat_bottom_edge[lat_max_i[-1]]}'

    elif product.lower() == 'soilgrids':

        # Return in format that's good to go for downloading (tuple)
        out = [coords[0], coords[3], coords[1], coords[2]]
    
    elif product.lower() == 'glclu2019':

        # Download edge values
        lon_left_edge   = np.arange(-180,180,10) # = array([-180,-170,..,160,170])
        lat_top_edge = np.arange(-40,90,10) # NOTE: latitudes (-90 to -50) and > 80N are NOT part of the GLCLU2019 domain

        # Indices if closest lowest
        lon_min_i = np.where(lon_left_edge <= domain_min_lon)[0]
        lon_max_i = np.where(lon_left_edge <= domain_max_lon)[0]
        lat_min_i = np.where(lat_top_edge <= domain_min_lat)[0]
        lat_max_i = np.where(lat_top_edge <= domain_max_lat)[0]

        # Convert to coordinate output (string)
        out = f'{lon_left_edge[lon_min_i[-1]]},{lon_left_edge[lon_max_i[-1]]},{lat_top_edge[lat_min_i[-1]]},{lat_top_edge[lat_max_i[-1]]}'
    
    else:
        print(f'WARNING: geospatial_coordinates_to_download_coordinates(): no code found to process {product}. Returning input as output.')
        out = coords

    print(f'Returning coordinates as type {type(out)} for use with {product} download code.')
    return out

def subset_tif(infile,outfile,subset_window):

    '''Subsets a GeoTIFF file'''

    infile = str(infile)
    outfile = str(outfile)
    
    tif_options = gdal.TranslateOptions(format='GTiff', projWin=subset_window, creationOptions=['COMPRESS=DEFLATE','BIGTIFF=YES']) 
    tif = gdal.Translate(outfile, infile, options=tif_options)
    tif = None
    
    return

# --- Soilgrids
def find_folders_on_webpage(url):
    
    # Send an HTTP GET request to the URL and get the HTML content.
    response = requests.get(url)
    
    # Check if the request was successful (status code 200).
    if response.status_code == 200:
        # Parse the HTML content using BeautifulSoup.
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Find and extract folder links. You'll need to inspect the HTML structure
        # of the webpage to determine the appropriate HTML tags and attributes.
        # For example, if the folder links are within <a> tags with a specific class:
        folder_links = [a['href'] for a in soup.find_all('a', href=True) if 'tile' in a['href']]
    
        # Convert relative URLs to absolute URLs.
        folder_links = [urljoin(url, link) for link in folder_links]
        
    else:
        print("Failed to retrieve the webpage. Status code:", response.status_code)
    
    return folder_links

def find_files_in_webpage_folder(url, extension='.tif'):
    
    # Send an HTTP GET request to the URL and get the HTML content.
    response = requests.get(url)
    
    # Check if the request was successful (status code 200).
    if response.status_code == 200:
        # Parse the HTML content using BeautifulSoup.
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Grab any hrefs that have the correct extension
        files = [a['href'] for a in soup.find_all('a', href=True) if extension in a['href']]
        
        # Convert relative URLs to absolute URLs.
        file_links = [urljoin(url, file) for file in files]
        
    else:
        print("Failed to retrieve the webpage. Status code:", response.status_code)        
    
    return file_links

def find_file_urls_in_webpage_folders(folders):

    file_urls = []
    for folder_url in folders:
        file_urls.append(find_files_in_webpage_folder(folder_url))

    # Ensure we return one list with values, instead of a list of lists
    flat_urls = [item for sublist in file_urls for item in sublist]

    return flat_urls

def find_files_in_folder(folder,extension=''):

    '''Searches folder for files, optionally only those with specified extension'''

    files = glob.glob(f'{folder}/*{extension}')

    return files

def process_soilgrids_tiles_into_single_geotiff(in_folder,out_folder,out_file,
                                                to_crs='',subset_window=''):

    # Handle inputs
    out = Path(out_folder)
    out.mkdir(parents=True, exist_ok=True)
    des = str(out/out_file)
    
    # Find the tiles
    files = find_files_in_folder(in_folder,extension='.tif')
    
    # Create a virtual dataset (VRT) of all individual GeoTIFF files that's not written to disk (argument '')
    vrt_options = gdal.BuildVRTOptions(resolution='highest')
    vrt = gdal.BuildVRT('', files, options=vrt_options)
        
    # Reproject into EPSG:4326
    if to_crs:
        # Convert the vrt to a different coordinate system
        driver = gdal.GetDriverByName('VRT')
        vrt_old = driver.CreateCopy('', vrt)
        vrt = None
        vrt = gdal.Warp('',vrt_old, format='vrt', dstSRS=to_crs)

    # Subset the area we want
    if subset_window:
        tif_options = gdal.TranslateOptions(format='GTiff', projWin=subset_window, creationOptions=['COMPRESS=DEFLATE','BIGTIFF=YES'])
    else:
        tif_options = gdal.TranslateOptions(format='GTiff', creationOptions=['COMPRESS=DEFLATE','BIGTIFF=YES'])

    # Save to GeoTIFF
    tif = gdal.Translate(des, vrt, options=tif_options)

    # Ensure merged GeoTIFF is actually written to disk
    # See: https://gdal.org/api/python_gotchas.html#saving-and-closing-datasets-datasources
    tif = None
    vrt = None
    
    return # nothing

def download_all_soilgrids_tiles_into_folder(url, dest):

    '''Follows a SOILGRIDS main webpage link and downloads all geotiffs available in sub-folders on the link'''
    
    # Determine which folders are present on the main download page - these each contain individual GeoTIFF files
    online_folders = find_folders_on_webpage(url)

    # Determine which files are present in the subfolders
    online_files = find_file_urls_in_webpage_folders(online_folders)

    # Download the files
    dest.mkdir(parents=True, exist_ok=True) # Ensure the destination exists
    for file in online_files:
        cs.download_url_into_folder(file,dest)
    
    return

# --- Global Land Class Land Use 2019 data
def convert_coordinates_to_glclu2019_download_lists(coords):
    
    '''Converts [coords] as (lon_min,lon_max,lat_min,lat_max) to lists that 
       can be used to download various Global Land Cover Land Use 2019 files for that area.'''

    # Convert area string into list
    coords = coords.split(',')

    # Store coordinates as floats in individual variables
    domain_min_lon = np.array(float(coords[0]))
    domain_max_lon = np.array(float(coords[1]))
    domain_min_lat = np.array(float(coords[2]))
    domain_max_lat = np.array(float(coords[3]))
    
    # Define the edges of the download areas
    lon_right_edge  = np.arange(-170,190,10) # = array([-170,-160,..,170,180])
    lon_left_edge   = np.arange(-180,180,10) # = array([-180,-170,..,160,170])
    lat_bottom_edge = np.arange(-50,80,10) # lat < 50S and > 80N are not part of the domain
    lat_top_edge    = np.arange(-40,90,10) 
    
    # Define the download variables
    lon_list = []
    for item in lon_left_edge:
        if item < 0:
            lon_list.append(f'{np.abs(item):03d}W')
        else:
            lon_list.append(f'{item:03d}E')
    lat_list = []
    for item in lat_top_edge:
        if item < 0:
            lat_list.append(f'{np.abs(item):02d}S')
        else:
            lat_list.append(f'{item:02d}N')
    
    dl_lon_all = np.array(lon_list)
    dl_lat_all = np.array(lat_list)
   
    # Find the upper-left corners of each download square
    dl_lons = dl_lon_all[(domain_min_lon < lon_right_edge) & (domain_max_lon >= lon_left_edge)]
    dl_lats = dl_lat_all[(domain_min_lat < lat_top_edge) & (domain_max_lat >= lat_bottom_edge)]

    return dl_lons,dl_lats

def download_glclu2019_grid(url, dest_folder, retries_max=10):
    
    # Extract the filename from the URL
    file_name = url.split('/')[-1].strip() # Get the last part of the url, strip whitespace and characters
    
    # Check if file already exists in destination
    if os.path.isfile(dest_folder / file_name):
        print('WARNING: download_glclu2019_grid: file {} already exists. Aborting download.'.format(dest_folder/file_name))
        return
        
    # Check if there is data for this specific location
    if not url_file_exists(url):
        print('WARNING: download_glclu2019_grid: Global Land Cover Land Use data does not contain data for {}. Aborting download.'.format(file_name))
        return
        
    # Download the file
    cs.download_url_into_folder(url,dest_folder)
    
    return

def url_file_exists(url):
    response = requests.head(url)
    return response.status_code == 200

def merge_glclu2019_files_into_one(merged_file, src_folder, des_folder, download_area):

    # Find the file names
    all_files = []
    for dir_path, dir_names, file_names in os.walk(src_folder):
        for file_name in file_names:
            if file_name.endswith('.tif'): # ensure we don't accidentally get .aux files from QGIS or something similar
                all_files.append(os.path.join(dir_path,file_name))

    # Ensure destination exists
    des_folder.mkdir(parents=True, exist_ok=True)

    # Convert subsetting area into a usable GDAL setting
    # subset_area = [lon_min, lon_max, lat_min, lat_max]
    # GDAL window = [ulx, uly, lrx, lry]; [upper left x, upper left y, lower right x, lower right y]
    # Mapping:
    #   ulx = lon_min = subset_area[0]
    #   uly = lat_max = subset_area[3]
    #   lrx = lon_max = subset_area[1]
    #   lry = lat_min = subset_area[2]
    subset_coor = download_area.split(',')
    window = [subset_coor[0], subset_coor[3], subset_coor[1], subset_coor[2]]

    # Merge into area of interest
    cs.merge_merit_downloads_into_area_of_interest(all_files, str(des_folder/merged_file), window) # originally built for MERIT, should work here too
    
    return

