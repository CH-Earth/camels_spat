'''Functions for processing of geospatial datasets'''

from bs4 import BeautifulSoup
import glob
import numpy as np
from osgeo import gdal
import requests
from urllib.parse import urljoin

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
        
    else:
        print(f'WARNING: geospatial_coordinates_to_download_coordinates(): no code found to process {product}. Returning input as output.')
        out = coords

    print(f'Returning coordinates as type {type(out)} for use with {product} download code.')
    return out

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