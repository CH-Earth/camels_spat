'''Functions for processing of geospatial datasets'''

from bs4 import BeautifulSoup
import geopandas as gpd
import glob
import numpy as np
import os
from osgeo import gdal
import rasterio
import re
import requests
from shapely.geometry import box
import shutil
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
    
    elif product.lower() == 'lgrip30':
        
        # Download edge values
        lon_left_edge   = np.arange(-180,180,10) # = array([-180,-170,..,160,170])
        lat_bottom_edge = np.arange(-60,70,10) # NOTE: latitudes (-90 to -60) and > 60N are NOT part of the LGRIP30 domain
        
        # Indices if closest lowest
        lon_min_i = np.where(lon_left_edge <= domain_min_lon)[0]
        lon_max_i = np.where(lon_left_edge <= domain_max_lon)[0]
        lat_min_i = np.where(lat_bottom_edge <= domain_min_lat)[0]
        lat_max_i = np.where(lat_bottom_edge <= domain_max_lat)[0]
        
        # Convert to coordinate output (string)
        out = f'{lon_left_edge[lon_min_i[-1]]},{lon_left_edge[lon_max_i[-1]]},{lat_bottom_edge[lat_min_i[-1]]},{lat_bottom_edge[lat_max_i[-1]]}'
    
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
def find_folders_on_webpage(url,product='soilgrids'):
    
    # Send an HTTP GET request to the URL and get the HTML content.
    response = requests.get(url)
    
    # Check if the request was successful (status code 200).
    if response.status_code == 200:
        # Parse the HTML content using BeautifulSoup.
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Find and extract folder links. You'll need to inspect the HTML structure
        # of the webpage to determine the appropriate HTML tags and attributes.
        # For example, if the folder links are within <a> tags with a specific class:
        if product.lower() == 'soilgrids':
            folder_links = [a['href'] for a in soup.find_all('a', href=True) if 'tile' in a['href']]
        elif product.lower() == 'mcd15a2h.061':
            pattern = r'\d{4}\.\d{2}\.\d{2}' # E.g., '2002.07.28'
            folder_links = [a['href'] for a in soup.find_all('a', href=True) if re.match(pattern, a['href'])]
    
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
        #files = [a['href'] for a in soup.find_all('a', href=True) if extension in a['href']]
        files = [a['href'] for a in soup.find_all('a', href=True) if a['href'].endswith(extension)]
        
        # Convert relative URLs to absolute URLs.
        file_links = [urljoin(url, file) for file in files]
        
    else:
        print("Failed to retrieve the webpage. Status code:", response.status_code)        
    
    return file_links

def find_file_urls_in_webpage_folders(folders, extension='.tif'):

    file_urls = []
    for folder_url in folders:
        file_urls.append(find_files_in_webpage_folder(folder_url, extension=extension))

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

# --- MODIS LAI
def hdf_tile_in_north_america(url):

    '''Checks if a HDF URL falls within the North America domain (sort of).
       See Fig. 2 here: https://lpdaac.usgs.gov/documents/926/MOD15_User_Guide_V61.pdf'''

    # Separate parts
    file_name = os.path.basename(url)
    file_parts = file_name.split('.')

    # Check if this matches expectations
    pattern = r'h\d{2}v\d{2}'
    if not re.match(pattern,file_parts[2]):
        print('WARNING: hdf_tile_in_north_america(): URL {url} does not match expected pattern')
        return False

    # Check h and v values against what we want
    h = int(file_parts[2][1:3])
    v = int(file_parts[2][4:6])

    flag = True
    if h > 15: flag = False
    if v > 8:  flag = False
    if h < 7:  flag = False

    return flag

def download_modis_into_day_folder(download_folder, url):

    # Make the names and paths
    folder_name = os.path.basename(os.path.dirname(url))
    folder = download_folder/folder_name
    folder.mkdir(parents=True, exist_ok=True)
    
    # Check if this file is in the domain we want
    if hdf_tile_in_north_america(url):
        cs.download_url_into_folder(url,folder)
    
    return folder

def process_daily_modis_hdf_to_tif(in_folder, out_folder,
                                   subdataset_front='HDF4_EOS:EOS_GRID',
                                   subdataset_back='MOD_Grid_MOD15A2H:Lai_500m',
                                   to_CRS='EPSG:4326',
                                   subset_window=[]):

    # Find files
    files = cs.find_files_in_folder(in_folder)

    # Transform to gdal.BuildVrt inputs, by specifying the subdataset we want
    hdf_inlist = [f'{subdataset_front}:"{file}":{subdataset_back}' for file in files if file.endswith('.hdf')]

    # Create a virtual dataset (VRT) of all individual GeoTIFF files that's not written to disk (argument '')
    vrt_options = gdal.BuildVRTOptions(resolution='highest')
    vrt = gdal.BuildVRT('', hdf_inlist, options=vrt_options)

    # Reproject the VRT
    if subset_window:
        warp_options = gdal.WarpOptions(format='VRT', 
                                        outputBounds=[subset_window[0],subset_window[3],subset_window[2],subset_window[1]], 
                                        dstSRS=to_CRS)
    else:
        warp_options = gdal.WarpOptions(format='VRT', dstSRS=to_CRS)
    vrt_4326 = gdal.Warp('', vrt, options=warp_options)

    # Subset and convert to geotiff
    if subset_window:
        tif_options = gdal.TranslateOptions(format='GTiff', projWin=subset_window, creationOptions=['COMPRESS=DEFLATE','BIGTIFF=YES'])
    else:
        tif_options = gdal.TranslateOptions(format='GTiff', creationOptions=['COMPRESS=DEFLATE','BIGTIFF=YES'])
    
    output_name = os.path.basename(in_folder).replace('.','') + '_' + subdataset_back.replace(':','_') + '.tif'
    gdal.Translate(str(out_folder/output_name), vrt_4326, options=tif_options).FlushCache()

    # Flush cache
    vrt = None
    vrt_4326 = None

def check_modis_interrupt_status(url,raw_folder):

    '''Checks if we already have a merged file for the URL we're about to download'''

    # Extract the date from URL, in the format (yyyymmdd) we use it in raw_folder
    file_date = os.path.basename(os.path.dirname(url)).replace('.','')

    # Find what we already have
    complete_files = cs.find_files_in_folder(raw_folder)

    # Check if we have already completed processing for the day this specific URL is part of
    is_complete = False # Start assuming that we did not complete the downloads for this day
    for complete in complete_files:
        if file_date in complete:
            is_complete = True

    return is_complete

# --- LGRIP30 agriculture
def convert_coordinates_to_lgrip_download_lists(coords):
    
    '''Converts [coords] as (lon_min,lon_max,lat_min,lat_max) to lists that 
       can be used to download various LGRIP30 files for that area.'''
    
    import numpy as np

    # Convert area string into list
    coords = coords.split(',')

    # Store coordinates as floats in individual variables
    domain_min_lon = np.array(float(coords[0]))
    domain_max_lon = np.array(float(coords[1]))
    domain_min_lat = np.array(float(coords[2]))
    domain_max_lat = np.array(float(coords[3]))
    
    # Define the edges of the download areas
    lon_right_edge  = np.arange(-170,190,10)
    lon_left_edge   = np.arange(-180,180,10)
    lat_bottom_edge = np.arange(-40,90,10) 
    lat_top_edge    = np.arange(-30,100,10) 
    
    # Define the download variables
    lon_list = []
    for item in lon_left_edge:
        if item < 0:
            lon_list.append(f'W{np.abs(item):02d}')
        else:
            lon_list.append(f'E{item:02d}')
    lat_list = []
    for item in lat_bottom_edge:
        if item < 0:
            lat_list.append(f'S{np.abs(item):02d}')
        else:
            lat_list.append(f'N{item:02d}')
    
    dl_lon_all = np.array(lon_list)
    dl_lat_all = np.array(lat_list)
   
    # Find the lower-left corners of each download square
    dl_lons = dl_lon_all[(domain_min_lon < lon_right_edge) & (domain_max_lon >= lon_left_edge)]
    dl_lats = dl_lat_all[(domain_min_lat < lat_top_edge) & (domain_max_lat >= lat_bottom_edge)]

    return dl_lons,dl_lats

def replace_data_value_in_geotiff(infile,outfile,old,new,large_file=False):

    '''Replaces _old_ value in _infile_ with _new_ in _outfile_'''

    if not large_file:
        with rasterio.open(infile, 'r') as src:
            # Read the data as a numpy array
            data = src.read(1)
    
            # Replace all old values with new
            data[data == old] = new
    
            # Get metadata from the source file
            profile = src.profile
    
           # Write the modified data to a new GeoTIFF file
            with rasterio.open(outfile, 'w', **profile) as dst:
                dst.write(data, 1)
    else:

        # Copy the source file to the output file, so we have all the geospatial info already and we only need to overwrite the data values
        shutil.copy(infile,outfile)
        
        with rasterio.open(infile, 'r') as src,  rasterio.open(outfile, 'r+') as dst:
            # Loop through the windows of the GeoTIFF
            for ij,window in src.block_windows(1):
                print(ij)
                
                # Read the data as a numpy array
                data = src.read(1, window=window)

                # Replace all old values with new
                data[data == old] = new

                # Write the modified data back to the output GeoTIFF
                #with rasterio.open(outfile, 'r+', driver='GTiff', 
                #                   width=window.width, height=window.height, count=1, dtype=data.dtype, 
                #                   transform=src.window_transform(window)) as dst:
                dst.write(data, 1, window=window)

    return

def compress_geotiff(infile,outfile):

    # Open the input GeoTIFF dataset
    input_dataset = gdal.Open(infile, gdal.GA_ReadOnly)

    # Create a copy of the input dataset with compression
    driver = gdal.GetDriverByName("GTiff")
    output_dataset = driver.CreateCopy(outfile, input_dataset, options=['COMPRESS=DEFLATE','BIGTIFF=YES'])

    # Flush cache
    input_dataset = None
    output_dataset = None
    return

# --- GLHYMPS
def process_glhymps(gdb_file,shp_file,download_area):

    '''Subsets the global GLHYMPS GeoDataBase to a bounding box, then saves to ESRI shapefile'''

    # Opent he GDB and convert to EPSG:4326
    gdf = gpd.read_file(gdb_file)
    gdf = gdf.to_crs('EPSG:4326')

    # Define the bounding box
    coords = [float(coord) for coord in download_area.split(',')]
    bounding_box = box(coords[0], coords[2], coords[1], coords[3])

    # Create the subset
    subset = gdf[gdf.geometry.intersects(bounding_box)]
    subset.to_file(shp_file)

    return # nothing