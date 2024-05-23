'''Functions for processing of geospatial datasets (forcing and parameters)'''

from bs4 import BeautifulSoup
import cdsapi
from datetime import datetime, timedelta
import geopandas as gpd
import glob
import netCDF4 as nc4
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import os
from osgeo import gdal
gdal.UseExceptions()
from pathlib import Path
import rasterio
import re
import requests
from shapely import make_valid
from shapely.geometry import box
import shutil
import sys
import time
from urllib.parse import urljoin
import warnings
import xarray as xr

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

def subset_geotiff_to_shapefile(src_file,src_shape,des_folder,
                                buffer=False,
                                out_no_data = None):

    # Input cleaning
    des_folder.mkdir(parents=True, exist_ok=True)
    des_file  = str(des_folder / os.path.basename(src_file))
    src_file  = str(src_file)
    src_shape = str(src_shape)

    # Handle the specific MERIT case - buffering breaks things here for no apparent reason
    # We can get away with no buffering because the shapes are derived from the MERIT DEM and
    #  thus align perfectly
    if 'merit' in src_file:
        buffer = False
    
    # Handle buffering of shapefile, if requested
    if buffer:
        tmp_shape = src_shape.replace('.shp','_TEMP.shp')
    
        # Find buffer distance
        src_tiff = gdal.Open(src_file, gdal.GA_ReadOnly)
        pixel_x  = src_tiff.GetGeoTransform()[1]
        pixel_y  = src_tiff.GetGeoTransform()[5]
        buffer   = 0.5*(pixel_x**2 + pixel_y**2)**(0.5) # I.e., half the maximum distance from center to edge of pixel
        src_tiff = None
        
        # Temporarily block warnings: 
        # gpd will tell us that buffering in EPSG:4326 is not accurate - this is fine because we're
        # buffering in lat/lon units
        with warnings.catch_warnings():
            warnings.simplefilter('ignore') 
    
            # Buffer the shapefile
            shp = gpd.read_file(src_shape)
            shp['geometry'] = shp.buffer(buffer)
            shp.geometry = shp.apply(lambda row: make_valid(row.geometry), axis=1)
            shp.to_file(tmp_shape)
    else:
        # Not using buffered shape, but code below still needs 'tmp_shape' to have a value
        tmp_shape = src_shape

    # Clip
    gdal.Warp(destNameOrDestDS = des_file,
              srcDSOrSrcDSTab  = src_file,
              cutlineDSName    = tmp_shape, # vector file
              cropToCutline    = True, # Select True
              copyMetadata     = True, # optional
              #dstAlpha         = True, # Dropping the alpha band saves half the file size
              dstNodata        = out_no_data,
              srcSRS           = 'EPSG:4326',
              dstSRS           = 'EPSG:4326',
              #resampleAlg      = "nearestneighbour"
             )
    
    # Remove buffered shapefile
    if buffer:
        remove_these = glob.glob(tmp_shape.replace('.shp','.*'))
        for file in remove_these:
            os.remove(file)

def subset_shapefile_to_shapefile(src,src_file,src_shape,des_folder):

    # Input cleaning
    des_folder.mkdir(parents=True, exist_ok=True)
    des_file  = str(des_folder / os.path.basename(src_file))
    src_file  = str(src_file)
    src_shape = str(src_shape)
   
    # Open the basin shapefile
    shp = gpd.read_file(src_shape)
    
    # Loop over the geometries to check if they intersect, then use this info to create a HydroLAKES subset
    # Note: this is cleaner than des = gpd.overlay(src,shp, how='intersection'), because this alternative
    #   approach may be faster, but it clips the lake polygons to the catchment extent. The consequence of
    #   this is that the lake area reported as part of the HydroLAKES subset may be inaccurate, in cases
    #   where that polygon was clipped by the catchment outline. With the current approach we return the
    #   complete lake polygons.
    if 'hydrolakes' in src_file.lower():
        src['mask'] = src.apply(lambda row: row.geometry.intersects(shp.geometry), axis=1)
        des = src[src['mask'] == True].copy().reset_index(drop=True)
        des = des.drop('mask', axis=1)
    elif 'glhymps' in src_file.lower():
        des = gpd.overlay(src,shp,how='intersection') # Faster than above, clipping these polygons is fine because they don't contain an 'area' field
    
    # To file
    if len(des) > 0:
        des.to_file(des_file) # write the shapefile to file
    else:
        txt_file = f'{os.path.splitext(des_file)[0]}.txt'
        with open(txt_file,'w') as f:
            f.write('Source data does not intersect with catchment shapefile')

    return # nothing

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
        elif product.lower() == 'mcd12q1.061':
            pattern = r'\d{4}\.\d{2}\.\d{2}' # E.g., '2002.07.28'
            folder_links = [a['href'] for a in soup.find_all('a', href=True) if re.match(pattern, a['href'])]
        elif product.lower() == 'rdrs':
            pattern = r'\d{4}' # E.g. '1980'
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

def check_geotiff_dtypes_and_convert(files):
    
    # Mapping from string to GDAL datatype constant
    datatype_mapping = {
        'Byte': gdal.GDT_Byte,
        'UInt16': gdal.GDT_UInt16,
        'Int16': gdal.GDT_Int16,
        'UInt32': gdal.GDT_UInt32,
        'Int32': gdal.GDT_Int32,
        'Float32': gdal.GDT_Float32,
        'Float64': gdal.GDT_Float64
    }

    # Check the data types in the files
    data_types = []
    for file in files:
        ds = gdal.Open(file)
        data_types.append(gdal.GetDataTypeName(ds.GetRasterBand(1).DataType))
        ds = None

    # Count occurrences in the list
    dtype_counts = {}
    for dtype in data_types:
        if dtype in dtype_counts:
            dtype_counts[dtype] += 1
        else:
            dtype_counts[dtype] = 1
    most_common_dtype = max(dtype_counts, key=dtype_counts.get)
    unique_dtypes = list(dtype_counts.keys())

    # Check if data types are consistent
    if len(unique_dtypes) != 1:
        print(f'--- WARNING: process_soilgrids_tiles_into_single_geotiff(): multiple data types found in input files: {unique_dtypes}.')
        print(f'gdalbuildvrt will not process these. Converting all files to majority type {most_common_dtype}.')
        
        # Convert data types
        desired_dtype = datatype_mapping[most_common_dtype]
        for file in files:
            
            ds = gdal.Open(file)
            current_dtype = gdal.GetDataTypeName(ds.GetRasterBand(1).DataType)
            if current_dtype != most_common_dtype:
                print(f'processing {file}')

                # Close the current file
                print(f'closing {file}')
                ds = None

                # Copy the file to a temporary location
                temp_file = file.replace('.tif','_TEMP.tif')
                print(f'Copying {temp_file}')
                shutil.copy(file,temp_file)

                # Use gdal.Translate to convert datatype
                print(f'Converting {temp_file}')
                ds = gdal.Open(temp_file)
                gdal.Translate(file, ds, format='GTiff', outputType= desired_dtype)
                ds = None

                # Remove the temporary file
                print(f'Removing {temp_file}')
                os.remove(temp_file)
            else:
                ds = None

    return # nothing

def select_uint16_only(files):

    # Check the data types in the files
    uint16_files = []
    for file in files:
        ds = gdal.Open(file)
        dtype = gdal.GetDataTypeName(ds.GetRasterBand(1).DataType)
        if dtype == 'UInt16':
            uint16_files.append(file)
        ds = None

    return uint16_files

def process_soilgrids_tiles_into_single_geotiff(in_folder,out_folder,out_file,
                                                to_crs='',subset_window=''):

    # Handle inputs
    out = Path(out_folder)
    out.mkdir(parents=True, exist_ok=True)
    des = str(out/out_file)
    
    # Find the tiles
    files = find_files_in_folder(in_folder,extension='.tif')
    
    # Ensure data types in files match
    #check_geotiff_dtypes_and_convert(files) # Note: even after converting data type we run into issues about mismatching projection between tiles 
    #files = select_uint16_only(files)
    
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

# --- Forcing folder setup ---
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
    
    if target == 'EM-Earth':
        
        # Round to EM-Earth 0.10 degree resolution
        rounded_lon = [math.floor(lon[0]*20)/20, math.ceil(lon[1]*20)/20]
        rounded_lat = [math.floor(lat[0]*20)/20, math.ceil(lat[1]*20)/20]

        # Find if we are still in the representative area of a different ERA5 grid cell
        if lat[0] > rounded_lat[0]+0.05:
            rounded_lat[0] += 0.10
        if lon[0] > rounded_lon[0]+0.05:
            rounded_lon[0] += 0.10
        if lat[1] < rounded_lat[1]-0.05:
            rounded_lat[1] -= 0.10
        if lon[1] < rounded_lon[1]-0.05:
            rounded_lon[1] -= 0.10
    
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

def find_first_day_next_month(time):

    '''Takes a datetime and find day 1 on the next month'''

    # Initial settings: day 1 of a month in the same year
    new_day = 1
    new_month = time.month+1
    new_year = time.year
    
    # Check if we're at the end of the year
    if new_month == 13: 
        new_month = 1
        new_year = time.year+1

    return time.replace(year=new_year, month=new_month, day=new_day)

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

        # Find what the 1st day of the next month is
        new_start = find_first_day_next_month(cur)

        # Find the end date of the current month
        cur_end = new_start - timedelta(days=1)

        # Ensure this does not step over our end date
        if cur_end >= end:
            cur_end = end

        # Add to end list
        end_l.append(cur_end)

        # Update variable for while loop
        cur = new_start

    return start_l,end_l

# --- ERA5 processing ---
def extract_ERA5_subset(infile, outfile, coords):
    
    '''Subsets an existing ERA5 forcing file by coordinates "latmax / lonmin / latmin / lonmax"'''

    # Source: https://github.com/CH-Earth/CWARHM/blob/main/0_tools/ERA5_subset_forcing_file_by_lat_lon.py

    # Notes:
    # 1. Works for North American continent
    # 2. Works for two specific ERA5 file layouts

    # Assumptions
    # 1. Latitude = [-90,90]
    # 2. Longitude = [-180,180]
    
    # Split coordinates
    coords = coords.split('/') # split string
    coords = [float(value) for value in coords] # string to array
    latmax = coords[0]
    lonmin = coords[1]
    latmin = coords[2]
    lonmax = coords[3]
    
    with xr.open_dataset(infile) as ds:

        # Handle specific cases
        if (ds['longitude'] > 180).any(): # convert ds longitude form 0/360 to -180/180
            lon = ds['longitude'].values
            lon[lon > 180] = lon[lon > 180] - 360
            ds['longitude'] = lon
        
        if (ds['latitude'] > 90).any(): # convert ds longitude form 0/180 to -90/90
            lat = ds['latitude'].values
            lat[lat > 90] = lat[lat > 90] - 180
            ds['latitue'] = lat

        # Subset
        ds_sub = ds.sel(latitude = slice(latmax, latmin), longitude = slice(lonmin, lonmax))
        ds_sub.to_netcdf(outfile)
        ds_sub.close()

def compare_forcing_data_and_shape_extents(save_path, nc_file, shp_file, nc_var, nc_time=0):

    '''Plots forcing grid and shapefile on top of each other'''

    catchment = gpd.read_file(shp_file)
    with xr.open_dataset(nc_file) as forcing:
        ax = plt.subplot()   
        if (len(forcing['latitude']) == 1) or (len(forcing['longitude']) == 1):
            # Case 1: in at least one direction, ERA5 file covers only a single grid cell. Needs special plotting attention
            lon = forcing['longitude'].values
            lat = forcing['latitude'].values
            grid = Rectangle( (lon.min()-0.125,lat.min()-0.125), lon.max()-lon.min()+0.25, lat.max()-lat.min()+0.25 )
            ax.add_patch(grid)
        else:
            # Case 2: more then a single ERA5 cell in either direction
            forcing[nc_var].isel(time=nc_time).plot(ax=ax,cbar_kwargs={'shrink': 0.95})
        catchment.plot(color='None', edgecolor='k', ax=ax);
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()

def compare_era5_netcdf_dimensions(surf_file, pres_file):

    '''Opens two netCDF files containing ERA5 data and checks if contents of latitude, longitude and time match'''

    # Extract space and time info for checks
    with nc4.Dataset(pres_file) as src1, nc4.Dataset(surf_file) as src2: # implicitly closes files
        pres_lat = src1.variables['latitude'][:]
        pres_lon = src1.variables['longitude'][:]
        pres_time = src1.variables['time'][:]
        surf_lat = src2.variables['latitude'][:]
        surf_lon = src2.variables['longitude'][:]
        surf_time = src2.variables['time'][:]

     # Update the pressure level coordinates
    pres_lat[pres_lat > 90] = pres_lat[pres_lat > 90] - 180
    pres_lon[pres_lon > 180] = pres_lon[pres_lon > 180] - 360

    # Perform the check
    if not all( [all(pres_lat == surf_lat), all(pres_lon == surf_lon), all(pres_time == surf_time)] ):
        err_txt = (f'ERROR: Dimension mismatch while merging {pressure_file} and {surface_file}. '
                    'Check latitude, longitude and time dimensions in both files. Continuing with next files.')
        return False,err_txt
    else:
        return True,'Variables match'

def merge_era5_surface_and_pressure_files(surf_file, pres_file, dest_file):

    '''Merges two ERA5 data files into dest_file'''

    # Transfer definitions
    dimensions_surf_transfer = ['longitude','latitude','time']
    variables_surf_transfer = ['msdwlwrf','msnlwrf','msdwswrf','msnswrf','mtpr','sp','mper']
    variables_pres_transfer = ['t','q','u','v']
    attr_names_expected = ['scale_factor','add_offset','_FillValue','missing_value','units','long_name','standard_name'] # these are the attributes we think each .nc variable has             
    loop_attr_copy_these = ['units','long_name','standard_name'] # we will define new values for _FillValue and missing_value when writing the .nc variables' attributes

    with nc4.Dataset(pres_file) as src1, nc4.Dataset(surf_file) as src2, nc4.Dataset(dest_file, 'w') as dest: # implicitly closes files

        # Set general attributes
        dest.setncattr('History','Created ' + time.ctime(time.time()))
        dest.setncattr('Reason','Merging separate ERA5 download files')
        dest.setncattr('Source','github.com/cH-Earth/camels_spat')

        # Copy meta data from the two source files
        if src1.getncattr('Conventions') == 'CF-1.6' and src2.getncattr('Conventions') == 'CF-1.6':
            dest.setncattr('Conventions','CF-1.6')
        dest.setncattr('History_source_file_1',src1.getncattr('history'))
        dest.setncattr('History_source_file_2',src2.getncattr('history'))

        # --- Dimensions: latitude, longitude, time
        # NOTE: we can use the lat/lon from the surface file (src2), because those are already in proper units. 
        # If there is a mismatch between surface and pressure we shouldn't have reached this point at all due to the earlier check
        for name, dimension in src2.dimensions.items():
            if dimension.isunlimited():
                dest.createDimension( name, None)
            else:
                dest.createDimension( name, len(dimension))

        # --- Get the surface level dimension variables (lat, lon, time)
        for name, variable in src2.variables.items():
    
            # Transfer lat, long and time variables because these don't have scaling factors
            if name in dimensions_surf_transfer:
                dest.createVariable(name, variable.datatype, variable.dimensions, fill_value = -999)
                dest[name].setncatts(src1[name].__dict__)
                dest.variables[name][:] = src2.variables[name][:]

        # ---  Transfer the surface level data first, for no particular reason
        # This should contain surface pressure (sp), downward longwave (msdwlwrf), downward shortwave (msdwswrf) and precipitation (mtpr)
        for name, variable in src2.variables.items():

            # Check that we are only using the names we expect, and thus the names for which we have the required code ready
            if name in variables_surf_transfer:
        
                # 0. Reset the dictionary that we keep attribute values in
                loop_attr_source_values = {name: 'n/a' for name in attr_names_expected}
        
                # 1a. Get the values of this variable from the source (this automatically applies scaling and offset)
                loop_val = variable[:]
        
                # 1b. Get the attributes for this variable from source
                for attrname in variable.ncattrs():
                    loop_attr_source_values[attrname] = variable.getncattr(attrname)
               
                # 2. Create the .nc variable 
                # Inputs: variable name as needed by SUMMA; data type: 'float'; dimensions; no need for fill value, because thevariable gets populated in this same script
                dest.createVariable(name, 'f4', ('time','latitude','longitude'), fill_value = False, zlib=True, shuffle=True)
        
                # 3a. Select the attributes we want to copy for this variable, based on the dictionary defined before the loop starts
                loop_attr_copy_values = {use_this: loop_attr_source_values[use_this] for use_this in loop_attr_copy_these}
        
                # 3b. Copy the attributes FIRST, so we don't run into any scaling/offset issues
                dest[name].setncattr('missing_value',-999)
                dest[name].setncatts(loop_attr_copy_values)
        
                # 3c. Copy the data SECOND
                dest[name][:] = loop_val

        # --- Transfer the pressure level variables next, using the same procedure as above
        for name, variable in src1.variables.items():
            if name in variables_pres_transfer:
        
                # 0. Reset the dictionary that we keep attribute values in
                loop_attr_source_values = {name: 'n/a' for name in attr_names_expected}
        
                # 1a. Get the values of this variable from the source (this automatically applies scaling and offset)
                loop_val = variable[:] 
        
                # 1b. Get the attributes for this variable from source
                for attrname in variable.ncattrs():
                    loop_attr_source_values[attrname] = variable.getncattr(attrname)
               
                # 2a. Create the .nc variable with the proper SUMMA name
                # Inputs: variable name as needed by SUMMA; data type: 'float'; dimensions; no need for fill value, because thevariable gets populated in this same script
                dest.createVariable(name, 'f4', ('time','latitude','longitude'), fill_value = False, zlib=True, shuffle=True)
        
                # 3a. Select the attributes we want to copy for this variable, based on the dictionary defined before the loop starts
                loop_attr_copy_values = {use_this: loop_attr_source_values[use_this] for use_this in loop_attr_copy_these}
        
                # 3b. Copy the attributes FIRST, so we don't run into any scaling/offset issues
                dest[name].setncattr('missing_value',-999)
                dest[name].setncatts(loop_attr_copy_values)
        
                # 3c. Copy the data SECOND
                dest[name][:] = loop_val

# --- ERA5 derived variables
def add_derived_variables(file):
    
    '''Adds various derived meteorological variables to a netCDF file with existing ERA5 forcing data'''

    with nc4.Dataset(file, 'r+') as f:
        #f = derive_reflected_shortwave_radiation(f) # disabled, because on reflection we don't want to provide this variable after all
        #f = derive_net_radiation(f) # as above
        f = derive_vapor_pressure(f)
        f = derive_relative_humidity(f)
        f = derive_mean_wind_speed(f)

def make_nc_variable(file, name, units, values, long_name='', standard_name='', history='', dims='lat/lon'):

    '''Makes a netcdf4 variable with dimensions (time,latitude,longitude) in file'''

    # Assumptions
    # - Dimension in input file are 'time', 'latitude', 'longitude'
    # - No fill value specified
    # - Compression options zlib=True, shuffle=True are active

    # 1. Create the .nc variable
    if dims.lower() == 'lat/lon':
        file.createVariable(name,'f4',('time','latitude','longitude'), fill_value = False, zlib=True, shuffle=True)
    elif dims.lower() == 'hru':
        file.createVariable(name,'f4',('time','hru'), fill_value = False, zlib=True, shuffle=True)
    else:
        print(f'!!! Warning: make_nc_variable(): dimensions option {dims} not recognized. Exiting.')
        return file

    # 2. Set the attributes FIRST, so we don't run into any scaling/offset issues
    file[name].setncattr('missing_value',-999)
    file[name].setncattr('units',units)
    if long_name: file[name].setncattr('long_name',long_name)
    if standard_name: file[name].setncattr('standard_name', standard_name)

    # 3. Copy the data SECOND
    file[name][:] = values

    # 4. Update history
    if history:
        current_history = file.getncattr('History')
        new_history = f'{current_history}{history}'
        file.setncattr('History', new_history)
    
    return file

def compute_wind_direction(u,v):
    return (180 + 180/np.pi*np.arctan2(u,v)) % 360

def derive_wind_direction(file, u_wind='u', v_wind='v', new_name='phi', test=False, dims='lat/lon'):

    '''Adds wind direction to a netCDF4 file that cantains u and v wind components'''

    # Note: np.arctan2() expects inputs as (y,x), so we would expect inputs as (v,u) (https://numpy.org/doc/stable/reference/generated/numpy.arctan2.html)
    #  We deliberately need to swap these to get the correct angles:
    # See also: https://confluence.ecmwf.int/pages/viewpage.action?pageId=133262398
    #  Here, the excel function atan2 expects inputs as (x,y) and ECMWF swaps these to (v,u)
    if test:
        us   = [ 10,  10, -10, -10,  10, -10,   0,   0]
        vs   = [ 10, -10,  10, -10,   0,   0,  10, -10]
        phis = [225, 315, 135,  45, 270,  90, 180,   0] # known
        for u_test,v_test,phi_test in zip(us,vs,phis):
            new_value = compute_wind_direction(u_test,v_test)
            assert new_value == phi_test, f'result {new_value} does not match known result {phi_test}, for u = {u_test}, v = {v_test}'
        print('-- derive_wind_direction() test completed successfully')

    # 0. Check if the variable already exists
    if new_name in file.variables:
        print(f'!!! Warning: derive_wind_direction(): variable {new_name} already exists in file. Exiting.')
        return file

    # 1. Get the values of this variable from the source (this automatically applies scaling and offset)
    u = file.variables[u_wind][:]
    v = file.variables[v_wind][:]

    # 2. Compute the new values
    values_new = compute_wind_direction(u,v)
    
    # 3. Define other inputs
    unit_new = 'degrees'
    long_name = 'Meteorological wind direction computed from U and V wind components, indicating angle where wind is coming from with North at 0 degrees and increasing clock-wise '
    new_history = f' On {time.ctime(time.time())}: derive_wind_direction().'

    # 4. Make the variable
    # By default, cs.make_nc_variable() uses dimensions (time,latitude,longitude)
    #   because we have only used this for gridded netcdf files before.
    #   This wind direction derivation is performed after netcdf regridding to
    #   polygons, meaning we now have to deal with the original gridded (time,lat,lon)
    #   dimensions, but also with the lumped and distributed cases where dimensions
    #   are (time,hru)
    file = cs.make_nc_variable(file, new_name, unit_new, values_new, long_name=long_name, history=new_history, dims=dims)

    return file

def derive_reflected_shortwave_radiation(file, incoming_shortwave='msdwswrf', net_shortwave='msnswrf', new_name='reflected_sw'):

    '''Adds reflected shortwave radiation to a netCDF4 file that cantains incoming and net shortwave radiation'''

    # 1. Get the values of this variable from the source (this automatically applies scaling and offset)
    ds = file.variables[incoming_shortwave][:]
    ns = file.variables[net_shortwave][:]

    # 2. Compute the values 
    values_new = ds - ns
    
    # 3. Define other inputs
    unit_new = 'W m**-2'
    long_name = 'Mean surface reflected shortwave radiation flux, computed from ERA5 mean surface downward short-wave radiation flux and mean surface net short-wave radiation flux'
    new_history = f' On {time.ctime(time.time())}: derive_reflected_shortwave_radiation().'

    # 4. Make the variable
    file = make_nc_variable(file, new_name, unit_new, values_new, long_name=long_name, history=new_history)
    
    return file

def derive_net_radiation(file, net_shortwave='msnswrf', net_longwave='msnlwrf', new_name='net_radiation'):

    '''Adds net radiation to a netCDF4 file that contains net shortwave and net longwave radiation'''

    # 1. Get the values of this variable from the source (this automatically applies scaling and offset)
    ns = file.variables[net_shortwave][:]
    nl = file.variables[net_longwave][:]

    # 2. Compute the values 
    values_new = ns + nl
    
    # 3. Define other inputs
    unit_new = 'W m**-2'
    long_name = 'Mean surface net radiation flux, computed from ERA5 mean surface net short-wave radiation flux and mean surface net long-wave radiation flux'
    new_history = f' On {time.ctime(time.time())}: derive_net_radiation().'

    # 4. Make the variable
    file = make_nc_variable(file, new_name, unit_new, values_new, long_name=long_name, history=new_history)
    
    return file

def derive_relative_humidity(file, air_temperature='t', vapor_pressure='e', new_name='rh'):

    '''Adds relative humidity to a netCDF4 file containing vapor pressure and air temperature'''

    # 0. Constants
    Rv = 461 # [J K-1 kg-1]
    T0 = 273.15 # [K]
    e0 = 0.6113 # [kPa]
    Lv = 2.5*10**6 # [J kg-1], latent heat of vaporization for liquid water
    Ld = 2.83*10**6 # [J kg-1], latent heat of deposition for ice

    # 1. Get the values of this variable from the source (this automatically applies scaling and offset)
    e = file.variables[vapor_pressure][:]
    T = file.variables[air_temperature][:]

    # 2. Compute the values
    is_ice = T < 273.15
    is_liq = T >= 273.15    
    es = T*0 # initialize the variable
    es[is_ice] = e0*np.exp(Ld/Rv * (1/T0 - 1/T[is_ice]))
    es[is_liq] = e0*np.exp(Lv/Rv * (1/T0 - 1/T[is_liq]))   
    values_new = e / es
    
    # 3. Define other inputs
    unit_new = 'kPa kPa**-1'
    long_name = 'relative humidity, computed from ERA5 temperature and derived vapor pressure'
    new_history = f' On {time.ctime(time.time())}: derive_relative_humidity().'

    # 4. Make the variable
    file = make_nc_variable(file, new_name, unit_new, values_new, long_name=long_name, history=new_history)
    return file

def derive_vapor_pressure(file, air_pressure='sp', specific_humidity='q', Pa_to_kPa=True, new_name='e', approximate=False):

    '''Adds vapor pressure to a netCDF4 file containing air pressure and specific humidity'''

    # Assumptions:
    # - New variable's attribute values are hard-coded
    
    # 0. Constants
    Rd = 2.871*10**(-4) # [kPa K-1 m3 kg-1]
    Rv = 4.61*10**(-4) # [kPa K-1 m3 kg-1]
    epsilon = Rd/Rv # Should be ~ 0.622
    
    # 1. Get the values of this variable from the source (this automatically applies scaling and offset)
    P = file.variables[air_pressure][:]
    q = file.variables[specific_humidity][:]

    # 2. Convert air pressure units if flagged
    if Pa_to_kPa: P = P / 1000

    # 3. Compute the values
    values_new = -1*(q*P)/(q*epsilon - q - epsilon)
    if approximate: values_new = q*P/(epsilon)

    # 4. Define other inputs
    unit_new = 'kPa'
    long_name = 'vapor pressure, computed from ERA5 specific humidity and surface pressure'
    new_history = f' On {time.ctime(time.time())}: derive_vapor_pressure().'

    # 5. Make the variable
    file = make_nc_variable(file, new_name, unit_new, values_new, long_name=long_name, history=new_history)
    return file

def derive_mean_wind_speed(file, var_1='u', var_2='v', new_name='w'):
    
    '''Adds mean wind speed to a netCDF4 file containing U and V wind direction data'''

    # Assumptions:
    # - New variable's attribute values are hard-coded
    
    # 1. Get the values of this variable from the source (this automatically applies scaling and offset)
    values_1 = file.variables[var_1][:]
    values_2 = file.variables[var_2][:]

    # 2. Create the variable attribute 'units' from the source data. This lets us check if the source units match (they should match)
    unit_1 = file.variables[var_1].getncattr('units')
    unit_2 = file.variables[var_2].getncattr('units')
    assert unit_1 == unit_2, f'WARNING: units of source variables do not match: variable {var1} in {unit_1}, varialbe {var_2} in {unit_2}'
    unit_new = '(({})**2 + ({})**2)**0.5'.format(unit_1,unit_2)

    # 3. Compute the values
    values_new = ((values_1**2)+(values_2**2))**0.5
    
    # 4. Define other inputs
    long_name = 'mean wind speed at the measurement height, computed from ERA5 U and V component of wind'
    standard_name = 'wind_speed'
    new_history = f'. On {time.ctime(time.time())}: derive_mean_wind_speed().'

    # 5. Make the variable
    file = make_nc_variable(file, new_name, unit_new, values_new, long_name=long_name, standard_name=standard_name, history=new_history)
    return file

# --- EM-Earth processing
def compare_em_earth_netcdf_dimensions(p_file, t_file):

    '''Opens two netCDF files containing EM-Earth data and checks if contents of latitude, longitude and time match'''

    # Extract space and time info for checks
    with nc4.Dataset(p_file) as src1, nc4.Dataset(t_file) as src2: # implicitly closes files
        p_lat = src1.variables['lat'][:]
        p_lon = src1.variables['lon'][:]
        p_time = src1.variables['time'][:]
        t_lat = src2.variables['lat'][:]
        t_lon = src2.variables['lon'][:]
        t_time = src2.variables['time'][:]

    # Perform the check
    if not all( [all(p_lat == t_lat), all(p_lon == t_lon), all(p_time == t_time)] ):
        err_txt = (f'ERROR: Dimension mismatch while merging {pressure_file} and {surface_file}. '
                    'Check latitude, longitude and time dimensions in both files. Continuing with next files.')
        return False,err_txt
    else:
        return True,'Variables match'

def merge_em_earth_prcp_and_tmean_files(p_file, t_file, dest_file):

    '''Merges two EM-Earth data files into dest_file'''

    # Transfer definitions
    dimensions_t_transfer = ['lon','lat','time']
    variables_t_transfer = ['tmean']
    variables_p_transfer = ['prcp']
    attr_names_expected = ['scale_factor','add_offset','_FillValue','missing_value','units','long_name','standard_name'] # these are the attributes we think each .nc variable has
    loop_attr_copy_these = ['units','long_name','standard_name'] # we will define new values for _FillValue and missing_value when writing the .nc variables' attributes

    with nc4.Dataset(p_file) as src1, nc4.Dataset(t_file) as src2, nc4.Dataset(dest_file, 'w') as dest: # implicitly closes files

        # Set general attributes
        dest.setncattr('History','Created ' + time.ctime(time.time()))
        dest.setncattr('Reason','Merging separate EM-Earth download files')
        dest.setncattr('Source','github.com/cH-Earth/camels_spat')

        # Copy meta data from the two source files
        dest.setncattr('Conventions','CF-1.6') # Hopefully - EM-Earth does not specify these
        dest.setncattr('History_source_file_1',src1.getncattr('history'))
        dest.setncattr('History_source_file_2',src2.getncattr('history'))

        # --- Dimensions: latitude, longitude, time
        # Using src2 here because this makes adapting the code from the ERA5 function simpler
        # If there is a mismatch between surface and pressure we shouldn't have reached this point at all due to the earlier check
        for name, dimension in src2.dimensions.items():

            # Replace the names to be consistent with ERA5 dimensions - this is better for user friendliness
            if name.lower() == 'lat': name_new = 'latitude'
            if name.lower() == 'lon': name_new = 'longitude'
            if name.lower() == 'time': name_new = 'time' # just to assign 'time' to name_new so we can always use name_new below

            # Create the dimensions
            if dimension.isunlimited():
                dest.createDimension( name_new, None)
            else:
                dest.createDimension( name_new, len(dimension))

        # --- Get the tmean-file dimension variables (lat, lon, time)
        for name, variable in src2.variables.items():
    
            # Transfer lat, long and time variables because these don't have scaling factors
            if name in dimensions_t_transfer:

                # Replace the names to be consistent with ERA5 dimensions - this is better for user friendliness
                if name.lower() == 'lat': name_new = 'latitude'
                if name.lower() == 'lon': name_new = 'longitude'
                if name.lower() == 'time': name_new = 'time' # just to assign 'time' to name_new so we can always use name_new below

                new_dims = tuple(
                    'latitude' if x == 'lat' else 'longitude' if x == 'lon' else x
                    for x in variable.dimensions
                    )
                
                dest.createVariable(name_new, variable.datatype, new_dims, fill_value = -999)
                dest[name_new].setncatts(src2[name].__dict__)
                dest.variables[name_new][:] = src2.variables[name][:]

        # ---  Transfer the tmean-file data first, for no particular reason
        for name, variable in src2.variables.items():

            # Check that we are only using the names we expect, and thus the names for which we have the required code ready
            if name in variables_t_transfer:
        
                # 0. Reset the dictionary that we keep attribute values in
                loop_attr_source_values = {name: 'n/a' for name in attr_names_expected}
        
                # 1a. Get the values of this variable from the source (this automatically applies scaling and offset)
                loop_val = variable[:]
        
                # 1b. Get the attributes for this variable from source
                for attrname in variable.ncattrs():
                    loop_attr_source_values[attrname] = variable.getncattr(attrname)
               
                # 2. Create the .nc variable 
                # Inputs: variable name as needed by SUMMA; data type: 'float'; dimensions; no need for fill value, because thevariable gets populated in this same script
                dest.createVariable(name, 'f4', ('time','latitude','longitude'), fill_value = False, zlib=True, shuffle=True)
        
                # 3a. Select the attributes we want to copy for this variable, based on the dictionary defined before the loop starts
                loop_attr_copy_values = {use_this: loop_attr_source_values[use_this] for use_this in loop_attr_copy_these}
        
                # 3b. Copy the attributes FIRST, so we don't run into any scaling/offset issues
                dest[name].setncattr('missing_value',-999)
                dest[name].setncatts(loop_attr_copy_values)
        
                # 3c. Copy the data SECOND
                dest[name][:] = loop_val

        # --- Transfer the pressure level variables next, using the same procedure as above
        for name, variable in src1.variables.items():
            if name in variables_p_transfer:
        
                # 0. Reset the dictionary that we keep attribute values in
                loop_attr_source_values = {name: 'n/a' for name in attr_names_expected}
        
                # 1a. Get the values of this variable from the source (this automatically applies scaling and offset)
                loop_val = variable[:] 
        
                # 1b. Get the attributes for this variable from source
                for attrname in variable.ncattrs():
                    loop_attr_source_values[attrname] = variable.getncattr(attrname)
               
                # 2a. Create the .nc variable with the proper SUMMA name
                # Inputs: variable name as needed by SUMMA; data type: 'float'; dimensions; no need for fill value, because thevariable gets populated in this same script
                dest.createVariable(name, 'f4', ('time','latitude','longitude'), fill_value = False, zlib=True, shuffle=True)
        
                # 3a. Select the attributes we want to copy for this variable, based on the dictionary defined before the loop starts
                loop_attr_copy_values = {use_this: loop_attr_source_values[use_this] for use_this in loop_attr_copy_these}
        
                # 3b. Copy the attributes FIRST, so we don't run into any scaling/offset issues
                dest[name].setncattr('missing_value',-999)
                dest[name].setncatts(loop_attr_copy_values)
        
                # 3c. Copy the data SECOND
                dest[name][:] = loop_val

def update_em_earth_units(file):
    
    '''Updates units in existing EM_Earth netcdf files'''
    
    with nc4.Dataset(file, 'r+') as f:
        update_em_earth_p(f)
        update_em_earth_t(f)

def update_em_earth_t(file, temperature='tmean'):
    
    '''Updates temperature units from [degree C] to [K]'''
    
    var = file.variables[temperature]
    var[:] = var[:] + 273.15 # [C] + 273.15 = [K]
    var.setncattr('units', 'K')
    
    new_history = f' On {time.ctime(time.time())}: update_em_earth_t().'
    if 'History' in file.ncattrs():
        current_history = file.getncattr('History')
        new_history = f'{current_history}{new_history}'
    file.setncattr('History', new_history)
    
    return file

def update_em_earth_p(file, precipitation='prcp'):
    
    '''Updates precipitation units from [mm hr-1] to [kg m-2 s-1]'''
    
    var = file.variables[precipitation]
    var[:] = var[:] / 3600 # (x) [mm hr-1] * (0.001) [m mm-1] * (1000) [kg m-3] * (1/3600) [s hr-1] = [kg m-2 s-1]
    var.setncattr('units', 'kg m**-2 s**-1')
    
    new_history = f' On {time.ctime(time.time())}: update_em_earth_p().'
    if 'History' in file.ncattrs():
        current_history = file.getncattr('History')
        new_history = f'{current_history}{new_history}'
    file.setncattr('History', new_history)
    
    return file
