'''Contains functions to remap geospatial data onto (sub)basins.'''

import easymore
import geopandas as gpd
import netCDF4 as nc4
from pathlib import Path
import shapefile # this is installed as "PyShp" library
import sys
import xarray as xr

def check_remap_need(file, lat_var='latitude', lon_var='longitude'):
    
    '''Opens a netcdf file and checks if latitude and longitude dimensions are larger than 1'''
    
    ds = xr.open_dataset(file)
    needs_remap = False
    if (len(ds[lat_var]) > 1) or (len(ds[lon_var]) > 1): 
        needs_remap = True
    return needs_remap

def run_easymore_to_make_remap_file(nc_files,esmr_objects):

    '''Loops over two lists of input netCDF4 files and EASYMORE objects to create a remap file for each combination'''

    for nc_file,esmr in zip(nc_files,esmr_objects):
        esmr.source_nc = nc_file
        esmr.nc_remapper()
        #print(f'Mapping {nc_file} with {esmr}')
    return

def add_crs_to_shapefile(file, crs='EPSG:4326'):

    '''Adds CRS definition to a shapefile'''

    shp = gpd.read_file(file)
    shp.crs = crs
    shp.to_file(file)

    return

def get_easymore_settings(data, case, grid_shp, basin_shp, temp_folder, out_folder):

    '''Creates an EASYMORE object for remapping, with variable names specific to ERA5 or EM-Earth dataset'''

    # Initialize an EASYMORE object
    esmr = easymore.easymore()

    # General settings
    esmr.author_name = 'CAMELS-spat workflow'
    esmr.case_name = data + '_' + case # Case name, used in EASYMORE-generated file names
    esmr.save_csv  = False # Flag that we do not want the data stored in .csv in addition to .nc
    esmr.remap_csv = '' # Initialize with 'no remap file available' - we'll update this later
    esmr.sort_ID = False # Enforce that we want our HRUs returned in the order we put them in

    # Folders
    esmr.temp_dir = str(temp_folder) + '/' # Path() to string; ensure the trailing '/' EASYMORE wants
    esmr.output_dir = str(out_folder) + '/' # Path() to string; ensure the trailing '/' EASYMORE wants

    # Shapefiles
    esmr.source_shp = grid_shp
    esmr.source_shp_lat = 'latitude'
    esmr.source_shp_lon = 'longitude'
    
    esmr.target_shp = basin_shp
    if case == 'lumped': esmr.target_shp_ID  = 'FID' # name of the HRU ID field
    if case == 'dist': esmr.target_shp_ID  = 'COMID' # name of the HRU ID field
    #esmr.target_shp_lat = 'latitude'
    #esmr.target_shp_lon = 'longitude'

    # Input netcdf
    esmr.var_lat   = 'latitude'  # name of the latitude dimensions
    esmr.var_lon   = 'longitude' # name of the longitude dimension
    esmr.var_time  = 'time'      # name of the time dimension

    # Output netcdf
    esmr.remapped_dim_id = 'hru'     # name of the non-time dimension; prescribed by SUMMA
    esmr.remapped_var_id = 'hruId'   # name of the variable associated with the non-time dimension
    esmr.format_list     = ['f4']    # variable type to save forcing as. Single entry here will be used for all variables
    esmr.fill_value_list = ['-9999'] # fill value
    
    # Dataset-specific settings
    if data == 'ERA5':
        esmr.var_names = ['msdwlwrf', 'msnlwrf', 'msdwswrf', 'msnswrf', 'mtpr', 'sp',
                          'mper', 't', 'q', 'u', 'v', 'reflected_sw', 'net_radiation',
                          'e', 'rh', 'w'] # variable names of forcing data - hardcoded because we prescribe them during ERA5 merging
    if data == 'EM-Earth':
        esmr.var_names = ['tmean', 'prcp'] 

    # EASYMORE uses a default name for the remap CSV file based on esmr.case_name:
    remap_csv = esmr.case_name + '_remapping.csv'
    
    return esmr, remap_csv

def prepare_easymore_temp_folder(country,basin_id,data_path):

    '''Prepares temporary folder for EASYMORE-based remapping'''

    # Make folder name
    full_id = country + '_' + basin_id
    
    # Construct the paths
    main_folder = Path(data_path) / 'basin_data' / (country + '_' + basin_id) / 'forcing' / 'TEMP_easymore'
    main_folder.mkdir(parents=True, exist_ok=True)

    return main_folder

def prepare_forcing_grid_shapefiles(country,basin_id,data_path):

    '''Crates folder and file names for ERA5 and EM-Earth forcing grid shapefiles'''

    # Make folder name
    full_id = country + '_' + basin_id
    
    # Construct the paths
    main_folder = Path(data_path) / 'basin_data' / (country + '_' + basin_id) / 'shapefiles' / 'forcing_grids'
    main_folder.mkdir(parents=True, exist_ok=True)

    # Make the output file paths
    era5_grid_file = main_folder / ('ERA5_grid_' + full_id + '.shp')
    eme_grid_file = main_folder / ('EM_Earth_grid_' + full_id + '.shp')

    return era5_grid_file, eme_grid_file

def make_forcing_grid_shapefile(infile, outfile, ndec=4):

    '''Reads grid spacing from infile (netCDF4) and stores grid shapefile in outfile'''  

    # ndec: number of decimals to round shape spacing to

    # Open the file and get the dimensions
    with nc4.Dataset(infile, 'r') as src:
        lat = src.variables['latitude'][:]
        lon = src.variables['longitude'][:]
        
    # Find the spacing - round to a few decimals so we don't get issues with float precision
    if len(lat) == 1 and len(lon) == 1:
        # Special case: 1x1 grid cell in forcing file
        if 'ERA5' in str(infile):
            half_dlat = 0.25/2
            half_dlon = 0.25/2
        elif 'EM_Earth' in str(infile):
            half_dlat = 0.10/2
            half_dlon = 0.10/2
    elif (len(lat) == 1) != (len(lon) == 1):
        print(' WARNING: make_forcing_grid_shapefile(): special case of 1-by-X or X-by-1 forcing grid not implemented yet. Exiting.')
        sys.exit(0)
    else:
        half_dlat = round(abs(lat[1] - lat[0])/2, ndec)
        half_dlon = round(abs(lon[1] - lon[0])/2, ndec)

    # create a new shapefile object
    with shapefile.Writer(str(outfile)) as w:
        w.autoBalance = 1 # turn on function that keeps file stable if number of shapes and records don't line up
        w.field("ID",'N') # create (N)umerical attribute fields, integer
        w.field('latitude','F',decimal=ndec) # float
        w.field('longitude','F',decimal=ndec)
        ID = 0 # start ID counter of empty
        
        for i in range(0,len(lon)):
            for j in range(0,len(lat)):
                ID += 1
                center_lon = lon[i]
                center_lat = lat[j]
                vertices = []
                parts = []
                vertices.append([center_lon-half_dlon, center_lat])
                vertices.append([center_lon-half_dlon, center_lat+half_dlat])
                vertices.append([center_lon          , center_lat+half_dlat])
                vertices.append([center_lon+half_dlon, center_lat+half_dlat])
                vertices.append([center_lon+half_dlon, center_lat])
                vertices.append([center_lon+half_dlon, center_lat-half_dlat])
                vertices.append([center_lon          , center_lat-half_dlat])
                vertices.append([center_lon-half_dlon, center_lat-half_dlat])
                vertices.append([center_lon-half_dlon, center_lat])
                parts.append(vertices)
                w.poly(parts)
                w.record(ID, center_lat, center_lon)

    # Add the CRS
    add_crs_to_shapefile(outfile)

    return

def add_geopotential_to_era5_grid(invariant_file, grid_file, g = 9.80665):

    '''Reads geopotential values from invariant_file and stores in the grid shapefile'''

    # Open the files
    geo = xr.open_dataset( invariant_file ).isel(time=0)
    shp = gpd.read_file( grid_file )

    # Add new columns to shapefile
    shp = shp.assign(geop = -999)  # insert a placeholder value [m**2 s**-2]
    shp = shp.assign(elev_m = -999)  # insert a placeholder value [m]

    # For each row in the shapefile, match its ERA5 lat/lon coordinates 
    # with those in the 'geo' file and extract the appropriate geopotential
    for index, row in shp.iterrows():
           
        # Find elevation
        geop = geo['z'].sel(latitude = row['latitude'], longitude=row['longitude']).values.flatten()
        elev = geop / g
        
        # Add elevation into shapefile
        shp.at[index,'geop'] = geop[0]
        shp.at[index,'elev_m'] = elev[0]

    # Overwrite the existing shapefile
    shp.to_file( grid_file )
    geo.close()

    return