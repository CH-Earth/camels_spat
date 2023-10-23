'''Contains functions to remap geospatial data onto (sub)basins.'''

import easymore
import geopandas as gpd
import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
import os
from pathlib import Path
import shapefile # this is installed as "PyShp" library
from shapely.geometry import Polygon
import sys
import xarray as xr

def check_can_remap_as_is(file, lat_var='latitude', lon_var='longitude'):
    
    '''Opens a netcdf file and checks if latitude and longitude dimensions are larger than 1'''
    
    ds = xr.open_dataset(file)
    can_remap = False
    if (len(ds[lat_var]) > 2) and (len(ds[lon_var]) > 2):
        can_remap = True
    return can_remap

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
    
    # Define how we want to calculate the grid space
    if len(lat) == 1:
        if 'ERA5' in str(infile): half_dlat = 0.25/2
        if 'EM_Earth' in str(infile): half_dlat = 0.10/2
    else:
        half_dlat = round(abs(lat[1] - lat[0])/2, ndec)
        
    if len(lon) == 1:
        if 'ERA5' in str(infile): half_dlon = 0.25/2
        if 'EM_Earth' in str(infile): half_dlon = 0.10/2
    else:
        half_dlon = round(abs(lon[1] - lon[0])/2, ndec)
    
#    # Find the spacing - round to a few decimals so we don't get issues with float precision
#    if len(lat) == 1 and len(lon) == 1:
#        # Special case: 1x1 grid cell in forcing file
#        if 'ERA5' in str(infile):
#            half_dlat = 0.25/2
#            half_dlon = 0.25/2
#        elif 'EM_Earth' in str(infile):
#            half_dlat = 0.10/2
#            half_dlon = 0.10/2
#    elif (len(lat) == 1) != (len(lon) == 1):
#        print(' WARNING: make_forcing_grid_shapefile(): special case of 1-by-X or X-by-1 forcing grid not implemented yet. Exiting.')
#        sys.exit(0)
#    else:
#        half_dlat = round(abs(lat[1] - lat[0])/2, ndec)
#        half_dlon = round(abs(lon[1] - lon[0])/2, ndec)

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

def add_empty_grid_cells_around_single_cell_netcdf(file, 
                                                   grid_spacing=0.25, lat_dim='latitude', lon_dim='longitude', tim_dim='time', 
                                                   to_file=''):

    '''Takes an existing netCDF4 file with latitude and longitude dimensions of size 1, and adds np.nan cells around this'''
    
    with xr.open_dataset(file) as ds:

        # Get existing latitude and longitude values, and define the extended coordinates
        old_lat = ds[lat_dim].values
        old_lon = ds[lon_dim].values

        # Derive new dimensions from existing ones
        grids_lat = len(old_lat) + 2
        grids_lon = len(old_lon) + 2

        # Make the new coordinates
        new_lat = np.linspace(min(old_lat) - grid_spacing, max(old_lat) + grid_spacing, grids_lat)
        new_lon = np.linspace(min(old_lon) - grid_spacing, max(old_lon) + grid_spacing, grids_lon)
        
        # Workaround - 
        # With the code above (np.linspace), we sometimes get floating point errors (i.e., longitude might be 
        # -66.55000000001 rather than the -66.55 we expect. This leads to issues below, where we want to match
        # time, lat and lon coordinates in the old DataSet (ds) with those in the new DataSet (new_ds), where
        # ds contains the "exact" values and new_ds the ones with floating point errors. The rounding below
        # (hopefully) ensures that we get the "exact" values in new_ds so the matching works, but this is likely
        # quite brittle. Don't rely on it too much.
        # Discovered this while running the code without the rounding. It seems to work as expected in most/sometimes
        # cases, but in some others we run into this floating point error. No idea what makes these cases different. 
        new_lat = np.round(new_lat, 2)
        new_lon = np.round(new_lon, 2)
   
        # Create a grid_cells-by-grid_cells-by-time set of NaNs
        new_data = np.empty( (len(ds[tim_dim]), grids_lat, grids_lon) ) # Note that dim order must match what's used below
        new_data[:] = np.nan
    
        # Create a new xarray dataset with the expanded grid and data
        new_ds = xr.Dataset( coords={lon_dim: new_lon, lat_dim: new_lat, tim_dim: ds[tim_dim]})
        for variable in ds.variables:
            if variable not in [lat_dim, lon_dim, tim_dim]:
                new_var = xr.DataArray(new_data, dims=(tim_dim, lat_dim, lon_dim), name=variable) # Note that order of dimensions here must match NaNs above
                new_ds[variable] = new_var.copy() # This copy() is ESSENTIAL - without it, whenever you update the variable values below, you update these values WHEREVER this new_var is inserted in the dataset, i.e., everywhere
                new_ds[variable].attrs = ds[variable].attrs
                new_ds[variable].loc[ {tim_dim: ds[tim_dim], lon_dim: ds[lon_dim], lat_dim: ds[lat_dim]} ] = ds[variable][:] # Copy existing values
                #assert all(new_ds[variable].isel(latitude=1,longitude=1) == ds[variable]), f'{variable} values not correctly copied' # slow
                assert (new_ds[variable].sel(longitude=ds['longitude'], latitude=ds['latitude'], time=ds['time']) - ds[variable]).sum() == 0,\
                    f'{variable} values not correctly copied' # faster

        if to_file:
            new_ds.to_netcdf(to_file)
            return
        else:    
            return new_ds

def easymore_workflow(data, case, esmr_temp, grid_shp, basin_shp, out_folder, infiles):

    ''' Container for repeated tasks needed to run EASYMORE for ERA5 or EM-Earth inputs'''

    # Initiate EASYMORE object
    esmr, remap_file = get_easymore_settings(data, case, grid_shp, basin_shp, esmr_temp, out_folder)

    # Create the remap file
    run_easymore_to_make_remap_file([infiles[0]], [esmr])

    # Update the EASYMORE object now we have remap files
    esmr.remap_csv = str(esmr_temp / remap_file)

    # Remap the remaining files
    for file in infiles[1:]:
        esmr.source_nc = file
        esmr.nc_remapper()

    return esmr # EASYMORE object, containing all info used during run

def easymore_workflow_with_cell_padding(data, case, esmr_temp, grid_shp, basin_shp, out_folder, infiles, grid_spacing=0.25):

    '''Container to add a step to easymore_workflow() where we apply padding to input netCDF file(s)'''
    
    # Initiate EASYMORE object
    esmr, remap_file = get_easymore_settings(data, case, grid_shp, basin_shp, esmr_temp, out_folder)

    # Loop over the files and process
    for ix, file_nc in enumerate(infiles):

        # Pad the input file we're about to use
        temp_nc = str( esmr_temp / os.path.basename(file_nc) )
        if not os.path.isfile(temp_nc): # This lets us skip over making the temporary file if we've already done so for the lumped case; Needed due to weird file closing issues
            add_empty_grid_cells_around_single_cell_netcdf(file_nc, grid_spacing=grid_spacing, to_file=temp_nc)

        # See if we need to do the one-off processing
        if ix == 0:
            
            # Create a new temporary shapefile for the padded netcdf
            temp_shp = str(esmr_temp/'padded.shp')
            make_forcing_grid_shapefile(temp_nc, temp_shp)
            esmr.source_shp = temp_shp

            # Create the remap file
            esmr.source_nc = temp_nc
            esmr.nc_remapper()
            
            # Update the EASYMORE object now we have remap files
            esmr.remap_csv = str(esmr_temp / remap_file)
        
        else:
            esmr.source_nc = temp_nc
            esmr.nc_remapper()

    # Insert the source file into the EASYMORE object for plotting
    esmr.source_nc = file_nc
    esmr.source_shp = grid_shp

    return esmr # EASYMORE object, containing all info used during run

def era5_eme_easymore_plotting_loop(esmr_list, temp_folder, save_here=''):

    '''Loops over ERA5 and EM-Earth EASYMORE objects to plot'''

    fig,axs = plt.subplots(2,2,figsize=(20,20))
    axs = axs.ravel()
    plt.rcParams['font.size'] = 16
    for ix, esmr, plot_var in zip([0,1,2,3],
                                  esmr_list,
                                  ['t','t','tmean','tmean']):
        ax = axs[ix]
        plot_averaging_check(esmr,temp_folder, plot_var, ax=ax)
    plt.tight_layout()
    plt.savefig(save_here, bbox_inches='tight')
    plt.close()

def plot_averaging_check(esmr, temp_folder, plot_var, ax='', save_here=''):

    '''Takes an EASYMORE object and plots a summary of results'''

    # Construct the names we don't yet have
    int_name = esmr.temp_dir + esmr.case_name + '_intersected_shapefile.shp'
    yyyy_mm = os.path.basename(esmr.source_nc).split('_')[-1].replace('.nc','')
    out_name = esmr.output_dir + esmr.case_name + '_remapped_' + yyyy_mm + '-01-00-00-00.nc'
    
    # Load the files we need
    int = gpd.read_file(int_name)
    src = gpd.read_file(esmr.target_shp)
    with xr.open_dataset(esmr.source_nc).isel(time=0) as raw, xr.open_dataset(out_name).isel(time=0) as out: # Get a single timestep
    
        # Map new netcdf values onto source shapefile
        out_df = out[[plot_var, esmr.remapped_var_id]].to_dataframe()
        src = src.merge(out_df, left_on=esmr.target_shp_ID, right_on=esmr.remapped_var_id)
    
        # Get plotting min/max values
        vmin = raw[plot_var].min().values
        vmax = raw[plot_var].max().values
        if vmin == vmax:
            vmin *= 0.9
            vmax *= 1.1
    
        # Figure out if we need to handle the case where 'raw' has a spatial dimension of length 1
        raw_flag = False
        if len(raw['latitude']) == 1 or len(raw['longitude']) == 1:
            raw_flag = True
            new_shp = gpd.read_file(esmr.source_shp)
            new_shp[plot_var] = raw[plot_var].values.flatten()
            lbl = f'{raw[plot_var].long_name} [{raw[plot_var].units}]'
            #new_shp,lbl = create_square_geodataframe_from_netcdf(raw,esmr.case_name,plot_var)
        
        # figure
        cb_fmt = '%0.2f'
        if not ax: ax = plt.subplot() # make plotting axis if none was provided
        if raw_flag:
            new_shp.plot(ax=ax, edgecolor='None', column=plot_var, vmin=vmin, vmax=vmax, 
                         legend=True, legend_kwds={'label': lbl, 'orientation': 'vertical', 'shrink': 0.75, 'format': cb_fmt})
        else:
            subplot = raw[plot_var].plot.pcolormesh(ax=ax, vmin=vmin, vmax=vmax, cbar_kwargs={'shrink': 0.75, 'format': cb_fmt}) # plot netcdf values
        src.plot(ax=ax, edgecolor='None', column=plot_var, vmin=vmin, vmax=vmax) # plot source shapefile (color=new values)
        int.plot(ax=ax, color='None',edgecolor='k') # plot intersected shapefile (borders only)
        ax.set_title(esmr.case_name)
        if save_here: plt.save_fig(save_here)

