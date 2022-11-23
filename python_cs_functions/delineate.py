'''Contains functions to do basin delineation.'''

def read_delineation_coords(df,i):
    
    '''Reads the station or ourlet location from CAMELS-spat metadata file.'''
    
    # Preferentially use the outlet location if one is specified
    if df['Outlet_lat'].iloc[i] > -999 and df['Outlet_lon'].iloc[i] > -999:
        lat,lon = df['Outlet_lat'].iloc[i],df['Outlet_lon'].iloc[i]
    else:
        lat,lon = df['Station_lat'].iloc[i],df['Station_lon'].iloc[i]
    
    return lat,lon

def determine_pysheds_data_loading_window(lat,lon,src_file,
                                          lat_extent=10, lon_extent=15): 
    
    '''Defines a subset of the full data (flow direction, acumulation .tifs) grid to load.
    
    Inputs:
    - lat:  latitude of gauge location of interest
    - lon:  longitude of gauge location of interest
    - src_file: string with path and name of .tif file from which subsetting will occur. Neded to define outer bounds. 
    
    Optional inputs:
    - lat_extent: half length of window in latitude direction [degrees]. Default = 10 
    - lon_extent: half length of window in longitude direction [degrees]. Default = 15
    
    Note: optional input defaults determined through visual assessment of approximate basin outlines and gauge 
          locations of Canada's Reference Hydrometric Basin Network. The Canadian RHBN database includes larger
          basins than the CAMELS-US data, so these settings should work for both. 
    
    Returns:
    window: tuple of (lon_min,lat_min,lon_max,lat_max) with data window defined as 
            (lon-lon_extent, lat-lat_extent, lon+lon_extent, lat+lat_extent), limited by actual grid extent
    '''
    
    from osgeo import gdal
    
    # Open the source file and find the bounds
    data = gdal.Open(src_file)
    geoTransform = data.GetGeoTransform()
    minx = geoTransform[0]
    maxy = geoTransform[3]
    maxx = minx + geoTransform[1] * data.RasterXSize
    miny = maxy + geoTransform[5] * data.RasterYSize
    data = [] # implicitly close the file
    
    # Define the window
    window = (max(minx,lon-lon_extent),
              max(miny,lat-lat_extent),
              min(maxx,lon+lon_extent),
              min(maxy,lat+lat_extent))
    
    return window

def subset_tifs_around_gauge(data_window, acc_file, fdir_file,
                             temp_dir = './tmp_basin_delineation'):
    
    '''Subsets flow accumulation and direction .tifs to a smaller region.
    
    Inputs:
    - data_window: Subsetting region. (lon_min, lat_min, lon_max, lat_max)
    - acc_file: string with path to and name of full flow accumulation .tif file. '/path/to/file.tif'
    - fdir_file: string with path to and name of full flow direction .tif file. '/path/to/file.tif'
    
    Optional inputs:
    - temp_dir: location of the temporary directory where the smaller .tifs will be saved
    
    Returns:
    - temp_acc: path to and name of smaller flow accumulation tif
    - temp_fdir: path to and name of smaller flow direction tif
    '''
    
    from osgeo import gdal
    from pathlib import Path
        
    # Make the temporary directory
    temp_path = Path(temp_dir)
    temp_path.mkdir(parents=True, exist_ok=True)
    
    # Subset the .tifs
    small_files = ['small_acc.tif','small_fdir.tif'] # temp file names
    for ix,infile in enumerate([acc_file,fdir_file]):
    
        # Define where the small tif needs to go
        outfile = str(temp_path/small_files[ix])

        # Open the full tif
        ds = gdal.Open(infile)

        # Subset to datawindow and immediately write to temporary small file
        ds = gdal.Translate(outfile, ds, projWin = [data_window[0], # upper left x
                                                    data_window[3], # upper left y
                                                    data_window[2], # lower right x 
                                                    data_window[1]])# lower right y
    
        # Close
        ds = None
    
    # Return the temp locations
    temp_acc = temp_path/'small_acc.tif'
    temp_fdir = temp_path/'small_fdir.tif'
    
    return temp_acc,temp_fdir

def load_tifs_with_pysheds(acc_file,fdir_file):
    
    '''Loads .tifs to prepare PySheds usage.
    
    Inputs:
    - temp_acc: string with path to and file name of flow accumulation. '/path/to/file.tif'
    - temp_fdir: string with path to and file name of flow accumulation. '/path/to/file.tif'
    
    Returns:
    - grid: Grid() object containing spatial info
    - acc: flow accumulation grid
    - fdir: flow direction grid
    '''
    
    from pysheds.grid import Grid
    
    grid = Grid.from_raster( acc_file )
    acc  = grid.read_raster( acc_file )
    fdir  = grid.read_raster( fdir_file )
    
    return grid,acc,fdir

def get_merit_hydro_accumulated_upstream_area(grid,acc,lon,lat):
    
    '''Reads the accumulated upstream area for a given (lon,lat) location.
    
    Inputs:
    - grid: Pysheds grid object with grid definition
    - acc:  Pysheds accumulated area raster on 'grid'
    - lon:  longitude of location of interest
    - lat:  latitude of location of interest
        
    Returns:
    - accumulated upstream area at (lon,lat) in [km^2]
    '''
    
    # Find the grid closest to the pour point
    col,row = grid.nearest_cell(lon,lat, snap='center')
    
    # Extract values in the search range
    area = acc[row,col]    
    
    return area

def prepare_delineation_outputs(df,i,data_path):
    
    '''Prepares output folders and file paths for lumped and distributed dilneation outcomes'''
    
    from pathlib import Path
    
    # Get identifiers
    country = df.iloc[i].Country
    basin_id = df.iloc[i].Station_id
    
    # Construct the paths
    main_folder = Path(data_path) / 'basin_data' / (country + '_' + basin_id) / 'shapefiles'
    lump_folder = main_folder / 'lumped'
    dist_folder = main_folder / 'distributed'
    refr_folder = main_folder / 'reference'
    
    # Make the paths
    lump_folder.mkdir(parents=True, exist_ok=True)
    dist_folder.mkdir(parents=True, exist_ok=True)
    
    # Make the output file paths
    lump_file = lump_folder / ('lumped_' + basin_id +'.shp')
    dist_file = dist_folder / ('distributed_' + basin_id +'_{}.shp')
    refr_file = refr_folder / ('reference_' + basin_id +'.shp')
    plot_file = main_folder / ('delineation_results_' + basin_id + '.png')
    
    return lump_file, dist_file, refr_file, plot_file

def delineate_catchment_with_pysheds(grid, lon, lat, fdir, shapefile_path,
                                     snap='center'):
    
    '''Loads the necessary data to delineate a basin with Pysheds and saves the result as a shapefile.
    
    Inputs:
    - grid: Pysheds grid showing full data extent 
    - lon: Longitude of pour point
    - lat: Latitude of pour point
    - fdir: flow direction raster loaded with PySheds.
    - shapefile_path: path and name of where to save the resulting .shp. Example: /path/to/shape.shp
    
    Optional inputs:
    - snap: Snapping approach used by Pysheds. Default 'center'. Alternative 'corner'
    
    Note that the snapping option 'corner' delineates the basin from a grid that contains the provided 
       (lon,lat) point as the top-left corner of the pixel from which delineation starts. This "derivation"
       grid DOES NOT necessarily match your provided grid - unless this happens by chance!
       
    Returns:
    - Polygon containing delineated basin geometry
    '''
      
    # Delineate basin
    basin = grid.catchment(x=lon, y=lat, fdir=fdir, xytype='coordinate', snap=snap) 
        
    # Clip raster to cover basin only
    grid.clip_to(basin) 
    
    # Write the shape
    basin_shp = write_pysheds_grid_to_shp(grid, shapefile_path)
        
    return basin_shp

def write_pysheds_grid_to_shp(grid,file_path,
                              connectivity = 8, 
                              driver = 'ESRI Shapefile',
                              schema = {'geometry': 'Polygon',
                                        'properties': {'LABEL': 'float:16'}}):
    
    '''Writes a Pysheds grid showing a basin outline to a shapefile.
    
    Inputs:
    - grid: Pysheds grid containing basin outline
    - file_path: location and file name of where to save
    
    Optional inputs:
    - connectivity: Pixel connectivity (4 or 8 directions). Default 8 prevents "orphan" pixels that 
                    are not connected to main shape
    - driver: OGR format driver used to open file - see Fiona docs for more details (Section 1.4.1)
    - schema: Definition of data record - see Fiona docs for more detail (Section 1.4.1)
    
    Fiona docs: https://fiona.readthedocs.io/en/latest/manual.html
    
    Returns:
    - Resulting polygon
    '''
    
    import fiona
    from fiona.crs import from_epsg
    import geopandas as gpd
    
    # Convert grid to polygons
    shapes = grid.polygonize(connectivity=connectivity) 
    
    # Write shapefile
    with fiona.open(str(file_path), 'w',
                    driver=driver,
                    #crs=grid.crs.srs,
                    crs=from_epsg(4326), # This is what the Merit Hydro .tifs in the grid are in, but somehow this doesn't work
                    schema=schema) as c:
        i = 0
        for shape, value in shapes:
            rec = {}
            rec['geometry'] = shape
            rec['properties'] = {'LABEL' : str(value)}
            rec['id'] = str(i)
            c.write(rec)
            i += 1
    
    # Open the shapefile, fix geometry, save and return
    shp = gpd.read_file(file_path)
    shp = fix_geom(shp)
    shp.to_file(str(file_path))
    
    return shp

def fix_geom(in_feature):
    
    '''Fixes polygon geometries if invalid
       
       Source: https://stackoverflow.com/a/71231092
    '''
    
    from shapely.validation import make_valid
    
    # avoid changing original geodf
    in_feature = in_feature.copy(deep=True)    
        
    # drop any missing geometries
    in_feature = in_feature[~(in_feature.is_empty)]
    
    # Repair broken geometries
    for index, row in in_feature.iterrows(): # Looping over all polygons
        if row['geometry'].is_valid:
            next
        else:
            fix = make_valid(row['geometry'])

            try:
                in_feature.loc[[index],'geometry'] =  fix # issue with Poly > Multipolygon
            except ValueError:
                in_feature.loc[[index],'geometry'] =  in_feature.loc[[index], 'geometry'].buffer(0)
                
    return in_feature

def subset_merit_hydro_to_basin(basins, rivers, mask, shapefile_path, lat, lon,
                                id_column='COMID', area_column='unitarea', crs='ESRI:102008'):
    
    '''Subsets the larger MERIT Hydro shapes (basins and rivers) to the delineated lumped basin outline. 
    
    This /should/ work because both lumped and existing MERIT Hydro basin delineations are based on the 
    same DEM, and therefore lumped basin outline should map almost perfectly on the existing MERIT Hydro
    outlines. As an added benefit we automatically subset the most downstream MERIT Hydro sub-basin to
    end exactly at the station location.
    
    Inputs:
    - basins: GeoDataframe with MERIT Hydro basin delineation in EPSG:4326
    - rivers: GeoDataframe with MERIT Hydro river delineation in EPSG:4326
    - mask: GeoDataframe with with umped catchment outline delineated with Pysheds in EPSG:4326
    - shapefile_path: path and name of where to save the resulting .shp, with a wildcard to add basin/river specification.
                      Example: /path/to/shape_{}.shp
    - lat: outlet latitude. Needed for subsetting of larger MERIT Hydro basin GeoDataframe in cases where we don't have a river
    - lon: outlet longitude. See 'lat'.
    
    Optional inputs:
    - id_column: 'COMID'. Name of column that contains basin/river IDs in 'basins' and 'rivers' shapes respectively.
                 Assumed to be the same name in both geodataframes.
    - area_column: 'unitarea'. Name of column that contains basin area in 'basins'.
    - crs: 'ESRI:102008'. Coordinate Reference System to calculate basin area in. Default: Albers Equal Area Conic
    
    Returns:
    - GeoDataframe containing the distributed basin outlines
    - GeoDataframe containing the distributed river locations
    '''
    
    import geopandas as gpd
    
    # 1. Clip the river first, so that we can use its IDs to create a preliminary subset of the basins
    dist_river = gpd.clip(rivers, mask, keep_geom_type=True)
    
    # 2. Basins
    # Create a preliminary geodataframe with the correct MERIT Hydro basins, not yet accounting for the gauge location in 
    #  the most downstream basin. We do this instead of clipping the large basin shape to avoid a specific issue:
    #  Using gpd.clip() on the large shape returns the inner basin outlined by 'mask', but also very tiny intersections all
    #  along the 'mask' outline where 'mask' edges overlap with the edges of polygons outside 'mask'. By subsetting 'basins'
    #  before we clip we remove those polygons outside 'mask' and thus avoid this problem.
    #selection = basins[id_column].isin(dist_river[id_column]) # TO DO: handle case where no delineated river was available. Move the selection into a dedicated function for this
    selection = select_merit_basin_ids(basins, dist_river, id_column, lat, lon)
    dist_basin_full = basins[selection].copy().reset_index(drop=True) 
    
    # Clip to 'mask' so that the most downstream basin polygon ends at the gauge
    dist_basin = gpd.clip(dist_basin_full, mask, keep_geom_type=True)
    
    # Update the basin areas since we changed the size of the most downstream basin
    #  Do this before we dissolve, because multipolygons coming out of dissolve sometimes don't play nicely with .area
    new_area = dist_basin.to_crs(crs).area / 10**6 # Calculate areas in [km^2]
    dist_basin[area_column] = new_area
    
    # Clean up the clipped shape to return a single (multi)polygon per basin ID
    dist_basin = dist_basin.dissolve(by='COMID', aggfunc='sum').reset_index()     
    
    # Make sure the geometry is valid
    dist_basin = fix_geom(dist_basin)
    
    # 3. Save
    dist_river.to_file( str(shapefile_path).format('river') )
    dist_basin.to_file( str(shapefile_path).format('basin') )
    
    return dist_basin, dist_river

def select_merit_basin_ids(basins, rivers, column, lat, lon):
    
    '''Creates a boolean mask of 'basin' based on IDs present in field 'column' in 'river'. If no river exists, returns 
    a mask of 'basin' based on which polygon contains the point (lon,lat).
    
    Input:
    - basins: GeoDataframe with MERIT Hydro basin delineation in EPSG:4326
    - rivers: GeoDataframe with MERIT Hydro river delineation in EPSG:4326
    - column: Name of the column that has basin/river IDs respectively. Assumed to be identical in both GeoDataframes.
    - lat: outlet latitude 
    - lon: outlet longitude
    
    Return:
    - mask: Seleciton in 'basins' of IDs present in 'rivers'    
    '''
    
    from shapely.geometry import Point  
    
    # Handle not-headwater & headwater cases
    if len(rivers) > 0:
        # We have a delineated river, so at least 2 basins (i.e. not a headwater basin w/o discretized river network)
        mask = basins[column].isin(rivers[column])
    else:
        # Headwater basin with no discretized river network in the MERIT Hydro data. Use outlet location instead of IDs
        outlet = Point(lon,lat)
        mask = basins.contains(outlet)
        
    if not any(mask):
        print('ERROR: select_merit_basin_ids: could not identify basin IDs')
    
    return mask
