'''Contains functions to do basin delineation.'''

# File paths and names
# ------------------------------------------------------------------------------------------------------------------------
def prepare_delineation_outputs(df,i,data_path):
    
    '''Prepares output folders and file paths for lumped and distributed dilneation outcomes'''
    
    from pathlib import Path
    
    # Get identifiers
    country = df.iloc[i].Country
    basin_id = df.iloc[i].Station_id
    full_id = country + '_' + basin_id
    
    # Construct the paths
    main_folder = Path(data_path) / 'basin_data' / (country + '_' + basin_id) / 'shapefiles'
    lump_folder = main_folder / 'lumped'
    dist_folder = main_folder / 'distributed'
    refr_folder = main_folder / 'reference'
    
    # Make the paths
    lump_folder.mkdir(parents=True, exist_ok=True)
    dist_folder.mkdir(parents=True, exist_ok=True)
    
    # Make the output file paths
    lump_file = lump_folder / (full_id + '_lumped.shp')
    dist_file = dist_folder / (full_id + '_distributed_{}.shp')
    refr_file = refr_folder / (full_id + '_reference.shp')
    plot_file = main_folder / (full_id + '_delineation_results.png')
    
    return full_id, lump_file, dist_file, refr_file, plot_file

# Station data
# ------------------------------------------------------------------------------------------------------------------------
def read_delineation_coords(df,i):
    
    '''Reads the station or ourlet location from CAMELS-spat metadata file.'''
    
    # Use a manual location if provided
    if df['Manual_outlet_location'].iloc[i] == 'yes':
        lat,lon = df['Manual_lat'].iloc[i],df['Manual_lon'].iloc[i]
    # If not, use the mapped locations 
    elif df['Mapped_lat'].iloc[i] > -999 and df['Mapped_lon'].iloc[i] > -999:
        lat,lon = df['Mapped_lat'].iloc[i],df['Mapped_lon'].iloc[i]
    # If no mapped location is known (i.e. when we're about to do the mapping), preferentially use the outlet location if we have it
    elif df['Outlet_lat'].iloc[i] > -999 and df['Outlet_lon'].iloc[i] > -999:
        lat,lon = df['Outlet_lat'].iloc[i],df['Outlet_lon'].iloc[i]
    # If no outlet location is provided, use the station location instead
    else:
        lat,lon = df['Station_lat'].iloc[i],df['Station_lon'].iloc[i]
    
    return lat,lon

# Delineation prep
# ------------------------------------------------------------------------------------------------------------------------
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
    
    import os
    import os.path
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
        
        # If a temporary file already exists, remove it
        if os.path.isfile(outfile):
            os.remove(outfile)

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

def find_subbasin_containing_point(shp,point):
    
    '''Finds which polygon in a shape contains a given (lon,lat) coordinate,
    
    Inputs:
    - shp:   shapefile with polygon geometries
    - point: Shapely Point() with lon,lat coordinate
    
    Returns:
    - Masked shapefile showing the row of the polygon that contains point
    '''
    
    mask = shp.contains(point)
    if not any(mask):
        print('ERROR: find_subbasin_containing_point: shape does not contain point')
    
    return shp[mask]

def extract_shape_subset(large_shp, id_name, ids):
    
    '''Extracts a sub-domain from a larger shape.
    
    Inputs:
    - large_shp: shapefile containing the full domain, with at minimum an ID column
    - id_name:   name of the ID column in the shapefile
    - ids:       IDs to be extracted
    
    Returns:
    - Masked shapefile showing only the selected IDs.
    '''
    
    # Extract the shape
    small_shp = large_shp.loc[large_shp[id_name].isin(ids)]
        
    return small_shp

# Lumped basin delineation
# ------------------------------------------------------------------------------------------------------------------------
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
                    crs=from_epsg(4326), # This is what the Merit Hydro .tifs in the grid are in, but somehow this doesn't work automatically
                    schema=schema) as c:
        i = 0
        for shape, value in shapes:
            rec = {}
            rec['geometry'] = shape
            rec['properties'] = {'LABEL' : str(value)} # This also results in an empty column called 'LABEL' that we don't need, but without this line the shapefile writer trips up. Drop the column later
            rec['id'] = str(i)
            c.write(rec)
            i += 1
    
    # Open the shapefile, fix geometry, save and return
    shp = gpd.read_file(file_path)
    shp = fix_geom(shp)
    shp = shp.drop('LABEL', axis=1)
    shp.to_file(str(file_path))
    
    return shp

def add_area_to_shape(shp, column='area', crs='ESRI:102008'):
    
    '''Computes the area of shp using crs and stores in column'''
    
    areas = shp.to_crs(crs).area / 10**6 # [m^2] > [km^2]
    shp[column] = areas
    
    return shp

# Distributed basin delineation
# ------------------------------------------------------------------------------------------------------------------------
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

# Stats 
# ------------------------------------------------------------------------------------------------------------------------
def calculate_basin_and_reference_overlap(basin, ref_file, crs='ESRI:102800'):
    
    '''Calculates areal overlap between delineated basin and reference shape, if a reference shape exists.
    
    Input:
    - basin: GeoDataframe with first shapefile
    - ref_file: GeoDataframe with reference shapefile
    
    Optional input:
    - crs: Coordinate Reference System to perform area comparison in
    
    Return:
    - overlap: fractional overlap between both shapes
    '''
    
    import os.path
    import geopandas as gpd
    
    overlap = 'n/a'
    if os.path.isfile(ref_file):
        ref_shp = gpd.read_file(ref_file)
        basin = basin.dissolve() # Create a single polygon for intersection 
        overlap_1 = (ref_shp.intersection(basin, align=False).to_crs(crs).area / ref_shp.to_crs(crs).area)[0]
        overlap_2 = (basin.intersection(ref_shp, align=False).to_crs(crs).area / basin.to_crs(crs).area)[0]
        overlap = min(overlap_1,overlap_2)
    
    return overlap

def get_reference_areas(df,i):
    
    '''Reads reference areas from metadata file
    
    Input:
    - df: dataframe with CAMELS-spat metadata file
    - i: row index of basin under investigations
    
    Return:
    - out: dictionary with {Reference area source: reference area [km^2]}
    '''
    
    out = [[df['Ref_area_1_src'].iloc[i],   df['Ref_area_1_km2'].iloc[i]],
           [df['Ref_area_2_src'].iloc[i],   df['Ref_area_2_km2'].iloc[i]],
           [df['Ref_shape_source'].iloc[i], df['Ref_shape_area_km2'].iloc[i]]]
    
    return out

# Plotting
# ------------------------------------------------------------------------------------------------------------------------
def prepare_plotting_stats(ref_areas,area_lump,area_dist,overlap_lump,overlap_dist):
    
    '''Creates a list with statistics about the delineationthat is input to the main plotting function'''
    
    stats = ref_areas.copy()
    stats.append(['Lumped basin', area_lump])
    stats.append(['Distributed basin', area_dist])
    stats.append(['Lumped basin', overlap_lump])
    stats.append(['Distributed basin', overlap_dist])
    stats
    
    return stats

def prepare_plotting_legend(handles,labels,label,**kwargs):
    
    '''Adds a new legend item to handles and labels list to use with GeoPandas PAtchCollection plots.
    
    Inputs:
    - handles: list of existing handle objects (can be empty)
    - labels: list of existing legend labels (can be empty)
    - label: name of legend item to be added
    - **kwargs: keyword arguments for matplotlib.lines.Line2D. Useful here:
        - Marker shape and size:
            - 'marker'
            - 'markersize'
        - For color, use either:
            - 'color'
            - 'markerfacecolor' and 'markeredgecolor'
    
    Returns:
    - handles: updated list with handle objects
    - labels: updated list with labels
    '''
    
    from matplotlib.lines import Line2D
    
    # Create a legend item with the specified kwargs
    handles.append(Line2D([0],[0], linestyle="none", **kwargs))
    
    # Update the label
    labels.append(label)
        
    return handles, labels

def add_statistics_to_axis(ax,basin_id,stats):
    
    import numpy as np
    
    # Check inputs - 
    # if no reference shape exists, stats[5][1], stats[6][1] are 'n/a'
    if stats[5][1] == 'n/a': stats[5][1] = np.nan
    if stats[6][1] == 'n/a': stats[6][1] = np.nan 
    if stats[2][1] == -999:  stats[6][1] = np.nan 
    
    # Make a string
    txt = ('{}\n'
           'Area comparison {:>31}\n'
           'Ref 1: {:<32}: {:.2f}\n'
           'Ref 2: {:<32}: {:.2f}\n'
           'Ref 3: {:<32}: {:.2f}\n'
           '       {:<32}: {:.2f}\n'
           '       {:<32}: {:.2f}\n\n'
           'Fractional overlap with reference shape {:>4}\n'
           '       {:<32}: {:.2f}\n'
           '       {:<32}: {:.2f}\n'
           ''.format(basin_id,'[km^2]',
                     stats[0][0],stats[0][1],
                     stats[1][0],stats[1][1],
                     stats[2][0],stats[2][1],
                     stats[3][0],stats[3][1],
                     stats[4][0],stats[4][1],
                     '[-]',
                     stats[5][0],stats[5][1],
                     stats[6][0],stats[6][1]))
    
    # Write data onto axis                
    ax.text(-0.1,0.85,txt, va='top', transform=ax.transAxes, family='monospace', fontsize=10)
    
    return

def plot_discretization_results(basin_id, lump_shp, basin_shp, river_shp, ref_file, lat, lon, statistics_text, save_path, to_screen=False):
    
    '''Plots lumped and distributed discretization outcomes and reference shape (if available).
    
    Inputs:
    - basin_id: string to use as plot title, e.g. 'CAN_01AD002'
    - lump_shp: GeoDataframe with lumped basin outline
    - basin_shp: GeoDataframe with distributed basin outlines
    - river_shp: GeoDataframe with river network
    - ref_file: Path() to reference file. Not used if the specified file does not exist
    - lat: latitude of delineation outlet
    - lon: longitude of delineation outlet
    - save_path: Path() to location where figure should be saved
    
    Optional inputs:
    - to_screen: flag indicating if figure should be shown on screen or only saved to file. Default: False
    '''
    
    import os.path
    import geopandas as gpd
    import matplotlib.pyplot as plt
    
    # Handle reference shape
    have_ref = False
    if os.path.isfile(ref_file):
        have_ref = True
        ref_shp = gpd.read_file(ref_file)
    
    # Settings
    ms = 5 # marker size in legends
    
    # Colors suitable for all (?) kinds of color blindness
    # Source: https://colorbrewer2.org/#type=diverging&scheme=RdYlBu&n=5
    criv = [215/256, 25/256, 28/256] # rivers, red // This works best because it's the highest contrast with yellow & black
    temp = [253/256,174/256, 97/256] # unused, orange
    cfac = [255/256,255/256,191/256] # lumped & distributed shape facecolor, yellow 
    cbac = [171/256,217/256,233/256] # plot background, light blue
    cref = [ 44/120,123/256,182/256] # reference shape, dark blue 
    cedg = [  0/256,  0/256,  0/256] # lumped & distributed shape edgecolor, black
    
    # Create the plot
    fig, axs = plt.subplots(1,3,figsize=(15,7), width_ratios=[2,2,1])
    
    # Subplot 1: lumped shape
    # ---------------------------------------------------------
    ax = axs[0]
    ax.set_facecolor(cbac)
    
    # Shapes 
    lump_shp.plot(ax=ax, facecolor=cfac, edgecolor=cedg)
    if have_ref: ref_shp.boundary.plot(ax=ax, color=cref)
    ax.plot(lon,lat, marker='o', markeredgecolor='k', markerfacecolor='w', markersize=10)
    
    # Legend
    lines = []; labels = []
    lines,labels = prepare_plotting_legend(lines, labels, 'Lumped shape', marker='s', 
                                                                          markersize=ms, 
                                                                          markerfacecolor=cfac,
                                                                          markeredgecolor=cedg)
    if have_ref: lines,labels = prepare_plotting_legend(lines, labels, 'Reference shape', marker='_', 
                                                                                          markersize=ms,
                                                                                          markeredgewidth=2,
                                                                                          color=cref)
    lines,labels = prepare_plotting_legend(lines, labels, 'Outlet', marker='o', 
                                                                    markersize=ms, 
                                                                    markeredgecolor='k',
                                                                    markerfacecolor='w')
    ax.legend(lines, labels, loc='upper left')
    
    # Chart junk
    ax.set_title(basin_id)
    ax.set_xlabel('Longitude [degrees]')
    ax.set_ylabel('Latitude [degrees]')
    
    # Subplot 2: Distributed shape
    # ---------------------------------------------------------
    ax = axs[1]
    ax.set_facecolor(cbac)
    
    # Shapes
    basin_shp.plot(ax=ax, facecolor=cfac, edgecolor=cedg) 
    if len(river_shp) > 0: river_shp.plot(ax=ax, color=criv)
    if have_ref: ref_shp.boundary.plot(ax=ax, color=cref)
    ax.plot(lon,lat, marker='o', markeredgecolor='k', markerfacecolor='w', markersize=10)
    
    # Legend
    lines = []; labels = []
    lines,labels = prepare_plotting_legend(lines, labels, 'Distributed shape', marker='s', 
                                                                               markersize=ms, 
                                                                               markerfacecolor=cfac,
                                                                               markeredgecolor=cedg)
    lines,labels = prepare_plotting_legend(lines, labels, 'River network', marker='_', 
                                                                           markersize=ms,
                                                                           markeredgewidth=2,
                                                                           color=criv)    
    if have_ref: lines,labels = prepare_plotting_legend(lines, labels, 'Reference shape', marker='_', 
                                                                                          markersize=ms,
                                                                                          markeredgewidth=2,
                                                                                          color=cref)
    lines,labels = prepare_plotting_legend(lines, labels, 'Outlet', marker='o', 
                                                                    markersize=ms, 
                                                                    markeredgecolor='k',
                                                                    markerfacecolor='w')
    ax.legend(lines, labels, loc='upper left')
    
    # Chart junk
    ax.set_title(basin_id)
    ax.set_xlabel('Longitude [degrees]')
    
    # Subplot 3: text
    # ---------------------------------------------------------
    ax = axs[2]
    add_statistics_to_axis(ax,basin_id,statistics_text)
    ax.axis('off')
    
    # Save the plot
    # ---------------------------------------------------------
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    
    # Close if we don't want it displayed
    # ---------------------------------------------------------
    if not to_screen:
        plt.close()
    
    return # nothing