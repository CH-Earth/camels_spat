'''Contains functions to process reference shapes from either CAMELS-US data or Water Survey of Canada downloads'''

def process_camels_us_ref_shape(ref_shps, mask):
    
    '''Processes the CAMELS-US reference shape dataset to extract a single basin and prepare that for CAMELS-spat'''
    
    # Inputs:
    # ref_shapes: GeoPandas dataframe with 671 HCDN CAMELS-US basins
    # mask: boolean array from intersection ref_shapes['hru_id'] and CAMELS-spat basin ID
    #
    # Output
    # out: GeoPandas dataframe with the basin of interest, in EPSG:4326. Columns: 'Station_id', 'Area_km2'
    # src: string with reference area source
    
    out = ref_shps[mask].copy() # Make a copy, not a view
    out = out.reset_index(drop=True) # drop=True: don't include old index as new column
    out = out.drop(['ann_P','ave_T','july_T','Perimeter','elev_mean','lon_cen','lat_cen'], axis=1) # remove columns we don't want
    out['AREA'] = out['AREA'] / 10**6 # [m^2] to [km^2]
    out = out.rename(columns={'hru_id': 'Station_id', 'AREA': 'Area_km2'}) # Rename hru_id for consistency with CAMELS-spat meta
    out = out.to_crs('EPSG:4326') # CAMELS-spat standard
    out['Country'] = 'USA'
    out = out[['Country','Station_id','Area_km2','geometry']] # Re-order columns
    src = 'CAMELS-US data set (HCDN)'
    
    return out,src

def process_wsc2022_ref_shape(ref_shps,mask):
    
    '''Processes the WSC2022 reference shape dataset to extract a single basin and prepare that for CAMELS-spat'''
    
    # Inputs:
    # ref_shapes: GeoPandas dataframe with 1008 WSC2022 RHBN basins
    # mask: boolean array from intersection ref_shapes['StationNum'] and CAMELS-spat basin ID
    #
    # Output
    # out: GeoPandas dataframe with the basin of interest, in EPSG:4326. Columns: 'Station_id', 'Area_km2'
    # src: string with reference area source
    
    out = ref_shps[mask].copy() # Make a copy, not a view
    out = out.reset_index(drop=True) # drop=True: don't include old index as new column
    out = out.drop(['NameNom','Status','Etat','Aire_km2','Version','Date'], axis=1) # remove columns we don't want
    out = out.rename(columns={'StationNum': 'Station_id'}) 
    out = out.to_crs('EPSG:4326') # CAMELS-spat standard
    out['Country'] = 'CAN'
    out = out[['Country','Station_id','Area_km2','geometry']] # Re-order columns
    src = 'WSC 2022 data set'
    
    return out,src

def process_wsc2016_ref_shape(ref_shps,mask):
    
    '''Processes the WSC2016 reference shape dataset to extract a single basin and prepare that for CAMELS-spat'''
    
    # Inputs:
    # ref_shapes: GeoPandas dataframe with 876 WSC2016 RHBN basins
    # mask: boolean array from intersection ref_shapes['Station'] and CAMELS-spat basin ID
    #
    # Output
    # out: GeoPandas dataframe with the basin of interest, in EPSG:4326. Columns: 'Station_id', 'Area_km2'
    # src: string with reference area source
    
    out = ref_shps[mask].copy() # Make a copy, not a view
    out = out.reset_index(drop=True) # drop=True: don't include old index as new column
    out = out.drop(['StationNam','Stn_UID','Shp_Perime','Shape_Leng','Shape_Area','HydexArea'], axis=1) # remove columns we don't want
    out = out.rename(columns={'Station': 'Station_id', 'Shp_Area': 'Area_km2'}) 
    out = out.to_crs('EPSG:4326') # CAMELS-spat standard
    out['Country'] = 'CAN'
    out = out[['Country','Station_id','Area_km2','geometry']] # Re-order columns
    src = 'WSC 2016 data set'
    
    return out,src

def process_wsc2022_location_info(station_id, ref_shp):
    
    '''Processes the WSC2022 reference shapes to extract station and outlet locations for CAMELS-spat metadata'''
    
    # Inputs:
    # ref_shp: GeoPandas dataframe with 1008 WSC2022 RHBN locations
    #          We don't need separate functions to extract station and outlet info, because they have the same columns
    # station_id: 
    #
    # Output
    # lat: latitude coordinate of the Geometry in EPSG:4326
    # lon: longitude coordinate of the Geometry in EPSG:4326
    
    # Select the data
    mask = ref_shp['StationNum'] == station_id
    out = ref_shp[mask].copy() # Make a copy, not a view
    
    # Check the geometry
    if out.crs != 'EPSG:4326':
        out = out.to_crs('EPSG:4326')
    
    # Get the coordinates
    lon = ref_shp[mask].geometry.x
    lat = ref_shp[mask].geometry.y
    
    return lat,lon    