import baseflow
import geopandas as gpd
import glob
import numpy as np
import os
import pandas as pd
import rasterio
import xarray as xr
import warnings
from dateutil.relativedelta import relativedelta
from osgeo import gdal, osr
from rasterstats import zonal_stats
from scipy.stats import circmean, circstd, skew, kurtosis
from scipy.optimize import curve_fit

## ------- Collection functions
def attributes_from_soilgrids(geo_folder, dataset, shp_str, l_values, l_index, case='lumped'):

    '''Calculates attributes from SOILGRIDS maps'''

    # File specifiction
    sub_folders = ['bdod',     'cfvo',       'clay',    'sand',    'silt',    'soc',      'porosity']
    units       = ['cg cm^-3', 'cm^3 dm^-3', 'g kg^-1', 'g kg^-1', 'g kg^-1', 'dg kg^-1', '-']
    depths = ['0-5cm', '5-15cm', '15-30cm', '30-60cm', '60-100cm', '100-200cm']
    fields = ['mean', 'uncertainty']
    stats = ['mean', 'min', 'max', 'std']
    
    # Loop over the files and calculate stats
    for sub_folder, unit in zip(sub_folders, units):
        for depth in depths:
            for field in fields:
                tif = str(geo_folder / dataset / 'raw' / f'{sub_folder}' / f'{sub_folder}_{depth}_{field}.tif')
                if not os.path.exists(tif): continue # porosity has no uncertainty field
                zonal_out = zonal_stats(shp_str, tif, stats=stats, all_touched=True)
                scale,offset = read_scale_and_offset(tif)
                l_values = update_values_list(l_values, stats, zonal_out, scale, offset, case=case)
                l_index += [('Soil', f'{sub_folder}_{depth}_{field}_min',  f'{unit}', 'SOILGRIDS'),
                            ('Soil', f'{sub_folder}_{depth}_{field}_mean', f'{unit}', 'SOILGRIDS'),
                            ('Soil', f'{sub_folder}_{depth}_{field}_max',  f'{unit}', 'SOILGRIDS'),
                            ('Soil', f'{sub_folder}_{depth}_{field}_std',  f'{unit}', 'SOILGRIDS')]
           
    # --- Specifc processing for the conductivity fields
    # For conductivity we want to have a harmonic mean because that should be more
    # representative as a spatial average. rasterstats doesn't have that, so we 
    # need a custom function. We're splitting out the processing from the rest for
    # clarity. Could also have done this with a bunch of if-statements, but this 
    # seems cleaner, also because conductivity doesn't have uncertainty maps.
    
    # Process the data fields
    sub_folder = 'conductivity'
    unit       = 'cm hr^-1'
    depths = ['0-5cm', '5-15cm', '15-30cm', '30-60cm', '60-100cm', '100-200cm']
    field  = 'mean' # no uncertainty maps for these
    stats = ['min', 'max', 'std']

    for depth in depths:
        tif = str(geo_folder / dataset / 'raw' / f'{sub_folder}' / f'{sub_folder}_{depth}_{field}.tif')
        if not os.path.exists(tif): continue # porosity has no uncertainty field
        zonal_out = zonal_stats(shp_str, tif, stats=stats, add_stats={'harmonic_mean': harmonic_mean}, all_touched=True)
        scale,offset = read_scale_and_offset(tif)
        l_values = update_values_list(l_values, ['min', 'max', 'std', 'harmonic_mean'], zonal_out, scale, offset, case=case)
        l_index += [('Soil', f'{sub_folder}_{depth}_{field}_min',  f'{unit}', 'SOILGRIDS'),
                    ('Soil', f'{sub_folder}_{depth}_{field}_harmonic_mean', f'{unit}', 'SOILGRIDS'),
                    ('Soil', f'{sub_folder}_{depth}_{field}_max',  f'{unit}', 'SOILGRIDS'),
                    ('Soil', f'{sub_folder}_{depth}_{field}_std',  f'{unit}', 'SOILGRIDS')]

    return l_values, l_index


def attributes_from_glhymps(geo_folder, dataset, basin_shp_path, l_values, l_index, equal_area_crs='ESRI:102008'):

    '''Calculates attributes from GLHYMPS'''

    # We need some special actions if we're dealing with the distributed case, so check that first
    case = 'lumped'
    if 'distributed' in basin_shp_path:
        case = 'distributed'
    
    # Load the geology file
    geol_str = str(geo_folder / dataset / 'raw' / 'glhymps.shp')
    geol = gpd.read_file(geol_str)

    # -- File updates
    # Rename the columns, because the shortened ones are not very helpful
    if ('Porosity' in geol.columns) and ('Permeabili' in geol.columns) \
    and ('Permeabi_1' in geol.columns) and ('Permeabi_2' in geol.columns):
        geol.rename(columns={'Porosity': 'porosity',   'Permeabili': 'logK_Ferr',
                             'Permeabi_1': 'logK_Ice', 'Permeabi_2': 'logK_std'}, inplace=True)
        geol.to_file(geol_str)

    # Clean up a few unsightly processing errors
    geol['Shape_Area'] = geol.to_crs(equal_area_crs).area
    geol = geol[['IDENTITY_','Shape_Area','porosity','logK_Ferr','logK_Ice','logK_std','geometry']]
    geol.to_file(geol_str)
    # -- End file updates

    # Ensure we have the correct areas to work with
    if 'New_area_m2' not in geol.columns:
        geol['New_area_m2'] = geol.to_crs(equal_area_crs).area
        #geol.to_file(geol_str)
    
    # Now handle the different cases
    if case == 'lumped':
    
        # Create areal averages and standard deviations
        # Stdev source: https://stats.stackexchange.com/a/6536
        porosity_mean = (geol['porosity']*geol['New_area_m2']).sum()/geol['New_area_m2'].sum()
        num_obs_coef = ((geol['New_area_m2'] > 0).sum()-1)/(geol['New_area_m2'] > 0).sum()
        porosity_std = np.sqrt((geol['New_area_m2'] * (geol['porosity']-porosity_mean)**2).sum() / geol['New_area_m2'].sum())
        
        # Porosity
        l_values.append( geol['porosity'].min() )
        l_index.append(('Geology', 'porosity_min',  '-', 'GLHYMPS'))
        l_values.append( porosity_mean ) # areal average
        l_index.append(('Geology', 'porosity_mean',  '-', 'GLHYMPS'))
        l_values.append( geol['porosity'].max() )
        l_index.append(('Geology', 'porosity_max',  '-', 'GLHYMPS'))
        l_values.append( porosity_std )
        l_index.append(('Geology', 'porosity_std',  '-', 'GLHYMPS'))
    
        # Create areal averages and standard deviations
        # Stdev source: https://stats.stackexchange.com/a/6536
        permeability_mean = (geol['logK_Ice']*geol['New_area_m2']).sum()/geol['New_area_m2'].sum()
        num_obs_coef = ((geol['New_area_m2'] > 0).sum()-1)/(geol['New_area_m2'] > 0).sum()
        permeability_std = np.sqrt((geol['New_area_m2'] * (geol['logK_Ice']-permeability_mean)**2).sum() / geol['New_area_m2'].sum())
        
        # Permeability
        l_values.append( geol['logK_Ice'].min() )
        l_index.append(('Geology', 'log_permeability_min',  'm^2', 'GLHYMPS'))
        l_values.append( permeability_mean )
        l_index.append(('Geology', 'log_permeability_mean',  'm^2', 'GLHYMPS'))
        l_values.append( geol['logK_Ice'].max() )
        l_index.append(('Geology', 'log_permeability_max',  'm^2', 'GLHYMPS'))
        l_values.append( permeability_std )
        l_index.append(('Geology', 'log_permeability_std',  'm^2', 'GLHYMPS'))

        # Set the remaining output
        l_comids_glhymps = None

    # Case 2
    elif case == 'distributed':
    
        # Load the basin shape 
        basin = gpd.read_file(basin_shp_path)

        # Loop over the individual polygons and create a new mini-HydroLAKES geodataframe for each polygon
        num_poly = len(basin)
        l_comids_glhymps = basin['COMID'].values
        for i_poly in range(num_poly):

            # Rest the storage lists
            tmp_values = []
            tmp_index = []

            # Subset the shape to the subbasin
            poly_glhymps = subset_glhymps_to_subbasin(geol, basin.iloc[i_poly]['geometry'], equal_area_crs)

            # Create areal averages and standard deviations
            # Stdev source: https://stats.stackexchange.com/a/6536
            porosity_mean = (poly_glhymps['porosity']*poly_glhymps['New_area_m2']).sum()/poly_glhymps['New_area_m2'].sum()
            num_obs_coef = ((poly_glhymps['New_area_m2'] > 0).sum()-1)/(poly_glhymps['New_area_m2'] > 0).sum()
            porosity_std = np.sqrt((poly_glhymps['New_area_m2'] * (poly_glhymps['porosity']-porosity_mean)**2).sum() / poly_glhymps['New_area_m2'].sum())
            
            # Porosity
            tmp_values.append( poly_glhymps['porosity'].min() )
            tmp_index.append(('Geology', 'porosity_min',  '-', 'GLHYMPS'))
            tmp_values.append( porosity_mean ) # areal average
            tmp_index.append(('Geology', 'porosity_mean',  '-', 'GLHYMPS'))
            tmp_values.append( poly_glhymps['porosity'].max() )
            tmp_index.append(('Geology', 'porosity_max',  '-', 'GLHYMPS'))
            tmp_values.append( porosity_std )
            tmp_index.append(('Geology', 'porosity_std',  '-', 'GLHYMPS'))
        
            # Create areal averages and standard deviations
            # Stdev source: https://stats.stackexchange.com/a/6536
            permeability_mean = (poly_glhymps['logK_Ice']*poly_glhymps['New_area_m2']).sum()/poly_glhymps['New_area_m2'].sum()
            num_obs_coef = ((poly_glhymps['New_area_m2'] > 0).sum()-1)/(poly_glhymps['New_area_m2'] > 0).sum()
            permeability_std = np.sqrt((poly_glhymps['New_area_m2'] * (poly_glhymps['logK_Ice']-permeability_mean)**2).sum() / poly_glhymps['New_area_m2'].sum())
            
            # Permeability
            tmp_values.append( poly_glhymps['logK_Ice'].min() )
            tmp_index.append(('Geology', 'log_permeability_min',  'm^2', 'GLHYMPS'))
            tmp_values.append( permeability_mean )
            tmp_index.append(('Geology', 'log_permeability_mean',  'm^2', 'GLHYMPS'))
            tmp_values.append( poly_glhymps['logK_Ice'].max() )
            tmp_index.append(('Geology', 'log_permeability_max',  'm^2', 'GLHYMPS'))
            tmp_values.append( permeability_std )
            tmp_index.append(('Geology', 'log_permeability_std',  'm^2', 'GLHYMPS'))
            

            # Add the new values for this subbasin into the main l_values list
            l_values[i_poly] = tmp_values # This works because we create a unique dataframe just for HydroLAKES results

        # End of subbasin loop, add the new index entries only once
        l_index = tmp_index # we don't use append() in the distributed case because we'll be storing this in a dedicated HydroLAKES attribute df
    
    return l_values, l_index, l_comids_glhymps

def attributes_from_pelletier(geo_folder, dataset, shp_str, l_values, l_index, case='lumped'):
    
    '''Calculates statistics for Pelletier maps'''
    # See: https://daac.ornl.gov/SOILS/guides/Global_Soil_Regolith_Sediment.html

    # General
    files = ['upland_hill-slope_regolith_thickness.tif',
             'upland_hill-slope_soil_thickness.tif',
             'upland_valley-bottom_and_lowland_sedimentary_deposit_thickness.tif',
             'average_soil_and_sedimentary-deposit_thickness.tif']
    attrs = ['regolith_thickness',    'soil_thickness',
             'sedimentary_thickness', 'average_thickness']
    units = ['m', 'm', 'm', 'm']
    stats = ['mean', 'min', 'max', 'std']

    for file,att,unit in zip(files,attrs,units):
        tif = str( geo_folder / dataset / 'raw' / file )
        check_and_set_nodata_value_on_pelletier_file(tif)
        zonal_out = zonal_stats(shp_str, tif, stats=stats, all_touched=True)
        zonal_out = check_zonal_stats_outcomes(zonal_out, new_val=0) # certain tifs have no data because there is no variable X, so we set these to 0
        scale,offset = read_scale_and_offset(tif)
        l_values = update_values_list(l_values, stats, zonal_out, scale, offset, case=case)
        l_index += [('Soil', f'{att}_min',  f'{unit}', 'Pelletier'),
                    ('Soil', f'{att}_mean', f'{unit}', 'Pelletier'),
                    ('Soil', f'{att}_max',  f'{unit}', 'Pelletier'),
                    ('Soil', f'{att}_std',  f'{unit}', 'Pelletier')]
        
    return l_values,l_index

def check_and_set_nodata_value_on_pelletier_file(tif, nodata=255):
    if not get_geotif_nodata_value(tif): # no no-data value set yet
        data = get_geotif_data_as_array(tif)
        if data.max() == nodata: # but we do need to set a new no-data value here
            data = None
            ds = gdal.Open(tif)
            band = ds.GetRasterBand(1)
            band.SetNoDataValue(255)
            ds.FlushCache()
            ds = None
            return

def attributes_from_streamflow(hyd_folder, dataset, basin_id, pre, row, l_values, l_index):
    
    '''Calculates various streamflow signatures'''
    
    # Constants
    seconds_per_minute = 60 # s min-1
    seconds_per_hour = 60 * seconds_per_minute # s hr-1
    seconds_per_day = 24 * seconds_per_hour # s d-1
    water_density = 1000 # kg m-3
    mm_per_m = 1000 # mm m-1
    m_per_km = 1000 # m km-1
    
    # Find the data source
    source = 'USGS/WSC'
    
    # Load observations
    hyd_file = hyd_folder / f'{basin_id}_daily_flow_observations.nc'
    hyd = xr.open_dataset(hyd_file)
    hyd = subset_dataset_to_max_full_years(hyd, res='day', water_year=True)
    
    # Convert observations to mm d-1
    area_km2 = row['Basin_area_km2']
    area_m2 = area_km2 * m_per_km**2 # [km2] * ([m km-1]^2 = [m2 km-2]) = [m2]
    hyd['q_obs'] = hyd['q_obs'] * seconds_per_day / area_m2 * mm_per_m # [m3 s-1] * [s d-1] / [m2] * [mm m-1] = [m d-1]
    
    # Convert precipitation into mm d-1
    pre = (pre * seconds_per_hour).resample(time='1D').sum() / water_density * mm_per_m # ([kg m-2 s-1] * [s hr-1]).resample(time='1D').sum() / [kg m-3] * [mm m-1] = [mm d-1]
    pre = subset_dataset_to_max_full_years(pre, res='day', water_year=True)
    
    # 2025-01-26 
    # We need to disable these checks because we're now mainly using RDRS for forcing, and RDRS has a more limited time span than ERA5 does
    # As a result, if we match the hydrologic data to the RDRS period we'll lose a lot of data for stations with obs starting before 1980
    # We now thus subset the hydrologic data only inside the calculate_signatures() function, for the two signatures that use precip data
    # Keeping the code below for posterity
    #
    ## Match times between hydrologic data and precipitation
    #pre = pre.sel(time=slice(hyd['time'][0].values, hyd['time'][-1].values))
    #assert hyd['time'][0].values == pre['time'][0].values, 'attributes_from_streamflow(): mismatch between precipitation and streamflow start timestamp'
    #assert hyd['time'][-1].values == pre['time'][-1].values, 'attributes_from_streamflow(): mismatch between precipitation and streamflow final timestamp'
    #
    # Gap-fill the streamflow series if needed, so that we have identical length time series to work with
    hyd = hyd.resample(time='D').asfreq()
    #assert len(hyd['time']) == len(pre['time']), 'attributes_from_streamflow(): different number of timesteps in precipitation and streamflow series'

    # Create a water-year time variable
    hyd['water_year'] = hyd['time'].dt.year.where(hyd['time'].dt.month < 10, hyd['time'].dt.year + 1)
    pre['water_year'] = pre['time'].dt.year.where(pre['time'].dt.month < 10, pre['time'].dt.year + 1)
    
    # Track the data years used
    num_years = len(hyd.groupby('time.year'))
    l_values.append(num_years)
    l_index.append( ('Hydrology', 'num_years_hyd ', 'years', '-') )
    
    # Signatures
    l_values, l_index = calculate_signatures(hyd, pre, source, l_values, l_index)
    
    return l_values, l_index

def attributes_from_hydrolakes(geo_folder, dataset, basin_shp_path, ea_crs, l_values, l_index):
    
    '''Calculates open water attributes from HydroLAKES'''

    # We need some special actions if we're dealing with the distributed case, so check that first
    case = 'lumped'
    if 'distributed' in basin_shp_path:
        case = 'distributed'
    
    # Define the standard file name
    lake_str = str(geo_folder / dataset / 'raw' / 'HydroLAKES_polys_v10_NorthAmerica.shp')

    # Now handle the different cases
    if case == 'lumped':

        # Check if the file exists (some basins won't have lakes)
        if os.path.exists(lake_str):
            lakes = gpd.read_file(lake_str)
            res_mask  = lakes['Lake_type'] == 2 # Lake Type 2 == reservoir; see docs (https://data.hydrosheds.org/file/technical-documentation/HydroLAKES_TechDoc_v10.pdf)
            num_lakes = len(lakes)
            num_resvr = res_mask.sum()

        else: # no lakes shapefile
            lakes = None # this tells get_open_water_stats() that we had no lake shapefile
            num_lakes = 0
            num_resvr = 0

        l_values.append(num_lakes)
        l_index.append(('Open water', 'open_water_number',  '-', 'HydroLAKES'))
        l_values.append(num_resvr)
        l_index.append(('Open water', 'known_reservoirs',  '-', 'HydroLAKES'))

        # Summary stats
        l_values, l_index = get_open_water_stats(lakes, 'Lake_area', 'all', l_values, l_index) # All open water
        l_values, l_index = get_open_water_stats(lakes, 'Vol_total', 'all', l_values, l_index)
        l_values, l_index = get_open_water_stats(lakes, 'Lake_area', 'reservoir', l_values, l_index) # Reservoirs only
        l_values, l_index = get_open_water_stats(lakes, 'Vol_total', 'reservoir', l_values, l_index)

        # set the remaing output
        l_comids_lakes = None

    elif case == 'distributed':

        # Load the basin shape 
        basin = gpd.read_file(basin_shp_path)
        
        # Load the lake shape if we have one
        if os.path.exists(lake_str):
            lakes = gpd.read_file(lake_str)
        else:
            lakes = None # we use this to indicate we had no lakes at all
            
        # Loop over the individual polygons and create a new mini-HydroLAKES geodataframe for each polygon
        num_poly = len(basin)
        l_comids_lakes = basin['COMID'].values
        for i_poly in range(num_poly):

            # Rest the storage lists
            tmp_values = []
            tmp_index = []

            # subset the 'lakes' gdf to each individual polygon if we have 'lakes'           
            if lakes is not None:
                poly_lakes = subset_hydrolakes_to_subbasin(lakes, basin.iloc[i_poly]['geometry'], ea_crs)

                # Now see if we have a lake at all, and act accordingly
                if len(poly_lakes) > 0:
                    res_mask  = lakes['Lake_type'] == 2 # Lake Type 2 == reservoir; see docs (https://data.hydrosheds.org/file/technical-documentation/HydroLAKES_TechDoc_v10.pdf)
                    num_lakes = len(lakes)
                    num_resvr = res_mask.sum()
        
                else: # no lakes in this subbasins
                    poly_lakes = None # this tells get_open_water_stats() that we had no lake shapefile
                    num_lakes = 0
                    num_resvr = 0
            elif lakes is None: # this means we didn't have a lakes shapefile for this basin at all
                poly_lakes = None # this tells get_open_water_stats() we had nothing
                num_lakes = 0
                num_resvr = 0

            # Stats
            tmp_values.append(num_lakes)
            tmp_index.append(('Open water', 'open_water_number',  '-', 'HydroLAKES'))
            tmp_values.append(num_resvr)
            tmp_index.append(('Open water', 'known_reservoirs',  '-', 'HydroLAKES'))
            tmp_values, tmp_index = get_open_water_stats(poly_lakes, 'Lake_area', 'all', tmp_values, tmp_index) # All open water
            tmp_values, tmp_index = get_open_water_stats(poly_lakes, 'Vol_total', 'all', tmp_values, tmp_index)
            tmp_values, tmp_index = get_open_water_stats(poly_lakes, 'Lake_area', 'reservoir', tmp_values, tmp_index) # Reservoirs only
            tmp_values, tmp_index = get_open_water_stats(poly_lakes, 'Vol_total', 'reservoir', tmp_values, tmp_index)

            # Add the new values for this subbasin into the main l_values list
            l_values[i_poly] = tmp_values # This works because we create a unique dataframe just for HydroLAKES results

        # End of subbasin loop, add the new index entries only once
        l_index = tmp_index # we don't use append() in the distributed case because we'll be storing this in a dedicated HydroLAKES attribute df

    return l_values, l_index, l_comids_lakes

def attributes_from_merit(geo_folder, dataset, shp_str, riv_str, row, l_values, l_index, equal_area_crs='ESRI:102008', case='lumped'):
    
    '''Calculates topographic attributes from MERIT data'''

    # We need some special actions if we're dealing with the distributed case, so check that first
    case = 'lumped'
    if 'distributed' in shp_str:
        case = 'distributed'

    # Load the shapefile
    basin = gpd.read_file(shp_str)

    ## Known values
    if case == 'lumped':
        clon = basin.to_crs(equal_area_crs).centroid.to_crs('EPSG:4326').x.iloc[0] # lon
        clat = basin.to_crs(equal_area_crs).centroid.to_crs('EPSG:4326').y.iloc[0] # lat       
        area = row['Basin_area_km2']
        slat = row['Station_lat']
        slon = row['Station_lon']
        src  = row['Station_source']
        l_values.append(clat)
        l_index.append(('Topography', 'centroid_lat',  'degrees', 'MERIT Hydro'))
        l_values.append(clon)
        l_index.append(('Topography', 'centroid_lon',  'degrees', 'MERIT Hydro'))
        l_values.append(area)
        l_index.append(('Topography', 'basin_area', 'km^2', 'MERIT Hydro'))
        l_values.append(slat)
        l_index.append(('Topography', 'gauge_lat',  'degrees', 'USGS/WSC'))
        l_values.append(slon)
        l_index.append(('Topography', 'gauge_lon',  'degrees', 'USGS/WSC'))
    
    elif case == 'distributed':
        num_poly = len(basin)
        for i_poly in range(num_poly):
            clon = basin.to_crs(equal_area_crs).centroid.to_crs('EPSG:4326').x.iloc[i_poly] # lon
            clat = basin.to_crs(equal_area_crs).centroid.to_crs('EPSG:4326').y.iloc[i_poly] # lat       
            area = (basin.to_crs(equal_area_crs).area / 10**6).iloc[i_poly]
            l_values[i_poly].append(clat)
            l_values[i_poly].append(clon)
            l_values[i_poly].append(area)    
        l_index.append(('Topography', 'centroid_lat',  'degrees', 'MERIT Hydro')) # these are outside the subbasin loop because we need them only once
        l_index.append(('Topography', 'centroid_lon',  'degrees', 'MERIT Hydro'))
        l_index.append(('Topography', 'basin_area', 'km^2', 'MERIT Hydro'))    

    ## RASTERS
    # Slope and elevation can use zonal stats built-ins
    files = [str(geo_folder / dataset / 'raw'    / 'merit_hydro_elv.tif'),
             str(geo_folder / dataset / 'slope'  / 'merit_hydro_slope.tif')]
    attrs = ['elev',     'slope']
    units = ['m.a.s.l.', 'degrees']
    stats = ['mean', 'min', 'max', 'std']
    for tif,att,unit in zip(files,attrs,units):
        zonal_out = zonal_stats(shp_str, tif, stats=stats, all_touched=True)
        scale,offset = read_scale_and_offset(tif)
        l_values = update_values_list(l_values, stats, zonal_out, scale, offset, case=case)
        l_index += [('Topography', f'{att}_min',  f'{unit}', 'MERIT Hydro'),
                    ('Topography', f'{att}_mean', f'{unit}', 'MERIT Hydro'),
                    ('Topography', f'{att}_max',  f'{unit}', 'MERIT Hydro'),
                    ('Topography', f'{att}_std',  f'{unit}', 'MERIT Hydro')]

    # Aspect needs circular stats as custom input
    stats = ['min', 'max']
    tif = str(geo_folder / dataset / 'aspect' / 'merit_hydro_aspect.tif')
    zonal_out = zonal_stats(shp_str, tif, stats=stats, add_stats={'circ_mean': calc_circmean,
                                                                  'circ_std':  calc_circstd}, all_touched=True)
    scale,offset = read_scale_and_offset(tif)
    l_values = update_values_list(l_values, ['min', 'max', 'circ_mean', 'circ_std'], zonal_out, scale, offset, case=case)
    l_index += [('Topography', 'aspect_min',   'degrees', 'MERIT Hydro'),
                ('Topography', 'aspect_mean',  'degrees', 'MERIT Hydro'),
                ('Topography', 'aspect_max',   'degrees', 'MERIT Hydro'),
                ('Topography', 'aspect_std',   'degrees', 'MERIT Hydro')]
    #l_values, l_index = csa.get_aspect_attributes(tif,l_values,l_index)
    
    ## VECTOR
    l_values, l_index, merit_comids = get_river_attributes(riv_str, shp_str, l_values, l_index, area, equal_area_crs=equal_area_crs)
    
    return l_values, l_index, merit_comids

def attributes_from_lgrip30(geo_folder, dataset, shp_str, l_values, l_index, case='lumped'):

    '''Calculates percentage occurrence of all classes in LGRIP30 map'''

    tif = geo_folder / dataset / 'raw' / 'lgrip30_agriculture.tif'
    zonal_out = zonal_stats(shp_str, tif, categorical=True, all_touched=True)
    check_scale_and_offset(tif)
    l_values,l_index = update_values_list_with_categorical(l_values, l_index, zonal_out, 'LGRIP30', prefix='lc3_', case=case)
    return l_values, l_index

def attributes_from_modis_land(geo_folder, dataset, shp_str, l_values, l_index, case='lumped'):

    '''Calculates percentage occurrence of all classes in MODIS IGBP map'''

    tif = geo_folder / dataset / 'raw' / '2001_2022_mode_MCD12Q1_LC_Type1.tif'
    zonal_out = zonal_stats(shp_str, tif, categorical=True, all_touched=True)
    check_scale_and_offset(tif)
    l_values,l_index = update_values_list_with_categorical(l_values, l_index, zonal_out, 'MCD12Q1.061', prefix='lc2_', case=case)
    return l_values, l_index

def attributes_from_glclu2019(geo_folder, dataset, shp_str, l_values, l_index, case='lumped'):

    '''Calculates percentage occurrence of all classes in GLCLU2019 map'''

    tif = geo_folder / dataset / 'raw' / 'glclu2019_map.tif'
    zonal_out = zonal_stats(shp_str, tif, categorical=True, all_touched=True)
    check_scale_and_offset(tif)
    l_values,l_index = update_values_list_with_categorical(l_values, l_index, zonal_out, 'GLCLU 2019', prefix='lc1_', case=case)
    return l_values, l_index


def attributes_from_rdrs(met_folder, shp_path, dataset, l_values, l_index, use_mfdataset=False):

    '''Calculates a variety of metrics from RDRS data'''

    # Define file locations, depending on if we are dealing with lumped or distributed cases
    if 'lumped' in shp_path:
        rdrs_folder = met_folder / 'lumped'
        case = 'lumped'
    elif 'distributed' in shp_path:
        rdrs_folder = met_folder / 'distributed'
        case = 'distributed'
    rdrs_files = sorted( glob.glob( str(rdrs_folder / 'RDRS_*.nc') ) )
    print(f'Running {case} case for RDRS data.')

    # Open the data
    if use_mfdataset or (case == 'distributed'):
        ds = xr.open_mfdataset(rdrs_files, combine="by_coords") # Don't use 'decode_cf=False' > this somehow loses most of the timesteps
    else:
        ds = xr.merge([xr.open_dataset(f) for f in rdrs_files])

    # --- Act according to case
    if case == 'lumped':

        # -- Load the data into memory
        ds = ds.load()
        ds = ds.isel(hru=0) # We need this so the lumped and distributed cases have the same dimensions inside the calculation function: (time,nbnds)
                            # Without this, the lumped case has an extra 'hru' dimension (length 1) and this complicates value extraction
        
        # -- Get the precipitation for hydrologic signature calculations later; this avoids having to reload the entire dataset another time
        ds_precip = ds['RDRS_v2.1_A_PR0_SFC'].copy()

        # --- Calculate the statistics
        l_values, l_index = calculate_rdrs_stats_from_ds(ds,l_values,l_index)
        return l_values, l_index, ds_precip, ds
    
    elif case == 'distributed':

        comid_order = [] # we need these later, to ensure we line the forcing attributes up with the geospatial attributes correctly

        # Loop over the HRUs
        num_hru = len(ds['hru'])
        for i in range(0,num_hru):

            print(f'Running subbasin {i}')

            # Load the data for this sub-basin
            ds_hru = ds.isel(hru=i) # example
            ds_hru.load() # Load data into memory; done in place

            # Specifically track the hruID (COMID) so we can ensure correct matches with the rest of the attributes later
            assert (ds_hru['hruId'].values == ds_hru['hruId'][0].values).all(), f"COMIDs not all identical {ds_hru['hruId'].values}"
            comid_order.append(ds_hru['hruId'][0].values)

            # calculate the stats
            tmp_values = [] # statistics for this subbasin - we need a single empty list so we can properly store in the nested list later
            tmp_index  = [] # we also store the index in a tmp variable, so we only get this once instead of appending num_hru times to the main index list
            tmp_values, tmp_index = calculate_rdrs_stats_from_ds(ds_hru,tmp_values,tmp_index) 
            l_values.append(tmp_values)

        # prep outputs
        comid_order = [arr.tolist() for arr in comid_order] # convert list of arrays to simple list
        l_index.append(tmp_index) # need this only once; tmp_index is already a list so this creates a nested list: [['Climate', 'num_years_rdrs', 'years', 'RDRS')]]
        return l_values, l_index[0], comid_order, ds # we need l_index[0] to return a non-nested list

def calculate_rdrs_stats_from_ds(ds,l_values,l_index):

    # Define various conversion constants
    water_density = 1000 # kg m-3
    mm_per_m = 1000 # mm m-1
    seconds_per_hour = 60*60 # s h-1
    seconds_per_day = seconds_per_hour*24 # s d-1
    days_per_month = np.array([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]).reshape(-1, 1).flatten() # d month-1 || NOTE: added the .flatten() for distributed case, not tested for lumped
    days_per_year = days_per_month.sum()
    flip_sign = -1 # -; used to convert PET from negative (by convention this indicates an upward flux) to positive (ERA5 only)
    kelvin_to_celsius = -273.15
    pa_per_kpa = 1000 # Pa kPa-1
    
    # Select whole years only
    #   This avoids issues in cases where we have incomplete whole data years
    #   (e.g. 2000-06-01 to 2007-12-31) in basins with very seasonal weather
    #   (e.g. all precip occurs in Jan, Feb, Mar). By using only full years
    #   we avoid accidentally biasing the attributes.
    ds = subset_dataset_to_max_full_years(ds)
    num_years = len(ds.groupby('time.year'))
    l_values.append(num_years)
    l_index.append( ('Climate', 'num_years_rdrs', 'years', 'RDRS') )

    # --- Annual statistics (P, PET, T, aridity, seasonality, temperature, snow)
    # P
    yearly_pr0 = ds['RDRS_v2.1_A_PR0_SFC'].resample(time='1Y').mean() * seconds_per_day * days_per_year * mm_per_m / water_density # kg m-2 s-1 > mm yr-1
    l_values.append(yearly_pr0.mean().values)
    l_index.append(('Climate', 'PR0_mean', 'mm', 'RDRS'))
    l_values.append(yearly_pr0.std().values)
    l_index.append(('Climate', 'PR0_std', 'mm', 'RDRS'))

    # PET
    yearly_pet = ds['pet'].resample(time='1Y').mean() * seconds_per_day * days_per_year * mm_per_m / water_density # kg m-2 s-1 > mm yr-1
    l_values.append(yearly_pet.mean().values)
    l_index.append(('Climate', 'pet1_mean', 'mm', 'RDRS'))
    l_values.append(yearly_pet.std().values)
    l_index.append(('Climate', 'pet1_std', 'mm', 'RDRS'))

    # T
    yearly_tt = ds['RDRS_v2.1_P_TT_1.5m'].resample(time='1Y').mean() + kelvin_to_celsius # K > C
    l_values.append(yearly_tt.mean().values)
    l_index.append(('Climate', 'TT_mean', 'C', 'RDRS'))
    l_values.append(yearly_tt.std().values)
    l_index.append(('Climate', 'TT_std', 'C', 'RDRS'))

    # Aridity
    yearly_ari  = yearly_pet / yearly_pr0
    l_values.append(yearly_ari.mean().values)
    l_index.append(('Climate', 'aridity1_mean', '-', 'RDRS'))
    l_values.append(yearly_ari.std().values)
    l_index.append(('Climate', 'aridity1_std', '-', 'RDRS'))

    # Snow
    ds['snow'] = xr.where(ds['RDRS_v2.1_P_TT_1.5m'] < 273.15, ds['RDRS_v2.1_A_PR0_SFC'],0)
    yearly_snow = ds['snow'].resample(time='1Y').mean() * seconds_per_day * days_per_year * mm_per_m / water_density
    yearly_fs = yearly_snow / yearly_pr0
    l_values.append(yearly_fs.mean().values)
    l_index.append(('Climate', 'fracsnow1_mean', '-', 'RDRS'))
    l_values.append(yearly_fs.std().values)
    l_index.append(('Climate', 'fracsnow1_std', '-', 'RDRS'))

    # Seasonality
    seasonality = find_climate_seasonality_rdrs(ds,use_typical_cycle=False)
    l_values.append(seasonality.mean())
    l_index.append(('Climate', 'seasonality1_mean', '-', 'RDRS'))
    l_values.append(seasonality.std())
    l_index.append(('Climate', 'seasonality1_std', '-', 'RDRS'))

    # --- Monthly attributes
    # Calculate monthly PET in mm
    #      kg m-2 s-1 / kg m-3
    # mm month-1 = kg m-2 s-1 * kg-1 m3 * s d-1 * d month-1 * mm m-1 * -
    monthly_pet = ds['pet'].resample(time='1M').mean().groupby('time.month')
    pet_m = monthly_pet.mean() / water_density * seconds_per_day * days_per_month * mm_per_m  # [kg m-2 s-1] to [mm month-1]; negative to indicate upward flux
    pet_s = monthly_pet.std() / water_density * seconds_per_day * days_per_month * mm_per_m  # [kg m-2 s-1] to [mm month-1]; negative to indicate upward flux
    l_values, l_index = process_monthly_means_to_lists(pet_m, 'mean', l_values, l_index, 'pet1', 'mm', source='RDRS')
    l_values, l_index = process_monthly_means_to_lists(pet_s, 'std', l_values, l_index, 'pet1', 'mm', source='RDRS')

    # Same for precipitation: [mm month-1]
    monthly_pr0 = ds['RDRS_v2.1_A_PR0_SFC'].resample(time='1M').mean().groupby('time.month')
    pr0_m = monthly_pr0.mean() / water_density * seconds_per_day * days_per_month * mm_per_m # [kg m-2 s-1] to [mm month-1]
    pr0_s = monthly_pr0.std() / water_density * seconds_per_day * days_per_month * mm_per_m # [kg m-2 s-1] to [mm month-1]
    l_values, l_index = process_monthly_means_to_lists(pr0_m, 'mean', l_values, l_index, 'PR0', 'mm', source='RDRS')
    l_values, l_index = process_monthly_means_to_lists(pr0_s, 'std', l_values, l_index, 'PR0', 'mm', source='RDRS')

    # Monthly temperature statistics [C]
    monthly_tavg = (ds['RDRS_v2.1_P_TT_1.5m'].resample(time='1D').mean().resample(time='1M').mean() + kelvin_to_celsius).groupby('time.month')
    tavg_m = monthly_tavg.mean()
    tavg_s = monthly_tavg.std()
    l_values, l_index = process_monthly_means_to_lists(tavg_m, 'mean', l_values, l_index, 'tdavg', 'C', source = 'RDRS')
    l_values, l_index = process_monthly_means_to_lists(tavg_s, 'std', l_values, l_index, 'tdavg', 'C', source = 'RDRS')

    monthly_tmin = (ds['RDRS_v2.1_P_TT_1.5m'].resample(time='1D').min().resample(time='1M').mean() + kelvin_to_celsius).groupby('time.month')
    tmin_m = monthly_tmin.mean()
    tmin_s = monthly_tmin.std()
    l_values, l_index = process_monthly_means_to_lists(tmin_m, 'mean', l_values, l_index, 'tdmin', 'C', source = 'RDRS')
    l_values, l_index = process_monthly_means_to_lists(tmin_m, 'std', l_values, l_index, 'tdmin', 'C', source = 'RDRS')
    
    monthly_tmax = (ds['RDRS_v2.1_P_TT_1.5m'].resample(time='1D').max().resample(time='1M').mean() + kelvin_to_celsius).groupby('time.month')
    tmax_m = monthly_tmax.mean()
    tmax_s = monthly_tmax.std()
    l_values, l_index = process_monthly_means_to_lists(tmax_m, 'mean', l_values, l_index, 'tdmax', 'C', source = 'RDRS')
    l_values, l_index = process_monthly_means_to_lists(tmax_s, 'std', l_values, l_index, 'tdmax', 'C', source = 'RDRS')

    # Monthly shortwave and longwave [W m-2]
    monthly_sw = ds['RDRS_v2.1_P_FB_SFC'].resample(time='1M').mean().groupby('time.month')
    sw_m = monthly_sw.mean()
    sw_s = monthly_sw.std()
    l_values, l_index = process_monthly_means_to_lists(sw_m, 'mean', l_values, l_index, 'FB', 'W m^-2', source = 'RDRS')
    l_values, l_index = process_monthly_means_to_lists(sw_s, 'std', l_values, l_index, 'FB', 'W m^-2', source = 'RDRS')
    
    monthly_lw = ds['RDRS_v2.1_P_FI_SFC'].resample(time='1M').mean().groupby('time.month')
    lw_m = monthly_lw.mean(dim='time')
    lw_s = monthly_lw.std(dim='time')
    l_values, l_index = process_monthly_means_to_lists(lw_m, 'mean', l_values, l_index, 'FI', 'W m^-2', source = 'RDRS')
    l_values, l_index = process_monthly_means_to_lists(lw_s, 'std', l_values, l_index, 'FI', 'W m^-2', source = 'RDRS')

    # Surface pressure [Pa]
    monthly_sp = ds['RDRS_v2.1_P_P0_SFC'].resample(time='1M').mean().groupby('time.month')
    sp_m = monthly_sp.mean() / pa_per_kpa # [Pa] > [kPa]
    sp_s = monthly_sp.std() / pa_per_kpa
    l_values, l_index = process_monthly_means_to_lists(sp_m, 'mean', l_values, l_index, 'P0', 'kPa', source = 'RDRS')
    l_values, l_index = process_monthly_means_to_lists(sp_s, 'std', l_values, l_index, 'P0', 'kPa', source = 'RDRS')

    # Humidity [-]
    monthly_q = ds['RDRS_v2.1_P_HU_1.5m'].resample(time='1M').mean().groupby('time.month') # specific
    q_m = monthly_q.mean()
    q_s = monthly_q.std()
    l_values, l_index = process_monthly_means_to_lists(q_m, 'mean', l_values, l_index, 'HU', 'kg kg^-1', source = 'RDRS')
    l_values, l_index = process_monthly_means_to_lists(q_s, 'std', l_values, l_index, 'HU', 'kg kg^-1', source = 'RDRS')
    
    monthly_rh = ds['RDRS_v2.1_P_HR_1.5m'].resample(time='1M').mean().groupby('time.month') # relative
    rh_m = monthly_rh.mean()
    rh_s = monthly_rh.std()
    l_values, l_index = process_monthly_means_to_lists(rh_m, 'mean', l_values, l_index, 'HR', 'kPa kPa^-1', source = 'RDRS')
    l_values, l_index = process_monthly_means_to_lists(rh_s, 'std', l_values, l_index, 'HR', 'kPa kPa^-1', source = 'RDRS')

    # Wind speed [m s-1]
    monthly_w = ds['RDRS_v2.1_P_UVC_10m'].resample(time='1M').mean().groupby('time.month')
    w_m = monthly_w.mean()
    w_s = monthly_w.std()
    l_values, l_index = process_monthly_means_to_lists(w_m, 'mean', l_values, l_index, 'UVC', 'm s^-1', source = 'RDRS')
    l_values, l_index = process_monthly_means_to_lists(w_s, 'std', l_values, l_index, 'UVC', 'm s^-1', source = 'RDRS')

    # Wind direction
    monthly_phi = ds['phi'].resample(time='1M').apply(circmean_group).groupby('time.month')
    phi_m = monthly_phi.apply(circmean_group)
    phi_s = monthly_phi.apply(circstd_group)
    l_values, l_index = process_monthly_means_to_lists(phi_m, 'mean', l_values, l_index, 'phi', 'degrees', source = 'RDRS')
    l_values, l_index = process_monthly_means_to_lists(phi_s, 'std', l_values, l_index, 'phi', 'degrees', source = 'RDRS')

    # aridity
    monthly_pet = ds['pet'].resample(time='1M').mean()
    monthly_pr0 = ds['RDRS_v2.1_A_PR0_SFC'].resample(time='1M').mean()
    if (monthly_pr0 == 0).any():
        print(f'--- WARNING: attributes_from_rdrs(): adding 1 mm to monthly precipitation to avoid divide by zero error in aridity calculation')
        monthly_pr0[(monthly_pr0 == 0)] = 1 / mm_per_m * water_density / (seconds_per_day * days_per_month.mean()) # [mm month-1] / [mm m-1] * [kg m-3] / ([s d-1] * [d month-1]) = [kg m-2 s-1]
    monthly_ari = (monthly_pet / monthly_pr0).groupby('time.month')
    ari_m = monthly_ari.mean()
    ari_s = monthly_ari.std()
    l_values, l_index = process_monthly_means_to_lists(ari_m, 'mean', l_values, l_index, 'aridity1', '-', source = 'RDRS')
    l_values, l_index = process_monthly_means_to_lists(ari_s, 'std', l_values, l_index, 'aridity1', '-', source = 'RDRS')

    # snow
    monthly_snow = ds['snow'].resample(time='1M').mean()
    monthly_pr0 = ds['RDRS_v2.1_A_PR0_SFC'].resample(time='1M').mean()
    if (monthly_pr0 == 0).any():
        print(f'--- WARNING: attributes_from_rdrs(): adding 1 mm to monthly precipitation to avoid divide by zero error in snow calculation. Note that by definition this cannot change the fraction snow result (if there is 0 precip, none of it will fall as snow)')
        monthly_pr0[(monthly_pr0 == 0)] = 1 / mm_per_m * water_density / (seconds_per_day * days_per_month.mean()) # [mm month-1] / [mm m-1] * [kg m-3] / ([s d-1] * [d month-1]) = [kg m-2 s-1]
    monthly_snow = (monthly_snow / monthly_pr0).groupby('time.month')
    fsnow_m = monthly_snow.mean()
    fsnow_s = monthly_snow.std()
    l_values, l_index = process_monthly_means_to_lists(fsnow_m, 'mean', l_values, l_index, 'fracsnow1', '-', source = 'RDRS')
    l_values, l_index = process_monthly_means_to_lists(fsnow_s, 'std', l_values, l_index, 'fracsnow1', '-', source = 'RDRS') 

    # --- High-frequency statistics (high/low duration/timing/magnitude)
    #  Everyone does precip. We'll add temperature too as a drought/frost indicator
    
    # -- LOW TEMPERATURE
    variable  = 'RDRS_v2.1_P_TT_1.5m'
    low_threshold = 273.15 # K, freezing point
    low_condition = ds[variable] < low_threshold
    l_values,l_index = calculate_temp_prcp_stats('temp',low_condition,'low',l_values,l_index, dataset='RDRS')

    # -- HIGH TEMPERATURE
    # WMO defines a heat wave as a 5-day or longer period with maximum daily temperatures 5C above 
    # "standard" daily max temperature (1961-1990; source:
    # https://www.ifrc.org/sites/default/files/2021-06/10-HEAT-WAVE-HR.pdf).
    # We define a "hot day" therefore as a day with a maximum temperature 5 degrees over the 
    # the long-term mean maximum temperature.
    #   Note: we don't have 1961-1990 data for some stations, so we stick with long-term mean.
    #   Note: this will in most cases slightly underestimate heat waves compared to WMO definition
    
    # First, we identify the long-term mean daily maximum temperature in a dedicated function
    var = 'RDRS_v2.1_P_TT_1.5m'
    high_threshold = create_mean_daily_max_series(ds,var=var)
    
    # Next, we check if which 't' values are 5 degrees above the long-term mean daily max 
    #  ("(ds['t'] > result_array + 5)"), and resample this to a daily time series 
    #  ("resample(time='1D')") filled with "True" if any value in that day was True.
    daily_flags = (ds[var] > high_threshold + 5).resample(time='1D').any()
    
    # Finally, we reindex these daily flags back onto the hourly time series by filling values
    high_condition = daily_flags.reindex_like(ds[var], method='ffill')
    
    # Now calculate stats like before
    l_values,l_index = calculate_temp_prcp_stats('temp',high_condition,'high',l_values,l_index, dataset='RDRS')

    # -- LOW PRECIPITATION
    variable = 'RDRS_v2.1_A_PR0_SFC'
    # We'll stick with the original CAMELS definition of low precipitation: < 1 mm day-1
    # It may not make too much sense to look at "dry hours" so we'll do this analysis at daily step
    low_threshold = 1 # [mm d-1]
    # Create daily precipitation sum (divided by density, times mm m-1 cancels out)
    # [kg m-2 s-1] * [s h-1] / [kg m-3] * [mm m-1] = [mm h-1]
    low_condition = (ds[variable] * seconds_per_hour).resample(time='1D').sum() < low_threshold
    l_values,l_index = calculate_temp_prcp_stats('prec',low_condition,'low',l_values,l_index,
                                             units='days', dataset='RDRS') # this 'units' argument prevents conversion to days inside the functiom
    
    # -- HIGH PRECIPITATION
    # CAMELS: > 5 times mean daily precip
    high_threshold = 5 * (ds[variable] * seconds_per_hour).resample(time='1D').sum().mean()
    high_condition = (ds[variable] * seconds_per_hour).resample(time='1D').sum() >= high_threshold
    l_values,l_index = calculate_temp_prcp_stats('prec',high_condition,'high',l_values,l_index,
                                                 units='days', dataset='RDRS')

    return l_values, l_index


def attributes_from_era5(met_folder, shp_path, dataset, l_values, l_index, use_mfdataset=False):

    '''Calculates a variety of metrics from ERA5 data'''

    # Define various conversion constants
    water_density = 1000 # kg m-3
    mm_per_m = 1000 # mm m-1
    seconds_per_hour = 60*60 # s h-1
    seconds_per_day = seconds_per_hour*24 # s d-1
    days_per_month = np.array([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]).reshape(-1, 1) # d month-1
    days_per_year = days_per_month.sum()
    flip_sign = -1 # -; used to convert PET from negative (by convention this indicates an upward flux) to positive
    kelvin_to_celsius = -273.15
    pa_per_kpa = 1000 # Pa kPa-1

    # Define file locations, depending on if we are dealing with lumped or distributed cases
    if 'lumped' in shp_path:
        era_folder = met_folder / 'lumped'
    elif 'distributed' in shp_path:
        era_folder = met_folder / 'distributed'
    era_files = sorted( glob.glob( str(era_folder / 'ERA5_*.nc') ) )

    # Open the data
    if use_mfdataset:
        ds = xr.open_mfdataset(era_files, engine='netcdf4')
        ds = ds.load() # load the whole thing into memory instead of lazy-loading
    else:
        ds = xr.merge([xr.open_dataset(f) for f in era_files])
        ds = ds.load()
    
    # -- Get the precipitation for hydrologic signature calculations
    # This avoids having to reload the entire dataset another time
    ds_precip = ds['mtpr'].copy()    
    
    # Select whole years only
    #   This avoids issues in cases where we have incomplete whole data years
    #   (e.g. 2000-06-01 to 2007-12-31) in basins with very seasonal weather
    #   (e.g. all precip occurs in Jan, Feb, Mar). By using only full years
    #   we avoid accidentally biasing the attributes.
    ds = subset_dataset_to_max_full_years(ds)
    num_years = len(ds.groupby('time.year'))
    l_values.append(num_years)
    l_index.append( ('Climate', 'num_years_era5', 'years', 'ERA5') )
    
    # --- Annual statistics (P, PET, T, aridity, seasonality, temperature, snow)
    # P
    yearly_mtpr = ds['mtpr'].resample(time='1YE').mean() * seconds_per_day * days_per_year * mm_per_m / water_density # kg m-2 s-1 > mm yr-1
    l_values.append(yearly_mtpr.mean().values)
    l_index.append(('Climate', 'mtpr_mean', 'mm', 'ERA5'))
    l_values.append(yearly_mtpr.std().values)
    l_index.append(('Climate', 'mtpr_std', 'mm', 'ERA5'))
    
    # PET
    yearly_mper = ds['mper'].resample(time='1YE').mean() * flip_sign * seconds_per_day * days_per_year * mm_per_m / water_density # kg m-2 s-1 > mm yr-1
    l_values.append(yearly_mper.mean().values)
    l_index.append(('Climate', 'mper_mean', 'mm', 'ERA5'))
    l_values.append(yearly_mper.std().values)
    l_index.append(('Climate', 'mper_std', 'mm', 'ERA5'))
    
    # T
    yearly_tavg = ds['t'].resample(time='1YE').mean() + kelvin_to_celsius # K > C
    l_values.append(yearly_tavg.mean().values)
    l_index.append(('Climate', 'tdavg_mean', 'C', 'ERA5'))
    l_values.append(yearly_tavg.std().values)
    l_index.append(('Climate', 'tdavg_std', 'C', 'ERA5'))
    
    # Aridity
    yearly_ari  = yearly_mper / yearly_mtpr
    l_values.append(yearly_ari.mean().values)
    l_index.append(('Climate', 'aridity1_mean', '-', 'ERA5'))
    l_values.append(yearly_ari.std().values)
    l_index.append(('Climate', 'aridity1_std', '-', 'ERA5'))

    # Snow
    ds['snow'] = xr.where(ds['t'] < 273.15, ds['mtpr'],0) # rates: kg m-2 s-1
    yearly_snow = ds['snow'].resample(time='1YE').mean() * seconds_per_day * days_per_year * mm_per_m / water_density # kg m-2 s-1 > mm yr-1
    yearly_fs   = yearly_snow / yearly_mtpr
    l_values.append(yearly_fs.mean().values)
    l_index.append(('Climate', 'fracsnow1_mean', '-', 'ERA5'))
    l_values.append(yearly_fs.std().values)
    l_index.append(('Climate', 'fracsnow1_std', '-', 'ERA5'))

    # Seasonality
    seasonality = find_climate_seasonality_era5(ds,use_typical_cycle=False)
    l_values.append(seasonality.mean())
    l_index.append(('Climate', 'seasonality1_mean', '-', 'ERA5'))
    l_values.append(seasonality.std())
    l_index.append(('Climate', 'seasonality1_std', '-', 'ERA5'))

    # --- Monthly attributes
    # Calculate monthly PET in mm
    #      kg m-2 s-1 / kg m-3
    # mm month-1 = kg m-2 s-1 * kg-1 m3 * s d-1 * d month-1 * mm m-1 * -
    monthly_mper = ds['mper'].resample(time='1ME').mean().groupby('time.month') 
    mper_m = monthly_mper.mean() / water_density * seconds_per_day * days_per_month * mm_per_m * flip_sign  # [kg m-2 s-1] to [mm month-1]; negative to indicate upward flux
    mper_s = monthly_mper.std() / water_density * seconds_per_day * days_per_month * mm_per_m  # [kg m-2 s-1] to [mm month-1]; negative to indicate upward flux
    l_values, l_index = process_monthly_means_to_lists(mper_m, 'mean', l_values, l_index, 'mper', 'mm')
    l_values, l_index = process_monthly_means_to_lists(mper_s, 'std', l_values, l_index, 'mper', 'mm')
        
    # Same for precipitation: [mm month-1]
    monthly_mtpr = ds['mtpr'].resample(time='1ME').mean().groupby('time.month')
    mtpr_m = monthly_mtpr.mean() / water_density * seconds_per_day * days_per_month * mm_per_m # [kg m-2 s-1] to [mm month-1]
    mtpr_s = monthly_mtpr.std() / water_density * seconds_per_day * days_per_month * mm_per_m # [kg m-2 s-1] to [mm month-1]
    l_values, l_index = process_monthly_means_to_lists(mtpr_m, 'mean', l_values, l_index, 'mtpr', 'mm')
    l_values, l_index = process_monthly_means_to_lists(mtpr_s, 'std', l_values, l_index, 'mtpr', 'mm')
    
    # Monthly temperature statistics [C]
    monthly_tavg = (ds['t'].resample(time='1D').mean().resample(time='1ME').mean() + kelvin_to_celsius).groupby('time.month')
    tavg_m = monthly_tavg.mean()
    tavg_s = monthly_tavg.std()
    l_values, l_index = process_monthly_means_to_lists(tavg_m, 'mean', l_values, l_index, 'tdavg', 'C')
    l_values, l_index = process_monthly_means_to_lists(tavg_s, 'std', l_values, l_index, 'tdavg', 'C')
    
    monthly_tmin = (ds['t'].resample(time='1D').min().resample(time='1ME').mean() + kelvin_to_celsius).groupby('time.month')
    tmin_m = monthly_tmin.mean()
    tmin_s = monthly_tmin.std()
    l_values, l_index = process_monthly_means_to_lists(tmin_m, 'mean', l_values, l_index, 'tdmin', 'C')
    l_values, l_index = process_monthly_means_to_lists(tmin_m, 'std', l_values, l_index, 'tdmin', 'C')
    
    monthly_tmax = (ds['t'].resample(time='1D').max().resample(time='1ME').mean() + kelvin_to_celsius).groupby('time.month')
    tmax_m = monthly_tmax.mean()
    tmax_s = monthly_tmax.std()
    l_values, l_index = process_monthly_means_to_lists(tmax_m, 'mean', l_values, l_index, 'tdmax', 'C')
    l_values, l_index = process_monthly_means_to_lists(tmax_s, 'std', l_values, l_index, 'tdmax', 'C')
    
    # Monthly shortwave and longwave [W m-2]
    monthly_sw = ds['msdwswrf'].resample(time='1ME').mean().groupby('time.month')
    sw_m = monthly_sw.mean()
    sw_s = monthly_sw.std()
    l_values, l_index = process_monthly_means_to_lists(sw_m, 'mean', l_values, l_index, 'msdwswrf', 'W m^-2')
    l_values, l_index = process_monthly_means_to_lists(sw_s, 'std', l_values, l_index, 'msdwswrf', 'W m^-2')
    
    monthly_lw = ds['msdwlwrf'].resample(time='1ME').mean().groupby('time.month')
    lw_m = monthly_lw.mean(dim='time')
    lw_s = monthly_lw.std(dim='time')
    l_values, l_index = process_monthly_means_to_lists(lw_m, 'mean', l_values, l_index, 'msdwlwrf', 'W m^-2')
    l_values, l_index = process_monthly_means_to_lists(lw_s, 'std', l_values, l_index, 'msdwlwrf', 'W m^-2')

    # Surface pressure [Pa]
    monthly_sp = ds['sp'].resample(time='1ME').mean().groupby('time.month')
    sp_m = monthly_sp.mean() / pa_per_kpa # [Pa] > [kPa]
    sp_s = monthly_sp.std() / pa_per_kpa
    l_values, l_index = process_monthly_means_to_lists(sp_m, 'mean', l_values, l_index, 'sp', 'kPa')
    l_values, l_index = process_monthly_means_to_lists(sp_s, 'std', l_values, l_index, 'sp', 'kPa')
    
    # Humidity [-]
    monthly_q = ds['q'].resample(time='1ME').mean().groupby('time.month') # specific
    q_m = monthly_q.mean()
    q_s = monthly_q.std()
    l_values, l_index = process_monthly_means_to_lists(q_m, 'mean', l_values, l_index, 'q', 'kg kg^-1')
    l_values, l_index = process_monthly_means_to_lists(q_s, 'std', l_values, l_index, 'q', 'kg kg^-1')
    
    monthly_rh = ds['rh'].resample(time='1ME').mean().groupby('time.month') # relative
    rh_m = monthly_rh.mean()
    rh_s = monthly_rh.std()
    l_values, l_index = process_monthly_means_to_lists(rh_m, 'mean', l_values, l_index, 'rh', 'kPa kPa^-1')
    l_values, l_index = process_monthly_means_to_lists(rh_s, 'std', l_values, l_index, 'rh', 'kPa kPa^-1')
    
    # Wind speed [m s-1]
    monthly_w = ds['w'].resample(time='1ME').mean().groupby('time.month')
    w_m = monthly_w.mean()
    w_s = monthly_w.std()
    l_values, l_index = process_monthly_means_to_lists(w_m, 'mean', l_values, l_index, 'w', 'm s^-1')
    l_values, l_index = process_monthly_means_to_lists(w_s, 'std', l_values, l_index, 'w', 'm s^-1')
    
    # Wind direction
    monthly_phi = ds['phi'].resample(time='1ME').apply(circmean_group).groupby('time.month')
    phi_m = monthly_phi.apply(circmean_group)
    phi_s = monthly_phi.apply(circstd_group)
    l_values, l_index = process_monthly_means_to_lists(phi_m, 'mean', l_values, l_index, 'phi', 'degrees')
    l_values, l_index = process_monthly_means_to_lists(phi_s, 'std', l_values, l_index, 'phi', 'degrees')
    
    # --- Long-term monthly statistics (aridity, seasonality, snow)
    # aridity
    monthly_mper = ds['mper'].resample(time='1ME').mean() * flip_sign
    monthly_mtpr = ds['mtpr'].resample(time='1ME').mean()
    if (monthly_mtpr == 0).any():
        print(f'--- WARNING: attributes_from_era5(): adding 1 mm to monthly precipitation to avoid divide by zero error in aridity calculation')
        monthly_mtpr[(monthly_mtpr == 0)] = 1 / mm_per_m * water_density / (seconds_per_day * days_per_month.mean()) # [mm month-1] / [mm m-1] * [kg m-3] / ([s d-1] * [d month-1]) = [kg m-2 s-1]
    monthly_ari = (monthly_mper / monthly_mtpr).groupby('time.month')
    ari_m = monthly_ari.mean()
    ari_s = monthly_ari.std()
    l_values, l_index = process_monthly_means_to_lists(ari_m, 'mean', l_values, l_index, 'aridity1', '-')
    l_values, l_index = process_monthly_means_to_lists(ari_s, 'std', l_values, l_index, 'aridity1', '-')

    # snow
    monthly_snow = ds['snow'].resample(time='1ME').mean()
    monthly_mtpr = ds['mtpr'].resample(time='1ME').mean()
    if (monthly_mtpr == 0).any():
        print(f'--- WARNING: attributes_from_era5(): adding 1 mm to monthly precipitation to avoid divide by zero error in snow calculation. Note that by definition this cannot change the fraction snow result (if there is 0 precip, none of it will fall as snow)')
        monthly_mtpr[(monthly_mtpr == 0)] = 1 / mm_per_m * water_density / (seconds_per_day * days_per_month.mean()) # [mm month-1] / [mm m-1] * [kg m-3] / ([s d-1] * [d month-1]) = [kg m-2 s-1]
    monthly_snow = (monthly_snow / monthly_mtpr).groupby('time.month')
    fsnow_m = monthly_snow.mean()
    fsnow_s = monthly_snow.std()
    l_values, l_index = process_monthly_means_to_lists(fsnow_m, 'mean', l_values, l_index, 'fracsnow1', '-')
    l_values, l_index = process_monthly_means_to_lists(fsnow_s, 'std', l_values, l_index, 'fracsnow1', '-')  
    
    # --- High-frequency statistics (high/low duration/timing/magnitude)
    #  Everyone does precip. We'll add temperature too as a drought/frost indicator
    #  ERA5 only
    
    # -- LOW TEMPERATURE
    variable  = 't'
    low_threshold = 273.15 # K, freezing point
    low_condition = ds[variable] < low_threshold
    l_values,l_index = calculate_temp_prcp_stats('temp',low_condition,'low',l_values,l_index)
    
    # -- HIGH TEMPERATURE
    # WMO defines a heat wave as a 5-day or longer period with maximum daily temperatures 5C above 
    # "standard" daily max temperature (1961-1990; source:
    # https://www.ifrc.org/sites/default/files/2021-06/10-HEAT-WAVE-HR.pdf).
    # We define a "hot day" therefore as a day with a maximum temperature 5 degrees over the 
    # the long-term mean maximum temperature.
    #   Note: we don't have 1961-1990 data for some stations, so we stick with long-term mean.
    #   Note: this will in most cases slightly underestimate heat waves compared to WMO definition
    
    # First, we identify the long-term mean daily maximum temperature in a dedicated function
    high_threshold = create_mean_daily_max_series(ds,var='t')
    
    # Next, we check if which 't' values are 5 degrees above the long-term mean daily max 
    #  ("(ds['t'] > result_array + 5)"), and resample this to a daily time series 
    #  ("resample(time='1D')") filled with "True" if any value in that day was True.
    daily_flags = (ds['t'] > high_threshold + 5).resample(time='1D').any()
    
    # Finally, we reindex these daily flags back onto the hourly time series by filling values
    high_condition = daily_flags.reindex_like(ds['t'], method='ffill')
    
    # Now calculate stats like before
    l_values,l_index = calculate_temp_prcp_stats('temp',high_condition,'high',l_values,l_index)
    
    # -- LOW PRECIPITATION
    variable = 'mtpr'
    # We'll stick with the original CAMELS definition of low precipitation: < 1 mm day-1
    # It may not make too much sense to look at "dry hours" so we'll do this analysis at daily step
    low_threshold = 1 # [mm d-1]
    # Create daily precipitation sum (divided by density, times mm m-1 cancels out)
    # [kg m-2 s-1] * [s h-1] / [kg m-3] * [mm m-1] = [mm h-1]
    low_condition = (ds[variable] * seconds_per_hour).resample(time='1D').sum() < low_threshold
    l_values,l_index = calculate_temp_prcp_stats('prec',low_condition,'low',l_values,l_index,
                                             units='days') # this 'units' argument prevents conversion to days inside the functiom
    
    # -- HIGH PRECIPITATION
    # CAMELS: > 5 times mean daily precip
    high_threshold = 5 * (ds[variable] * seconds_per_hour).resample(time='1D').sum().mean()
    high_condition = (ds[variable] * seconds_per_hour).resample(time='1D').sum() >= high_threshold
    l_values,l_index = calculate_temp_prcp_stats('prec',high_condition,'high',l_values,l_index,
                                                 units='days')

    return l_values, l_index, ds_precip, ds


def attributes_from_forest_height(geo_folder, dataset, shp_str, l_values, index, case='lumped'):

    '''Calculates mean, min, max and stdv for forest height 2000 and 2020 tifs'''

    # Year 2000 min, mean, max, stdev
    tif = str( geo_folder / dataset / 'raw' / 'forest_height_2000.tif' )
    stats = ['mean', 'min', 'max', 'std']
    zonal_out = zonal_stats(shp_str, tif, stats=stats, all_touched=True)
    scale,offset = read_scale_and_offset(tif)
    l_values = update_values_list(l_values, stats, zonal_out, scale, offset, case=case)
    index += [('Land cover', 'forest_height_2000_min',   'm', 'GLCLUC 2000-2020'),
              ('Land cover', 'forest_height_2000_mean',  'm', 'GLCLUC 2000-2020'),
              ('Land cover', 'forest_height_2000_max',   'm', 'GLCLUC 2000-2020'),
              ('Land cover', 'forest_height_2000_std',   'm', 'GLCLUC 2000-2020')]

    # Year 2020 mean, stdev
    tif = geo_folder / dataset / 'raw' / 'forest_height_2020.tif'
    stats = ['mean', 'min', 'max', 'std']
    zonal_out = zonal_stats(shp_str, tif, stats=stats, all_touched=True)
    scale,offset = read_scale_and_offset(tif)
    l_values = update_values_list(l_values, stats, zonal_out, scale, offset, case=case)
    index += [('Land cover', 'forest_height_2020_min',   'm', 'GLCLUC 2000-2020'),
              ('Land cover', 'forest_height_2020_mean',  'm', 'GLCLUC 2000-2020'),
              ('Land cover', 'forest_height_2020_max',   'm', 'GLCLUC 2000-2020'),
              ('Land cover', 'forest_height_2020_std',   'm', 'GLCLUC 2000-2020')]

    return l_values, index

def attributes_from_lai(geo_folder, dataset, temp_path, shp_str, l_values, index, case='lumped'):

    '''Calculates mean and stdv for tifs of monthly LAI values'''

    # Calculate monthly mean maps (e.g. mean Jan, Feb, etc.)
    lai_folder = geo_folder / dataset / 'raw' 
    lai_files = sorted( glob.glob(str(lai_folder / '*.tif')) ) # Find LAI files
    month_files = calculate_monthly_lai_maps(lai_files, temp_path) # Create 12 monthly maps

    # Monthly mean, stdev LAI; monthly mean, stdev GVF
    for month_file in month_files:
        stats = ['mean', 'std']
        zonal_out = zonal_stats(shp_str, month_file, stats=stats, all_touched=True)
        scale, offset = read_scale_and_offset(month_file)
        scale,offset = read_scale_and_offset(month_file)
        l_values = update_values_list(l_values, stats, zonal_out, scale, offset, case=case)
        month = os.path.basename(month_file).split('_')[2]
        index += [('Land cover', f'lai_mean_month_{month}',  'm^2 m^-2', 'MCD15A2H.061'),
                  ('Land cover', f'lai_std_month_{month}',   'm^2 m^-2', 'MCD15A2H.061')]
    
    # Clear temp folder
    files_to_remove = os.listdir(temp_path)
    for file_to_remove in files_to_remove:
        file_remove_path = os.path.join(temp_path, file_to_remove)
        if os.path.isfile(file_remove_path):
            os.remove(file_remove_path)
    
    return l_values, index

def attributes_from_worldclim(geo_folder, dataset, shp_str, l_values, l_index, case='lumped'):

    '''Calculates mean and stdv for tifs of monthly WorldClim values'''

    # Define file locations
    # Units source: https://www.worldclim.org/data/worldclim21.html
    clim_folder = geo_folder / dataset / 'raw'
    sub_folders =      ['prec', 'srad',   'tavg', 'tmax', 'tmin', 'vapr', 'wind',   'pet', 'aridity2', 'fracsnow2'] # aridity and fractionsnow have subscript 2 to distinguish them from ERA5 attributes
    sub_folder_units = ['mm',   'W m^-2', 'C',    'C',    'C',    'kPa',  'm s^-1', 'mm',  '-',       '-']  # srad original: kJ m^-2 d^-1. Converted below

    # Get the annual values
    l_values,l_index = get_annual_worldclim_attributes(clim_folder, shp_str, 'prec', 'tavg', 'pet', 'snow2', l_values, l_index)

    # Loop over the files and calculate the stats
    for sub_folder, sub_folder_unit in zip(sub_folders, sub_folder_units):
        month_files = sorted( glob.glob(str(clim_folder / sub_folder / '*.tif')) )
        for month_file in month_files:
            month_file = clim_folder / sub_folder / month_file # Construct the full path, because listdir() gives only files
            stats = ['mean', 'std']
            zonal_out = zonal_stats(shp_str, month_file, stats=stats, all_touched=True)
            
            scale, offset = read_scale_and_offset(month_file)
            if sub_folder == 'srad':
                zonal_out = zonal_stats_unit_conversion(zonal_out,stats,'srad', scale, offset)

            l_values = update_values_list(l_values, stats, zonal_out, scale, offset, case=case)
            
            month = os.path.basename(month_file).split('_')[3].split('.')[0]
            var = os.path.basename(month_file).split('_')[2]
            source = 'WorldClim'
            if var == 'pet': 
                source = 'WorldClim (derived, Oudin et al., 2005)'
                var = 'pet2' # overwrite so we can distinguish this from pet1 (RDRS)
            l_index += [('Climate', f'{var}_mean_month_{month}', f'{sub_folder_unit}',  source),
                        ('Climate', f'{var}_std_month_{month}', f'{sub_folder_unit}', source)]

    return l_values, l_index

## ------- Component functions
# Shapefile subsetting for distributed attributes
def subset_glhymps_to_subbasin(glhymps, subbasin, ea_crs):
    poly_glhymps = gpd.clip(glhymps,subbasin)
    poly_glhymps['New_area_m2'] = poly_glhymps.to_crs(ea_crs).area
    return poly_glhymps

def subset_hydrolakes_to_subbasin(lakes, subbasin, ea_crs):

    '''lakes: GeoDataframe | subbasin: geometry | ea_crs: string'''

    # Clip the lakes polygon to the subbasin
    poly_lakes = gpd.clip(lakes,subbasin)

    # If we are left with any lake polygons:
    if len(poly_lakes) > 0 :
    # Update the Lake_area [km2] and Vol_total [million~m^3] values
    # This is needed in cases where due to clipping the lake polygon we end up with a partial lake in this subbasin
    
        # Old values
        old_areas = poly_lakes['Lake_area'] # km2
        old_volumes = poly_lakes['Vol_total'] # million m3
        
        # Get new area
        new_areas = poly_lakes.to_crs(ea_crs).area / 10**6 # [m2] / 10^6 = [km2]
        new_areas_rounded = round(new_areas,2) # this matches the number of significant digits in the test case (CAN_01DJ005)
        
        # Scale volume by new area - this is a bit simplistic but we have no better way to estimate the volume
        new_volumes = old_volumes * (new_areas_rounded / old_areas)
        
        # Replace values
        poly_lakes['Lake_area'] = new_areas_rounded
        poly_lakes['Vol_total'] = new_volumes

    return poly_lakes

# harmonic mean for conductivity averages
def harmonic_mean(x):
    return np.ma.count(x) / (1/x).sum()

def check_zonal_stats_outcomes(zonal_out, new_val=np.nan):
    '''Checks for None value in zonal_out, and sets these to np.nan or user-defined value'''
    for ix in range(0,len(zonal_out)):
         for key,val in zonal_out[ix].items():
             if val is None:
                 zonal_out[ix][key] = new_val
    return zonal_out

# Sine functions for P and T as per Woods (2009), Eq. 2 and 3. Period = 1 year
def sine_function_temp(x, mean, delta, phase, period=1):
    result = mean + delta * np.sin(2*np.pi*(x-phase)/period)
    return result

def sine_function_prec(x, mean, delta, phase, period=1):
    result = mean * (1+ delta* np.sin(2*np.pi*(x-phase)/period))
    results = np.where(result < 0, 0, result)
    return result

# Combined function for ERA5
def fit_climate_sines(t,p):

    # Define the x-coordinate as daily steps as fractions of a year
    x  = np.arange(0,len(p))/len(p) # days as fraction of a year
    
    # Fit the temperature sine
    t_mean   = float(t.mean())
    t_delta  = float(t.max()-t.min())
    t_phase  = 0.5
    initial_guess = [t_mean, t_delta, t_phase]
    t_pars, _ = curve_fit(sine_function_temp, x, t, p0=initial_guess)

    # Fit the precipitation sine
    p_mean   = float(p.mean())
    p_delta  = float(p.max()-p.min()) / float(p.mean())
    p_phase  = 0.5
    initial_guess = [p_mean, p_delta, p_phase]
    p_pars, _ = curve_fit(sine_function_prec, x, p, p0=initial_guess)

    # Compute seasonality as per Eq. 14 in Woods (2009)
    return p_pars[1]*np.sign(t_pars[1])*np.cos(2*np.pi*(p_pars[2]-t_pars[2])/1)

# Individual functions for stacked arrays in WorldClim
def fit_temp_sines(t):

    # Short-circuit the computation if we have only masked values
    if np.ma.getmaskarray(t).all():
        return [np.nan, np.nan, np.nan]
    
    # Define the x-coordinate as daily steps as fractions of a year
    x  = np.arange(0,len(t))/len(t) # days as fraction of a year
    
    # Fit the temperature sine
    t_mean   = float(t.mean())
    t_delta  = float(t.max()-t.min())
    t_phase  = 0.5
    initial_guess = [t_mean, t_delta, t_phase]
    t_pars, _ = curve_fit(sine_function_temp, x, t, p0=initial_guess)

    return t_pars

def fit_prec_sines(p):

    # Short-circuit the computation if we have only masked values
    if np.ma.getmaskarray(p).all():
        return [np.nan, np.nan, np.nan]
    
    # Define the x-coordinate as daily steps as fractions of a year
    x  = np.arange(0,len(p))/len(p) # days as fraction of a year

    # Fit the precipitation sine
    p_mean   = float(p.mean())
    p_delta  = float(p.max()-p.min()) / float(p.mean())
    p_phase  = 0.5
    initial_guess = [p_mean, p_delta, p_phase]
    p_pars, _ = curve_fit(sine_function_prec, x, p, p0=initial_guess)

    return p_pars

def get_geotif_nodata_value(tif):
    with rasterio.open(tif) as src:
        nodata_value = src.nodata
    return nodata_value

def check_scale_and_offset(tif):
    scale,offset = read_scale_and_offset(tif) # just to check we don't have any scale/offset going on
    if scale is None: scale = 1 # If scale is undefined that means we simply multiply by 1
    if offset is None: offset = 0 # Undefined offset > add 0
    if not (scale == 1) and not (offset == 0):
        print(f'--- WARNING: check_scale_and_offset(): scale or offset not 1 or 0 respectively.')
    return   

def get_categorical_dict(source):
    '''Contains dictionaries for categorical variables'''
    
    if source == 'GLCLU 2019':
        cat_dict = {1: 'true_desert',
                    2: 'semi_arid',
                    3: 'dense_short_vegetation',
                    4: 'open_tree_cover',
                    5: 'dense_tree_cover',
                    6: 'tree_cover_gain',
                    7: 'tree_cover_loss',
                    8: 'salt_pan',
                    9: 'wetland_sparse_vegetation',
                   10: 'wetland_dense_short_vegetation',
                   11: 'wetland_open_tree_cover',
                   12: 'wetland_dense_tree_cover',
                   13: 'wetland_tree_cover_gain',
                   14: 'wetland_tree_cover_loss',
                   15: 'ice',
                   16: 'water',
                   17: 'cropland',
                   18: 'built_up',
                   19: 'ocean',
                   20: 'no_data'}

    if source == 'MCD12Q1.061':
        cat_dict = {1: 'evergreen_needleleaf_forest',
                    2: 'evergreen_broadleaf_forest',
                    3: 'deciduous_needleleaf_forest',
                    4: 'deciduous_broadleaf_forest',
                    5: 'mixed_forest',
                    6: 'closed_shrubland',
                    7: 'open_shrubland',
                    8: 'woody_savanna',
                    9: 'savanna',
                   10: 'grassland',
                   11: 'permanent_wetland',
                   12: 'cropland',
                   13: 'urban_and_built_up',
                   14: 'cropland_natural_mosaic',
                   15: 'permanent_snow_ice',
                   16: 'barren',
                   17: 'water',
                  255: 'unclassified'}

    if source == 'LGRIP30':
        cat_dict = {0: 'water',
                    1: 'non_cropland',
                    2: 'irrigated_cropland',
                    3: 'rainfed_cropland'}
    
    return cat_dict

def update_values_list_with_categorical(l_values, l_index, zonal_out, source, prefix='', case='lumped'):
    '''Maps a zonal histogram of categorical classes onto descriptions and adds to lists'''

    # Get the category definitions
    cat_dict = get_categorical_dict(source)    

    # Separately handle lumped and distributed cases
    if case == 'lumped':
    
        # Find the total number of classified pixels
        total_pixels = 0
        for land_id,count in zonal_out[0].items():
            total_pixels += count
        
        # Loop over all categories and see what we have in this catchment
        for land_id,text in cat_dict.items():
            land_prct = 0
            if land_id in zonal_out[0].keys():
                land_prct = zonal_out[0][land_id] / total_pixels
            l_values.append(land_prct)
            l_index.append(('Land cover', f'{prefix}{text}_fraction', '-', f'{source}'))

    # distributed case, multiple polygons
    elif case == 'distributed': 

        # confirm that l_values has as many nested lists as we have zonal stats outputs
        num_nested_lists = sum(1 for item in l_values if isinstance(item, list))
        assert num_nested_lists == len(zonal_out), f"zonal_out length does not match expected list length {num_nested_lists}. zonal_out: {zonal_out}"

        # now loop over the zonal outputs and append to relevant lists
        for i in range(0,num_nested_lists):

            # Find the total number of classified pixels
            total_pixels = 0
            for land_id,count in zonal_out[i].items():
                total_pixels += count
            
            # Loop over all categories and see what we have in this catchment
            tmp_index = [] # we need this so the index resets on each subbasin iteration, and we need that because we only need the index once
            for land_id,text in cat_dict.items():
                land_prct = 0
                if land_id in zonal_out[i].keys():
                    land_prct = zonal_out[i][land_id] / total_pixels
                l_values[i].append(land_prct)
                tmp_index.append(('Land cover', f'{prefix}{text}_fraction', '-', f'{source}'))

        # Add the index values only once
        for item in tmp_index:
            l_index.append(item)

    return l_values,l_index


def zonal_stats_unit_conversion(zonal_out, stat_to_convert, variable, scale, offset):
    '''Takes a zonal_stats output and converts the units of any variable listed in stat_to_convert'''

    # Constants
    j_per_kj = 1000 # [J kJ-1]
    seconds_per_day = 24*60*60 # [s day-1]

    # Keep track of scale and offset
    # Update scale and offset to usable values - we get None if scale and offset are 1 and 0 in the GeoTIFF
    if scale is None: scale = 1 # If scale is undefined that means we simply multiply by 1
    if offset is None: offset = 0 # Undefined offset > add 0

    #  We'll need code to handle this if these aren't 1 and 0 respectively
    if (scale != 1) or (offset !=0):
        print(f'--- ERROR: zonal_stats_unit_conversion(): code needed to deal with scale {scale} and offset {offset}')
        return -1

    # Select conversion factor
    if variable == 'srad':
        # From [kJ m-2 day-1] to [W m-2]:
        # [kJ m-2 day-1] * 1/[s day-1] * [J kJ-1] = [J m-2 s-1] = [W m-2]
        c_factor = 1/seconds_per_day * j_per_kj

    # loop over all list elements
    for list_id in range(0,len(zonal_out)):
        zonal_dict = zonal_out[list_id]

        # Loop over dictionary entries
        for key,val in zonal_dict.items():
            if key in stat_to_convert:
                zonal_out[list_id][key] = zonal_out[list_id][key] * c_factor

    return zonal_out

def update_values_list(l_values, stats, zonal_out, scale, offset, case='lumped'):

    # Update scale and offset to usable values
    if scale is None: scale = 1 # If scale is undefined that means we simply multiply by 1
    if offset is None: offset = 0 # Undefined offset > add 0

    # Deal with the occassional None that we get when raster data input to zonal-stats is missing
    # This fixes 9 occurrences of missing data in the soilgrids maps
    zonal_out = zonal_out_none2nan(zonal_out)
    
    # We loop through the calculated stats in a pre-determined order:
    # 1. min
    # 2. mean (harmonic mean: 'harmonic_mean', circular mean: 'circ_mean')
    # 3. max
    # 4. stdev (circular stdev:  'circ_std'
    # 5. ..
    if case == 'lumped': # lumped case
        if 'min' in stats:  l_values.append(zonal_out[0]['min']  * scale + offset)
        if 'mean' in stats: l_values.append(zonal_out[0]['mean'] * scale + offset)
        if 'harmonic_mean'  in stats: l_values.append(zonal_out[0]['harmonic_mean'] * scale + offset) # SOILGRIDS conductivity
        if 'circ_mean'  in stats: l_values.append(zonal_out[0]['circ_mean'] * scale + offset) # MERIT Hydro aspect
        if 'max' in stats:  l_values.append(zonal_out[0]['max']  * scale + offset)
        if 'std' in stats:  l_values.append(zonal_out[0]['std']  * scale + offset)
        if 'circ_std' in stats:  l_values.append(zonal_out[0]['circ_std']  * scale + offset)
    
    # distributed case, multiple polygons
    elif case == 'distributed':
        
        # confirm that l_values has as many nested lists as we have zonal stats outputs
        num_nested_lists = sum(1 for item in l_values if isinstance(item, list))
        assert num_nested_lists == len(zonal_out), f"zonal_out length does not match expected list length {num_nested_lists}. zonal_out: {zonal_out}"
        
        # now loop over the zonal outputs and append to relevant lists
        for i in range(0,num_nested_lists):
            if 'min' in stats:  l_values[i].append(zonal_out[i]['min']  * scale + offset)
            if 'mean' in stats: l_values[i].append(zonal_out[i]['mean'] * scale + offset)
            if 'harmonic_mean'  in stats: l_values[i].append(zonal_out[i]['harmonic_mean'] * scale + offset) # only here for soilgrids conductivity, in which case we don't have 'mean'
            if 'circ_mean'  in stats: l_values[i].append(zonal_out[i]['circ_mean'] * scale + offset) # MERIT Hydro aspect
            if 'max' in stats:  l_values[i].append(zonal_out[i]['max']  * scale + offset)
            if 'std' in stats:  l_values[i].append(zonal_out[i]['std']  * scale + offset)
            if 'circ_std' in stats:  l_values[i].append(zonal_out[i]['circ_std']  * scale + offset)

    return l_values
    

def zonal_out_none2nan(zonal_out):
    '''Loops over zonal_out list that contains dictionaries, and replaces any None with np.nan'''
    for zonal_dict in zonal_out:
        for key in zonal_dict:
            if zonal_dict[key] is None:
                zonal_dict[key] = np.nan
                print(f'Replaced None in {key} with NaN')
    return zonal_out  

def read_scale_and_offset(geotiff_path):
    # Enforce data type
    geotiff_path = str(geotiff_path)
    
    # Open the GeoTIFF file
    dataset = gdal.Open(geotiff_path)

    if dataset is None:
        raise FileNotFoundError(f"File not found: {geotiff_path}")

    # Get the scale and offset values
    scale = dataset.GetRasterBand(1).GetScale()
    offset = dataset.GetRasterBand(1).GetOffset()

    # Close the dataset
    dataset = None

    return scale, offset

def get_geotif_data_as_array(file, band=1):
    file = str(file)
    ds = gdal.Open(file) # open the file
    band = ds.GetRasterBand(band) # get the data band
    data = band.ReadAsArray() # convert to numpy array for further manipulation   
    return data

def enforce_data_range(data,min,max,replace_with='limit'):

    '''Clamps data at min and max values'''

    if replace_with =='limit':
        data[data<min] = min
        data[data>max] = max
    else:
        data[data<min] = replace_with
        data[data>max] = replace_with
    
    return data

def write_geotif_sameDomain(src_file,des_file,des_data, nodata_value=None):
    # Enforce data type
    src_file = str(src_file)
    
    # load the source file to get the appropriate attributes
    src_ds = gdal.Open(src_file)
    
    # get the geotransform
    des_transform = src_ds.GetGeoTransform()

    # Get the scale factor from the source metadata
    scale_factor = src_ds.GetRasterBand(1).GetScale()
    offset = src_ds.GetRasterBand(1).GetOffset()
    
    # get the data dimensions
    ncols = des_data.shape[1]
    nrows = des_data.shape[0]
    
    # make the file
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(des_file,ncols,nrows,1,gdal.GDT_Float32, options = [ 'COMPRESS=DEFLATE' ])
    dst_ds.GetRasterBand(1).WriteArray( des_data )

    # Set the scale factor in the destination band
    if scale_factor: 
        dst_ds.GetRasterBand(1).SetScale(scale_factor)
    if offset: 
        dst_ds.GetRasterBand(1).SetOffset(offset)
    
    # Set the nodata value
    if nodata_value is not None:
        dst_ds.GetRasterBand(1).SetNoDataValue(nodata_value)

    # Set the geotransform
    dst_ds.SetGeoTransform(des_transform)

    # Set the projection
    wkt = src_ds.GetProjection()
    srs = osr.SpatialReference()
    srs.ImportFromWkt(wkt)
    dst_ds.SetProjection( srs.ExportToWkt() )
    
    # close files
    src_ds = None
    des_ds = None

    return

def subset_dataset_to_max_full_years(ds, time='time', res='hour', water_year=False, debug=False) -> xr.Dataset:

    # Find final and end years
    start_year = ds[time][0].dt.year.values
    final_year = ds[time][-1].dt.year.values
    final_timestamp = ds[time][-1]

    # Set how a year is defined
    if water_year: # Oct-1 to Sep-30
        start_month = 10
        start_day   = 1
        final_month = 9
        final_day   = 30
    else: # Jan-1 to Dec-31
        start_month = 1
        start_day   = 1
        final_month = 12
        final_day   = 31

    # Find the first occurrence of Jan-01 00:00 as start of the subset
    for year in range(start_year,final_year):
        
        # Define the initial timestamp to check
        if res == 'hour':
            start_timestamp = pd.Timestamp(year,start_month,start_day,0,0,0)
        elif res == 'day':
            start_timestamp = pd.Timestamp(year,start_month,start_day)
        if debug: print(f'checking {start_timestamp}')

        # Select the subset of the dataset for the current duration
        # Note: if either start or final are not part of the time series,
        #  this will silently just use whatever is available
        subset_ds = ds.sel(time=slice(start_timestamp, final_timestamp))

        # Check if we actually selected the duration we requested
        subset_start = pd.Timestamp(subset_ds[time][0].values)
        if subset_start == start_timestamp:
            subset_start_year = year # Keep track of where we start
            break # stop searching. We've found the first occurrence of Jan-01

    # Find the last occurrence of Dec-31 23:00 as end of the subset
    for year in range(final_year,subset_start_year,-1):

        # Define the initial timestamp to check
        if res == 'hour':
            end_timestamp = pd.Timestamp(year,final_month,final_day,23,0,0)
        elif res == 'day':
            end_timestamp = pd.Timestamp(year,final_month,final_day)
        if debug: print(f'checking {end_timestamp}')

        # Select the subset of the dataset for the current duration
        subset_ds = ds.sel(time=slice(start_timestamp, end_timestamp))

        # Check if we actually selected the duration we requested
        subset_end = pd.Timestamp(subset_ds[time][-1].values)
        if subset_end == end_timestamp:
            subset_end_year = year # Keep track of where we start
            break # stop searching. We've found the first occurrence of Jan-01

    # Now check if we have selected a zero-year period
    # This would imply we have less than a full year of data
    # In this case, just return the original data set with a warning
    if len(subset_ds) == 0:
        print(f'--- WARNING: subset_dataset_to_max_full_years(): Found no full data years. Returning original DataArray')
        return ds
    else:   
        return subset_ds

# Define a function to map months to seasons
def get_season(month):
    if month in [12, 1, 2]:
        return 'djf'
    elif month in [3, 4, 5]:
        return 'mam'
    elif month in [6, 7, 8]:
        return 'jja'
    elif month in [9, 10, 11]:
        return 'son'

# Finds duration counts in a vector of True/False values
def find_durations(condition):
    '''Counts the duration(s) of True values in an xarray dataseries'''

    previous = False
    duration = 0
    durations = []
    for flag in condition.values:
        
        # Time step where we reset the count
        if not previous and flag:
            duration = 0 # New first timestep where condition = True, so duration = 0
            previous = True
    
        # Time step where we're in a sequence of condition = True
        if previous and flag:
            duration += 1 # Update duration, implicitly retain previous = True by not changing it
        
        # Time step where we reach the end of a condition = True duration
        if previous and not flag:
            durations.append(duration) # Save the maximum duration length to list
            previous = False # Update previous; duration will be reset next time we encounter a condition = True
    
        # Time step where we're in a continuation of condition = False
        if not previous and not flag:
            continue # do nothing
    
    return np.array(durations)

## ---- HYDROLOGY
def calculate_flow_period_stats(var, condition, hilo, l_values, l_index,
                                dataset='ERA5', units='hours', category='Climate'):
    
    '''Calculates frequency (mean), duration (mean, median, skew, kurtosis) and
        timing of periods identified with a certain condition'''
    
    # Constants. We want everything in [days] for consistency with original CAMELS
    hours_per_day = 24 # [hours day-1]
    days_per_year = 365.25 # [days year-1]
    
    # Calculate frequencies
    freq = condition.mean(dim='time') * days_per_year # [-] * [days year-1]
    l_values.append(float(freq.values))
    l_index.append( (f'{category}', f'{hilo}_{var}_freq', 'days year^-1', dataset) )
    
    # Calculate duration statistics
    durations = find_durations(condition) # [time steps]
    if len(durations) == 0:
        dur_mean = 0
        dur_med = 0
        dur_skew = 0
        dur_kur = 0
    else:
        dur_mean = np.mean(durations)
        dur_med = np.median(durations)
        dur_skew = skew(durations)
        dur_kur = kurtosis(durations)
    l_values.append(dur_mean) # [days]
    l_index.append((f'{category}', f'{hilo}_{var}_dur_mean', 'days', dataset) ) # Consistency with
    l_values.append(dur_med) # [days]
    l_index.append((f'{category}', f'{hilo}_{var}_dur_median', 'days', dataset) )
    l_values.append(dur_skew) # [-]
    l_index.append((f'{category}', f'{hilo}_{var}_dur_skew', '-', dataset) )
    l_values.append(dur_kur) # [-]
    l_index.append((f'{category}', f'{hilo}_{var}_dur_kurtosis', '-', dataset) )
    
    # Calculate timing statistic
    condition['season'] = ('time', 
        [get_season(month) for month in condition['time.month'].values]) # add seasons
    season_groups = condition.groupby('season')
    season_list   = list(season_groups.groups.keys())
    max_season_id = int(season_groups.sum().argmax(dim='season').values) # find season with most True values
    l_values.append(season_list[max_season_id]) # add season abbrev
    l_index.append( (f'{category}', f'{hilo}_{var}_timing', 'season', dataset) )
    
    return l_values, l_index

def calculate_signatures(hyd, pre, source, l_values, l_index):
    '''Calculates various signatures'''

    ## prep for signatures that require precip
    # find the section of 'hyd' and 'pre' with the largest overlap
    pre_s = pre['time'][0]
    pre_e = pre['time'][-1]
    hyd_s = hyd['time'][0]
    hyd_e = hyd['time'][-1]
    pre_ss = pre.sel(time=slice(max(pre_s,hyd_s), min(pre_e, hyd_e)))
    hyd_ss = hyd.sel(time=slice(max(pre_s,hyd_s), min(pre_e, hyd_e))) # subset streamflow to precip record length
    assert hyd_ss['time'][0].values == pre_ss['time'][0].values, 'attributes_from_streamflow(): mismatch between precipitation and streamflow start timestamp' # confirm time periods are the same
    assert hyd_ss['time'][-1].values == pre_ss['time'][-1].values, 'attributes_from_streamflow(): mismatch between precipitation and streamflow final timestamp'
    assert len(hyd_ss['time']) == len(pre_ss['time']), 'attributes_from_streamflow(): different number of timesteps in precipitation and streamflow series'
    assert all(hyd_ss['water_year'] == pre_ss['water_year']) 

    ## LONG-TERM STATISTICS
    # Mean daily discharge
    daily_mean_q = hyd['q_obs'].groupby(hyd['water_year']).mean() # .mean() of daily values gives us [mm d-1]
    mean_q_m = daily_mean_q.mean()
    mean_q_s = daily_mean_q.std()
    l_values.append(float(mean_q_m.values))
    l_index.append( ('Hydrology', 'daily_discharge_mean', 'mm d^-1', f'{source}') )
    l_values.append(float(mean_q_s.values))
    l_index.append( ('Hydrology', 'daily_discharge_std', 'mm d^-1', f'{source}') )
    
    # Mean monthly flows
    monthly_q = hyd['q_obs'].resample(time='1M').mean().groupby('time.month')
    monthly_m = monthly_q.mean()
    with warnings.catch_warnings():
        # This mutes a "RuntimeWarning: Degrees of freedom <= 0 for slice. 
        #  var = nanvar(a, axis=axis, dtype=dtype, out=out, ddof=ddof,"
        #  warning which we'll get when calculating std() on NaN months
        warnings.simplefilter("ignore", category=RuntimeWarning)
        monthly_s = monthly_q.std()
    l_values, l_index = process_monthly_means_to_lists(monthly_m, 'mean', l_values, l_index, 'daily_streamflow', 'mm day^-1', 'Hydrology', source)
    l_values, l_index = process_monthly_means_to_lists(monthly_s, 'std',  l_values, l_index, 'daily_streamflow', 'mm day^-1', 'Hydrology', source)

    # Runoff ratio
    daily_mean_p = pre_ss.groupby('water_year').mean()
    daily_mean_q_ss = hyd_ss['q_obs'].groupby(hyd_ss['water_year']).mean()
    yearly_rr = daily_mean_q_ss/daily_mean_p
    rr_m = yearly_rr.mean()
    rr_s = yearly_rr.std()
    l_values.append(float(rr_m.values))
    l_index.append( ('Hydrology', 'runoff_ratio_mean', 'mm d^-1 month-1', f'{source}, RDRS') )
    l_values.append(float(rr_s.values))
    l_index.append( ('Hydrology', 'runoff_ratio_std', 'mm d^-1 month-1', f'{source}, RDRS') )
    
    # Streamflow elasticity
    q_elas = np.nanmedian(((daily_mean_q_ss - daily_mean_q_ss.mean())/(daily_mean_p - daily_mean_p.mean()))*(daily_mean_p.mean()/daily_mean_q_ss.mean()))
    l_values.append(q_elas)
    l_index.append( ('Hydrology', 'streamflow_elasticity', '-', f'{source}, RDRS') )

    # Slope of FDC
    groups = hyd['q_obs'].groupby(hyd['water_year'])
    slopes = []
    for year,group in groups:
        flows = group.values.copy()
        if np.logical_or(np.isnan(flows),flows == 0).all():
            slopes.append(np.nan) # so we can ignore these years
        else:
            # Account for NaNs - we don't want these to influence the calculations
            flows = flows[~np.isnan(flows)]
            # Account for zero and NaN flows that mess with log calculations
            flows[flows == 0] = flows.mean()/1000 # add 0.1% of mean flow to zeroes      
            # Do the actual stuff
            flows.sort()
            flows = np.log(flows)
            #slope = (np.percentile(flows,66) - np.percentile(flows,33)) / (.66*len(flows) - .33*len(flows))
            slope = (np.percentile(flows,66) - np.percentile(flows,33)) / (.66 - .33)
            slopes.append(slope)
    slopes = np.array(slopes)
    slope_m = np.nanmean(slopes)
    slope_s = np.nanstd(slopes)
    
    l_values.append(slope_m)
    l_index.append( ('Hydrology', 'fdc_slope_mean', '-', f'{source}') )
    l_values.append(slope_s)
    l_index.append( ('Hydrology', 'fdc_slope_std', '-', f'{source}') )
    
    # Baseflow index
    nan_mask = np.isnan(hyd['q_obs'])
    q_filled = hyd['q_obs'].interpolate_na(dim='time', method='linear')
    rec = baseflow.separation(q_filled.to_dataframe(), method='Eckhardt') # Find baseflow with Eckhardt filter
    tmp = xr.DataArray(rec['Eckhardt']['q_obs'])
    tmp[nan_mask] = np.nan
    hyd['q_bas'] = tmp
    daily_mean_qbase = hyd['q_bas'].groupby(hyd['water_year']).mean()
    bfi_m = (daily_mean_qbase / daily_mean_q).mean()
    bfi_s = (daily_mean_qbase / daily_mean_q).std()
    l_values.append(float(bfi_m.values))
    l_index.append( ('Hydrology', 'bfi_mean', '-', f'{source}') )
    l_values.append(float(bfi_s.values))
    l_index.append( ('Hydrology', 'bfi_std', '-', f'{source}') )
    
    # Half-flow date
    dates = []
    for year, group in hyd['q_obs'].groupby(hyd['water_year']):
        # Calculate the values for this year
        tmp_cum_flow = group.cumsum() # Cumulative flow per water year
        tmp_sum_flow = group.sum() # Total flow per water year
        tmp_frc_flow = tmp_cum_flow / tmp_sum_flow # Fractional flow per water year   
        # Deal with NaNs
        if ((group.isnull()) | (group == 0)).all(): # happens in a few basins with year-long+ periods of zero flow, and isnull() accounts for missing values in such years
            dates.append(np.nan)
        else:
            hdf = group[tmp_frc_flow > 0.5][0]['time'].values
            dates.append(pd.Timestamp(hdf).dayofyear)
    
    dates = np.array(dates)
    dates = dates[~np.isnan(dates)] # Remove NaNs so we can use the circ stats
    hdf_m = circmean(dates, high=366)
    hdf_s = circstd(dates, high=366)
    l_values.append(hdf_m)
    l_index.append( ('Hydrology', 'hfd_mean', 'day of year', f'{source}') )
    l_values.append(hdf_s)
    l_index.append( ('Hydrology', 'hfd_std', 'days', f'{source}') )
    
    # Quantiles
    groups = hyd['q_obs'].groupby(hyd['water_year'])
    with warnings.catch_warnings():
        # This mutes a "RuntimeWarning: All-NaN slice encountered
        #  return _nanquantile_unchecked("
        #  warning which we'll get when calculating quantiles on NaN years
        warnings.simplefilter("ignore", category=RuntimeWarning)
        for quantile in [0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]:
            q_m = groups.quantile(quantile, skipna=True).mean()
            q_s = groups.quantile(quantile).std()
            l_values.append(float(q_m.values))
            l_index.append(('Hydrology', f'q{int(quantile*100)}_mean', 'mm day^-1', f'{source}'))
            l_values.append(float(q_s.values))
            l_index.append(('Hydrology', f'q{int(quantile*100)}_std', 'mm day^-1', f'{source}'))
    
    # Durations
    variable  = 'q_obs'
    no_flow_threshold = 0 # mm d-1
    no_flow_condition = hyd[variable] <= no_flow_threshold
    l_values,l_index  = calculate_flow_period_stats('flow',no_flow_condition,'no',l_values,l_index,dataset=source,units='days',category='Hydrology')
    
    low_flow_threshold = 0.2 * hyd['q_obs'].mean() # mm d-1
    low_flow_condition = hyd[variable] < low_flow_threshold
    l_values,l_index   = calculate_flow_period_stats('flow',low_flow_condition,'low',l_values,l_index,dataset=source,units='days',category='Hydrology')
    
    high_flow_threshold = 9 * hyd['q_obs'].median() # mm d-1
    high_flow_condition = hyd[variable] > high_flow_threshold
    l_values,l_index    = calculate_flow_period_stats('flow',high_flow_condition,'high',l_values,l_index,dataset=source,units='days',category='Hydrology')
    
    return l_values,l_index

## ---- HydroLAKES
def get_open_water_stats(gdf, att, mask, l_values, l_index):
    '''Calculates min, mean, max, std and total att ('Lake_area', 'Vol_total'), optionally using a reservoir mask'''

    # Initialization - this must happen before we check if we got a GDF,
    #   because otherwise these values won't be available for l_index
    if mask == 'all':
        water = 'open_water' # no reservoir mask, hence we're working with lakes
    elif mask == 'reservoir':
        water = 'reservoir'

    if att == 'Lake_area':
        units = 'km^2' # https://data.hydrosheds.org/file/technical-documentation/HydroLAKES_TechDoc_v10.pdf
        att_d = 'area'
    elif att == 'Vol_total':
        units = 'million~m^3' # https://data.hydrosheds.org/file/technical-documentation/HydroLAKES_TechDoc_v10.pdf
        att_d = 'volume'

    # Check if we were handed a lake polygon
    if gdf is not None:
        # Setup
        if mask == 'all':
            mask = gdf.index >= 0 # Selects everything if no mask is provided
        elif mask == 'reservoir':
            mask = gdf['Lake_type'] == 2
        
        # Get the values
        min_val = gdf[mask][att].min()
        mean_val = gdf[mask][att].mean()
        max_val = gdf[mask][att].max()
        std_val = gdf[mask][att].std()
        tot_val = gdf[mask][att].sum()# total
    else:
        min_val = 0
        mean_val = 0
        max_val = 0
        std_val = 0
        tot_val = 0

    # Stats
    l_values.append(min_val)
    l_index.append(('Open water', f'{water}_{att_d}_min',  f'{units}', 'HydroLAKES'))
    l_values.append(mean_val)
    l_index.append(('Open water', f'{water}_{att_d}_mean',  f'{units}', 'HydroLAKES'))
    l_values.append(max_val)
    l_index.append(('Open water', f'{water}_{att_d}_max',  f'{units}', 'HydroLAKES'))
    l_values.append(std_val)
    l_index.append(('Open water', f'{water}_{att_d}_std',  f'{units}', 'HydroLAKES'))
    l_values.append(tot_val) 
    l_index.append(('Open water', f'{water}_{att_d}_total',  f'{units}', 'HydroLAKES'))
    
    return l_values, l_index

## ---- MERIT
# Circular stats for aspect
# # need compressed() because scipy < 1.12.0 circmean/std don't work well with masked arrays. This flattens the array and keeps values only
def calc_circmean(x):
    return circmean(x.compressed(), high=360) 

def calc_circstd(x):
    return circstd(x.compressed(), high=360)

# deprecated, keeping for posterity
def get_aspect_attributes(tif,l_values,l_index):
    '''Calculates circular statistics for MERIT Hydro aspect'''

    # Get data as a masked array - we know there are no-data values outside the catchment boundaries
    aspect = get_geotif_data_as_array(tif)
    no_data = get_geotif_nodata_value(tif)
    masked_aspect = np.ma.masked_array(aspect, aspect == no_data)
    values = masked_aspect.compressed() # scipy < 1.12.0 circmean/std don't work well with masked arrays

    ## Calculate the statistics
    l_values.append(values.min())
    l_index.append(('Topography', 'aspect_min',  'degrees', 'MERIT Hydro'))
    
    l_values.append(circmean(values,high=360))
    l_index.append(('Topography', 'aspect_mean',  'degrees', 'MERIT Hydro'))

    l_values.append(circstd(values,high=360))
    l_index.append(('Topography', 'aspect_std',  'degrees', 'MERIT Hydro'))

    l_values.append(values.max())
    l_index.append(('Topography', 'aspect_max',  'degrees', 'MERIT Hydro'))
    
    return l_values,l_index

def get_river_attributes(riv_str, bas_str, l_values, l_index, area, equal_area_crs='ESRI:102008'):
    
    '''Calculates topographic attributes from a MERIT Hydro Basins river polygon'''

    # We need some special actions if we're dealing with the distributed case, so check that first
    case = 'lumped'
    if 'distributed' in bas_str:
        case = 'distributed'

    # Handle the cases
    if case == 'lumped': 
    
        # Initialize NaN for cases where we have no shapefile or an empty one
        stream_total = np.nan
        min_length = np.nan
        mean_length = np.nan
        max_length = np.nan
        std_length = np.nan
        riv_order = np.nan
        density = np.nan
        elongation = np.nan
        
        # Check if the file exists (for headwaters we won't have a river polygon)
        if os.path.exists(riv_str):
            
            # Load shapefiles
            river = gpd.read_file(riv_str)
            river = river.set_index('COMID')
            river = river[~river.index.duplicated(keep='first')] # Removes any duplicate river segments
    
            # Check if we actually have a river segment (empty shapefile is not the same as no shapefile)
            if len(river) > 0:
                # Raw data
                stream_lengths = []
                headwaters = river[(river['maxup'] == 0) | (river['maxup'].isna())] # identify reaches with no upstream
                for COMID in headwaters.index:
                    stream_length = 0
                    while COMID in river.index:
                        stream_length += river.loc[COMID]['new_len_km'] # Add the length of the current segment
                        COMID = river.loc[COMID]['NextDownID'] # Get the downstream reach
                    stream_lengths.append(stream_length) # If we get here we ran out of downstream IDs
            
                # Stats
                stream_total = river['new_len_km'].sum()
                stream_lengths = np.array(stream_lengths)
                min_length = stream_lengths.min()
                mean_length = stream_lengths.mean()
                max_length = stream_lengths.max()
                std_length = stream_lengths.std()
                riv_order = river['order'].max()
                density = stream_total/area
                elongation = 2*np.sqrt(area/np.pi)/max_length
    
        # Update lists
        l_values.append(min_length)
        l_index.append(('Topography', 'stream_length_min',  'km', 'MERIT Hydro Basins'))
        l_values.append(mean_length)
        l_index.append(('Topography', 'stream_length_mean', 'km', 'MERIT Hydro Basins'))
        l_values.append(max_length)
        l_index.append(('Topography', 'stream_length_max',  'km', 'MERIT Hydro Basins'))
        l_values.append(std_length)
        l_index.append(('Topography', 'stream_length_std',  'km', 'MERIT Hydro Basins'))
        l_values.append(stream_total)
        l_index.append(('Topography', 'segment_length_total', 'km', 'MERIT Hydro Basins'))
        
        # Order
        l_values.append(riv_order)
        l_index.append(('Topography', 'stream_order_max',  '-', 'MERIT Hydro Basins'))
    
        # Derived
        l_values.append(density)
        l_index.append(('Topography', 'stream_density',  'km^-1', 'MERIT Hydro, MERIT Hydro Basins'))
        l_values.append(elongation)
        l_index.append(('Topography', 'elongation_ratio','-', 'MERIT Hydro, MERIT Hydro Basins'))

        # Final empty output
        merit_comids = None

    # distributed case
    elif case == 'distributed':

        # Get the basin shape and figure out what we're dealing with
        basin = gpd.read_file(bas_str)
        num_poly = len(basin)
        merit_comids = basin['COMID'].values

        # Check if the file exists (for headwaters we won't have a river polygon)
        if os.path.exists(riv_str):
            
            # Load shapefiles
            river = gpd.read_file(riv_str)
            river = river.set_index('COMID')
            river = river[~river.index.duplicated(keep='first')] # Removes any duplicate river segments

            # Loop over the basins and process
            for i_poly in range(num_poly):

                # Find the river segment, if any, and get values for the statistics
                comid = basin.iloc[i_poly]['COMID']
                if comid in river.index:
                    riv_poly = river.loc[comid]
                    stream_length = riv_poly['new_len_km']
                    slope = riv_poly['slope']
                    upstream_area = riv_poly['uparea']
                    basin_area = basin.to_crs('ESRI:102008').area.iloc[0] / 10**6 # area in km2
                    density = stream_length/basin_area
                    elongation = 2*np.sqrt(basin_area/np.pi)/stream_length
                else:
                    stream_length = np.nan # no delineated river here, so no stats
                    slope = np.nan
                    upstream_area = np.nan
                    density = np.nan
                    elongation = np.nan
                
                # Append statistics to list
                l_values[i_poly].append(stream_length)
                l_values[i_poly].append(slope)
                l_values[i_poly].append(upstream_area)
                l_values[i_poly].append(density)
                l_values[i_poly].append(elongation)
                    
        else: # no river shapefile
            
            # just add nans everywhere
            stream_length = np.nan # no delineated river here, so no stats
            slope = np.nan
            upstream_area = np.nan
            density = np.nan
            elongation = np.nan
    
            for i_poly in range(num_poly):
                l_values[i_poly].append(stream_length)
                l_values[i_poly].append(slope)
                l_values[i_poly].append(upstream_area)
                l_values[i_poly].append(density)
                l_values[i_poly].append(elongation)

        # Update the index list
        l_index += [('Topography', 'stream_length',    'km',     'MERIT Hydro Basins'),
                    ('Topography', 'stream_slope',     'm m^-1', 'MERIT Hydro Basins'),
                    ('Topography', 'upstream_area',    'km^2',   'MERIT Hydro Basins'),
                    ('Topography', 'stream_density',   'km^-1',  'MERIT Hydro, MERIT Hydro Basins'),
                    ('Topography', 'elongation_ratio', '-',      'MERIT Hydro, MERIT Hydro Basins')]
    
    return l_values,l_index,merit_comids

## ---- RDRS
def find_climate_seasonality_rdrs(ds, use_typical_cycle=False):

    if not use_typical_cycle:
        # Resample the observations to daily, retain individual years
        daily_p_groups = ds['RDRS_v2.1_A_PR0_SFC'].resample(time='1D').mean().groupby('time.year')
        daily_t_groups = ds['RDRS_v2.1_P_TT_1.5m'].resample(time='1D').mean().groupby('time.year')

        seasonalities = []
        for zip_p, zip_t in zip(daily_p_groups,daily_t_groups):
            year = zip_p[0]
            yp = zip_p[1].values.flatten()
            yt = zip_t[1].values.flatten()
            seasonalities.append(fit_climate_sines(yt,yp))
        return np.array(seasonalities)
    
    else:
        # Resample to typical seasonal cycles at daily resolution
        daily_p = ds['RDRS_v2.1_A_PR0_SFC'].resample(time='1D').mean().groupby('time.dayofyear').mean()
        daily_t = ds['RDRS_v2.1_P_TT_1.5m'].resample(time='1D').mean().groupby('time.dayofyear').mean()

        # Time series of input and output
        yp = daily_p.values.flatten()
        yt = daily_t.values.flatten()

        # Fit the sine curves and return the seasonality coefficient
        return fit_climate_sines(yt,yp)



## ---- ERA5
def find_climate_seasonality_era5(ds, use_typical_cycle=False):

    if not use_typical_cycle:
        # Resample the observations to daily, retain individual years
        daily_p_groups = ds['mtpr'].resample(time='1D').mean().groupby('time.year')
        daily_t_groups = ds['t'].resample(time='1D').mean().groupby('time.year')

        seasonalities = []
        for zip_p, zip_t in zip(daily_p_groups,daily_t_groups):
            year = zip_p[0]
            yp = zip_p[1].values.flatten()
            yt = zip_t[1].values.flatten()
            seasonalities.append(fit_climate_sines(yt,yp))
        return np.array(seasonalities)
    
    else:
        # Resample to typical seasonal cycles at daily resolution
        daily_p = ds['mtpr'].resample(time='1D').mean().groupby('time.dayofyear').mean()
        daily_t = ds['t'].resample(time='1D').mean().groupby('time.dayofyear').mean()

        # Time series of input and output
        yp = daily_p.values.flatten()
        yt = daily_t.values.flatten()

        # Fit the sine curves and return the seasonality coefficient
        return fit_climate_sines(yt,yp)

def calculate_temp_prcp_stats(var, condition, hilo, l_values,l_index,
                              dataset='ERA5', units='hours'):
    
    '''Calculates frequency (mean) and duration (mean, median, skew, kurtosis) 
        of temperature/precipitation periods'''

    # Constants. We want everything in [days] for consistency with original CAMELS
    hours_per_day = 24 # [hours day-1]
    days_per_year = 365.25 # [days year-1]

    # Calculate frequencies
    freq = condition.mean(dim='time') * days_per_year # [-] * [days year-1]
    l_values.append(freq.values)
    l_index.append( ('Climate', f'{hilo}_{var}_freq', 'days year^-1', dataset) )
    
    # Calculate duration statistics
    durations = find_durations(condition) # [time steps]
    if units == 'hours':
        durations = durations / hours_per_day # [days] = [hours] / [hours day-1]
    l_values.append(np.mean(durations)) # [days]
    l_index.append( ('Climate', f'{hilo}_{var}_dur_mean', 'days', dataset) ) # Consistency with
    l_values.append(np.median(durations)) # [days]
    l_index.append( ('Climate', f'{hilo}_{var}_dur_median', 'days', dataset) )
    l_values.append(skew(durations)) # [-]
    l_index.append( ('Climate', f'{hilo}_{var}_dur_skew', '-', dataset) )
    l_values.append(kurtosis(durations)) # [-]
    l_index.append( ('Climate', f'{hilo}_{var}_dur_kurtosis', '-', dataset) )

    # Calculate timing statistic
    condition['season'] = ('time', [get_season(month) for month in condition['time.month'].values]) # add seasons
    season_groups = condition.groupby('season')
    season_list   = list(season_groups.groups.keys())
    max_season_id = int(season_groups.sum().argmax(dim='season').values) # find season with most True values
    l_values.append(season_list[max_season_id]) # add season abbrev
    l_index.append( ('Climate', f'{hilo}_{var}_timing', 'season', dataset) )

    return l_values, l_index

# Finds a long-term mean daily maximum temperature in a way that doesn't rely on 
#  time.dt.dayofyear, because the latter will have Dec-31 as day 365 in non-leap-years,
#  and as 366 in leap years. Hence the day-of-year means do not use the same dates for 
#  a given DoY in leap years as they do in regular years.
def create_mean_daily_max_series(ds,var='t'):
    '''Finds the long-term mean daily maximum value of a variable'''
    
    # Create an array of all the month-days we have (e.g. 1949-12-31 00:00 becomes 1231)
    month_days_all = ds.time.dt.month * 100 + ds.time.dt.day

    # Loop over the unique month-days we have, and find the mean daily maximum value for each
    month_days_unique = np.unique(month_days_all)
    mean_daily_max = []
    for month_day in month_days_unique:
        val = ds[var].sel(time=(month_days_all==month_day)).groupby('time.year').max().mean().values
        mean_daily_max.append(val)

    # Convert the list to an array for further processing
    mean_daily_max = np.array(mean_daily_max)

    # Extract month_day values from the long xarray DataArray
    month_day_values = month_days_all.values
    
    # Find the indices of each month_day in the unique_month_days array
    indices = np.searchsorted(month_days_unique, month_day_values)
    
    # Use the indices to extract the corresponding data values
    corresponding_data_values = mean_daily_max[indices]
    
    # Create a new DataArray with the corresponding data values
    result_array = xr.DataArray(corresponding_data_values, 
                                coords=month_days_all.coords, 
                                dims=month_days_all.dims)

    return result_array


# Processing function to update the two main lists we're populating
def process_monthly_means_to_lists(da, stat, l_values, l_index, var, unit, category='Climate', source='ERA5'):
    '''Takes an xarray data array with monthly statistics and processes into l_values and l_index lists'''
    for month in range(1,13):
        val = da.sel(month=month).values.flatten()[0]
        txt = (f'{category}', f'{var}_{stat}_month_{month:02}', f'{unit}', f'{source}')
        l_values += [val] # Needs to be this way because we're appending to a list
        l_index  += [txt]
    return l_values, l_index

# Define custom functions to apply circmean and circstd to Xarray group objects
# Without this, xarray chokes on dimensions when converting the circmean/circstd output 
#  back into something with the right month indices
def circmean_group(group):
    return xr.DataArray(circmean(group, high=360, low=0), name='phi')
def circstd_group(group):
    return xr.DataArray(circstd(group, high=360, low=0), name='phi')

## ---- LAI
def filter_lai_files_by_date(files, last_n_years=[], last_n_months=[], last_n_days=[],
                                    years=[], months=[], days=[]):

    '''Filters list of LAI file names by last n years/months/days and/or by year/month/day x.
       Assumes date is given as 'yyyymmdd_*.tif', as part of the filename.
       Use years/months/days (input as list) to subset further.'''

    # Check inputs
    if (last_n_years and last_n_months) or \
       (last_n_years and last_n_days) or \
       (last_n_months and last_n_days):
        print('WARNING: filter_lai_files_by_date(): specify only one of last_n_years, last_n_months, last_n_days')
        return

    # Create a DatetimeIndex from filenames
    dates = []
    for file in files:
        file_name = os.path.basename(file)
        yyyymmdd = file_name[0:8]
        dates.append(yyyymmdd)
    dti = pd.to_datetime(dates,format='%Y%m%d')

    # Set the first entry
    start_date = dti[0]
    
    # Find the last entry
    last_year  = dti[-1].year
    last_month = dti[-1].month
    last_day   = dti[-1].day
    
    # Select the last n entries
    if last_n_years:    start_date = dti[-1] - relativedelta(years = last_n_years)
    elif last_n_months: start_date = dti[-1] - relativedelta(months = last_n_months)
    elif last_n_days:   start_date = dti[-1] - relativedelta(days = last_n_days)
    last_n = (dti >= start_date) & (dti <= dti[-1])

    # Specify filters to include all if no specific years/months/days were requested
    if not years:  years  = list(set(dti.year))  # i.e. filter to include all unique years in dti, \
    if not months: months = list(set(dti.month)) #    else use user input
    if not days:   days   = list(set(dti.day))
    mask = dti.year.isin(years) & dti.month.isin(months) & dti.day.isin(days)

    # Return the filtered list
    return [file for file, bool1, bool2 in zip(files,last_n,mask) if bool1 and bool2]

def calculate_monthly_lai_maps(lai_files, des_path):
    des_files = []
    for month in range(1,13):

        # Define valid data range
        # See docs, Table 4: https://lpdaac.usgs.gov/documents/926/MOD15_User_Guide_V61.pdf
        modis_min = 0
        modis_max = 100
    
        # Get the files we have for this month, for the last n years
        #print(f'Processing month {month:02d}')
        month_files = filter_lai_files_by_date(lai_files, months=[month])
    
        # Remove the one file we know is incomplete, 2022-10-16
        month_files = [file for file in month_files if '20221016' not in file]
        
        # Load the data as numpy arrays, stack vertically, and find the mean value (ignoring nan)
        data = [get_geotif_data_as_array(file) for file in month_files] # Get data as uint8
        stacked = np.dstack(data) # Create a 3D stack
        stacked_msk = np.ma.masked_array(stacked, mask=(stacked<modis_min) | (stacked>modis_max)) # Retain valid values only
        mean_lai = np.ma.mean(stacked_msk, axis=2)
    
        # Define the no-data locations
        #mean_all = np.nanmean(stacked, axis=2) # Any pixel that consistently has no-data in the source files (>= 249) should have a >= 249 mean
        #mean_lai[mean_all >= 249] = mean_all[mean_all >= 249] # Place the no-data values in the new monthly-mean-lai file
        
        # Define output file name and write to disk
        src_file = month_files[0] # We use this to copy over domain, projection, data scaling, etc
        des_file = str( des_path / f'month_mean_{month:02d}_MOD_Grid_MOD15A2H_Lai_500m.tif' )
        write_geotif_sameDomain(src_file, des_file, mean_lai)
        des_files.append(des_file)
    return des_files

## --- WorldClim
def aridity_and_fraction_snow_from_worldclim(geo_folder, dataset, overwrite=False):
    
    '''Calculates aridity and fraction snow maps from WorldClim data'''

    # Find files
    clim_folder = geo_folder / dataset / 'raw'
    prc_files = sorted( glob.glob(str(clim_folder / 'prec' / '*.tif')) ) # [mm]
    pet_files = sorted( glob.glob(str(clim_folder / 'pet' / '*.tif')) ) # [mm]
    tmp_files = sorted( glob.glob(str(clim_folder / 'tavg' / '*.tif')) ) # [C]
    
    # Make the output locations
    ari_folder = clim_folder / 'aridity2'
    ari_folder.mkdir(parents=True, exist_ok=True)
    snow_folder = clim_folder / 'snow2'
    snow_folder.mkdir(parents=True, exist_ok=True)    
    fsnow_folder = clim_folder / 'fracsnow2'
    fsnow_folder.mkdir(parents=True, exist_ok=True)

    # Loop over files and calculate aridity
    for prc_file, pet_file, tmp_file in zip(prc_files, pet_files, tmp_files):

        # Define output file names
        ari_name = prc_file.replace('prec','aridity2')
        ari_file = str(ari_folder / ari_name)
        fsnow_name = prc_file.replace('prec','fracsnow2')
        fsnow_file = str(fsnow_folder / fsnow_name)
        snow_name = prc_file.replace('prec','snow2')
        snow_file = str(snow_folder / snow_name)

        # Check if we proceed
        if os.path.isfile(ari_file) and os.path.isfile(fsnow_file) and os.path.isfile(snow_file) and not overwrite:
            return # step out early

        # Define month
        month = prc_file.split('_')[-1].split('.')[0] # 'wc2.1_30s_prec_01.tif' > '01', ..., '12'
        month_ix = int(month)-1 # -1 to account for zero-based indexing: Jan value is at index 0, not 1

        # Load data
        prc_path = str(clim_folder / 'prec' / prc_file)
        pet_path = str(clim_folder / 'pet'  / pet_file)
        tmp_path = str(clim_folder / 'tavg' / tmp_file )     
        prc = get_geotif_data_as_array(prc_path) # [mm]
        pet = get_geotif_data_as_array(pet_path) # [mm]
        tmp = get_geotif_data_as_array(tmp_path) # [C]

        # Calculate variables
        snow = np.where(tmp < 0, prc, 0) # get snow first, because this needs precip and we'll (possibly) be updating the precip value below
        if (prc == 0).any():
            prc[prc == 0] = 1 # add 1 mm to avoid divide by zero errors
        ari = pet/prc # [-]
        frac_snow = snow/prc # [-]

        # Write to disk
        write_geotif_sameDomain(prc_path, ari_file, ari)
        write_geotif_sameDomain(prc_path, fsnow_file, frac_snow)
        write_geotif_sameDomain(prc_path, snow_file, snow)
    
    return

def oudin_pet_from_worldclim(geo_folder, dataset, debug=False, overwrite=False):

    '''Calculates PET estimates from WorldClim data, using the Oudin (2005; 10.1016/j.jhydrol.2004.08.026) formulation'''

    # Constants
    lh = 2.45 # latent heat flux, MJ kg-1
    rw = 1000 # rho water, kg m-3
    days_per_month = [31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31] # days month-1
    mm_per_m = 1000 # mm m-1
    
    # Find files
    clim_folder = geo_folder / dataset / 'raw'
    srad_files = sorted( glob.glob(str(clim_folder / 'srad' / '*.tif')) ) # website says [kJ m-2 day-1], but paper says [MJ m-2 day-1]
    tavg_files = sorted( glob.glob(str(clim_folder / 'tavg' / '*.tif')) ) # C

    # Make the output location
    pet_folder = clim_folder / 'pet'
    pet_folder.mkdir(parents=True, exist_ok=True)

    # Loop over files and calculate PET
    for srad_file, tavg_file in zip(srad_files, tavg_files):

        # Define output file name
        pet_name = srad_file.replace('srad','pet')
        pet_file = str(pet_folder / pet_name)

        # Check if we proceed
        if os.path.isfile(pet_file) and not overwrite:
            return # step out early

        # Define month
        month = srad_file.split('_')[-1].split('.')[0] # 'wc2.1_30s_srad_01.tif' > '01', ..., '12'
        month_ix = int(month)-1 # -1 to account for zero-based indexing: Jan value is at index 0, not 1
        
        # Load data
        srad_path = str(clim_folder / 'srad' / srad_file)
        tavg_path = str(clim_folder / 'tavg' / tavg_file)     
        srad = get_geotif_data_as_array(srad_path) / 1000 # [kJ m-2 day-1] / 1000 = [MJ m-2 day-1]
        tavg = get_geotif_data_as_array(tavg_path)
        
        # Oudin et al, 2005, Eq. 3
        pet = np.where(tavg+5 > 0, (srad / (lh*rw)) * ((tavg+5)/100) * mm_per_m, 0) # m day-1 > mm day-1
        pet_month = pet * days_per_month[month_ix] # mm month-1
        if debug: print(f'Calculating monthly PET for month {month} at day-index {month_ix}')

        # Write to disk
        write_geotif_sameDomain(srad_path, pet_file, pet_month)
    return

def get_annual_worldclim_attributes(clim_folder, shp_str, prec_folder, tavg_folder, pet_folder, snow_folder, l_values, l_index):

    '''Calculates annual WorldClim statistics'''

    # Create the output folder
    ann_folder = clim_folder / 'annual'
    ann_folder.mkdir(exist_ok=True, parents=True)

    # General settings
    stats = ['mean', 'std']

    # We need some special actions if we're dealing with the distributed case, so check that first
    case = 'lumped'
    if 'distributed' in shp_str:
        case = 'distributed'
        num_nested_lists = len(l_values)
    
    # --- P
    prec_files = sorted( glob.glob(str(clim_folder / prec_folder / '*.tif')) )
    annual_prec = create_annual_worldclim_map(prec_files, ann_folder, 'prec_sum.tif', 'prec')
    zonal_out = zonal_stats(shp_str, annual_prec, stats=stats, all_touched=True)

    # Handle scale and offset
    scale,offset = read_scale_and_offset(annual_prec)
    if scale is None: scale = 1 # If scale is undefined that means we simply multiply by 1
    if offset is None: offset = 0 # Undefined offset > add 0

    # Update lists
    for stat in stats:
        if case == 'lumped':
            l_values.append(zonal_out[0][stat] * scale + offset)
        elif case == 'distributed':
            for i in range(0,num_nested_lists):
                l_values[i].append(zonal_out[i][stat] * scale + offset)
        l_index.append(('Climate',f'prec_{stat}','mm', 'WorldClim'))

    # --- PET
    pet_files = sorted( glob.glob(str(clim_folder / pet_folder / '*.tif')) )
    annual_pet = create_annual_worldclim_map(pet_files, ann_folder, 'pet_sum.tif', 'pet')
    zonal_out = zonal_stats(shp_str, annual_pet, stats=stats, all_touched=True)
    
    scale,offset = read_scale_and_offset(annual_pet)
    if scale is None: scale = 1
    if offset is None: offset = 0
        
    for stat in stats:
        if case == 'lumped':
            l_values.append(zonal_out[0][stat] * scale + offset)
        elif case == 'distributed':
            for i in range(0,num_nested_lists):
                l_values[i].append(zonal_out[i][stat] * scale + offset)
        l_index.append(('Climate',f'pet2_{stat}','mm', 'WorldClim'))
        
    # --- T
    tavg_files = sorted( glob.glob(str(clim_folder / tavg_folder / '*.tif')) )
    annual_tavg = create_annual_worldclim_map(tavg_files, ann_folder, 't_avg.tif', 'tavg')
    zonal_out = zonal_stats(shp_str, annual_tavg, stats=stats, all_touched=True)
    
    scale,offset = read_scale_and_offset(annual_tavg)
    if scale is None: scale = 1
    if offset is None: offset = 0
        
    for stat in stats:
        if case == 'lumped':
            l_values.append(zonal_out[0][stat] * scale + offset)
        elif case == 'distributed':
            for i in range(0,num_nested_lists):
                l_values[i].append(zonal_out[i][stat] * scale + offset)
        l_index.append(('Climate',f'tavg_{stat}','C', 'WorldClim'))

    # --- Snow
    snow_files = sorted( glob.glob(str(clim_folder / snow_folder / '*.tif')) )
    annual_snow = create_annual_worldclim_map(snow_files, ann_folder, 'snow_sum.tif', 'snow')

    # --- Aridity
    annual_ari = derive_annual_worldclim_aridity(annual_prec, annual_pet, ann_folder, 'aridity.tif')
    zonal_out = zonal_stats(shp_str, annual_ari, stats=stats, all_touched=True)
    
    scale,offset = read_scale_and_offset(annual_ari)
    if scale is None: scale = 1
    if offset is None: offset = 0
        
    for stat in stats:
        if case == 'lumped':
            l_values.append(zonal_out[0][stat] * scale + offset)
        elif case == 'distributed':
            for i in range(0,num_nested_lists):
                l_values[i].append(zonal_out[i][stat] * scale + offset)
        l_index.append(('Climate',f'aridity2_{stat}','-', 'WorldClim'))

    # --- Seasonality
    annual_seas = derive_annual_worldclim_seasonality(prec_files, tavg_files, ann_folder, 'seasonality.tif')
    zonal_out = zonal_stats(shp_str, annual_seas, stats=stats, all_touched=True)
    
    scale,offset = read_scale_and_offset(annual_seas)
    if scale is None: scale = 1
    if offset is None: offset = 0
        
    for stat in stats:
        if case == 'lumped':
            l_values.append(zonal_out[0][stat] * scale + offset)
        elif case == 'distributed':
            for i in range(0,num_nested_lists):
                l_values[i].append(zonal_out[i][stat] * scale + offset)
        l_index.append(('Climate',f'seasonality2_{stat}','-', 'WorldClim'))
    
    # --- Snow fraction
    annual_fs = derive_annual_worldclim_fracsnow(annual_prec, annual_snow, ann_folder, 'fracsnow.tif')
    zonal_out = zonal_stats(shp_str, annual_fs, stats=stats, all_touched=True)
    
    scale,offset = read_scale_and_offset(annual_fs)
    if scale is None: scale = 1
    if offset is None: offset = 0
        
    for stat in stats:
        if case == 'lumped':
            l_values.append(zonal_out[0][stat] * scale + offset)
        elif case == 'distributed':
            for i in range(0,num_nested_lists):
                l_values[i].append(zonal_out[i][stat] * scale + offset)
        l_index.append(('Climate',f'fracsnow2_{stat}','-', 'WorldClim'))
    
    return l_values, l_index

def create_annual_worldclim_map(files, output_folder, output_name, var):
    '''Creates an annual map from monthly WorldClim data'''

    # Load the data and mask the no-data values
    datasets = [get_geotif_data_as_array(file) for file in files]
    no_datas = [get_geotif_nodata_value(file) for file in files]
    masked_data = [np.ma.masked_array(data, mask=data==no_data) for data,no_data in zip(datasets,no_datas)]

    # Create stacked array for computations along the third dimension
    stacked_array = np.ma.dstack(masked_data)

    # Create the annual value depending on what we're after
    if var in ['prec','pet','snow']:
        res_array = np.ma.sum(stacked_array, axis=2) # Sum along the third dimension
    if var in ['tavg']:
        res_array = np.ma.mean(stacked_array, axis=2) # Sum along the third dimension

    # Apply no-data value and write to file
    res_array = res_array.filled(no_datas[0])
    output_path = str(output_folder/output_name)
    write_geotif_sameDomain(files[0], output_path, res_array, nodata_value=no_datas[0])
    
    return output_path

def derive_annual_worldclim_aridity(annual_prec, annual_pet, output_folder, output_name):
    '''Creates an annual aridity map from annual P and PET maps'''

    # Load the data and mask the no-data values
    p_data = get_geotif_data_as_array(annual_prec)
    no_data = get_geotif_nodata_value(annual_prec)
    masked_p = np.ma.masked_array(p_data, mask=p_data==no_data)
    
    pet_data = get_geotif_data_as_array(annual_pet)
    no_data = get_geotif_nodata_value(annual_pet)
    masked_pet = np.ma.masked_array(pet_data, mask=pet_data==no_data)

    # Build in a failsafe for zero P
    masked_p = np.where(masked_p == 0, 1, masked_p)
    
    # Calculate aridity
    ari = np.ma.divide(masked_pet,masked_p)

    # Apply no-data value and write to file
    ari = ari.filled(no_data)
    output_path = str(output_folder/output_name)
    write_geotif_sameDomain(annual_prec, output_path, ari, nodata_value=no_data)
    
    return output_path

def derive_annual_worldclim_fracsnow(annual_prec, annual_snow, output_folder, output_name):
    '''Creates an annual aridity map from annual P and snow maps'''

    # Load the data and mask the no-data values
    p_data = get_geotif_data_as_array(annual_prec)
    no_data = get_geotif_nodata_value(annual_prec)
    masked_p = np.ma.masked_array(p_data, mask=p_data==no_data)
    
    snow_data = get_geotif_data_as_array(annual_snow)
    no_data = get_geotif_nodata_value(annual_snow)
    masked_snow = np.ma.masked_array(snow_data, mask=snow_data==no_data)

    # Build in a failsafe for zero P
    masked_p = np.where(masked_p == 0, 1, masked_p)
    
    # Calculate aridity
    fsnow = np.ma.divide(masked_snow,masked_p)

    # Apply no-data value and write to file
    fsnow = fsnow.filled(no_data)
    output_path = str(output_folder/output_name)
    write_geotif_sameDomain(annual_prec, output_path, fsnow, nodata_value=no_data)
    
    return output_path

def derive_annual_worldclim_seasonality(prec_files, tavg_files, output_folder, output_name):

    '''Calculates a map of seasonality values using WorldClim precipitation and temperature'''

    # Get the precip files
    datasets = [get_geotif_data_as_array(file) for file in prec_files]
    no_datas = [get_geotif_nodata_value(file) for file in prec_files]
    masked_data = [np.ma.masked_array(data, mask=data==no_data) for data,no_data in zip(datasets,no_datas)]
    stacked_prec = np.ma.dstack(masked_data)

    # Get the tavg files
    datasets = [get_geotif_data_as_array(file) for file in tavg_files]
    no_datas = [get_geotif_nodata_value(file) for file in tavg_files]
    masked_data = [np.ma.masked_array(data, mask=data==no_data) for data,no_data in zip(datasets,no_datas)]
    stacked_tavg = np.ma.dstack(masked_data)

    # Fit the sine curves in the third dimension
    t_pars = np.apply_along_axis(fit_temp_sines, axis=2, arr=stacked_tavg)
    p_pars = np.apply_along_axis(fit_prec_sines, axis=2, arr=stacked_prec)

    # Calculate the dimensionless seasonality index in space as per Eq. 14 in Woods (2009)
    seasonality = p_pars[:,:,1]*np.sign(t_pars[:,:,1])*np.cos(2*np.pi*(p_pars[:,:,2]-t_pars[:,:,2])/1)

    # Apply no-data value and write to file
    output_path = str(output_folder/output_name)
    write_geotif_sameDomain(prec_files[0], output_path, seasonality, nodata_value=np.nan)
     
    return output_path