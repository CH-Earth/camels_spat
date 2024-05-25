import os
import numpy as np
import rasterio
import xarray as xarray
from pathlib import Path
from osgeo import gdal

# Data location
daymet_folder = Path("/project/gwf/gwf_cmt/wknoben/camels_spat/temp/daymet/")
dayl_file = 'daymet_v4_daily_na_dayl_1980.nc'
elev_file = 'daymet_v4_na_elev.nc'
dem_folder = Path('/project/gwf/gwf_cmt/wknoben/camels_spat/geospatial_temp/merit/raw/')
dem_file = 'merit_hydro_elv.tif'
lcc_file = 'merit_hydro_elv_lcc.tif'

# --- Prepare the DEM in LCC

# Open the NetCDF file
ds = xarray.open_dataset(daymet_folder / dayl_file)

# Extract the projection information
proj_attr = ds['lambert_conformal_conic'].attrs

# Construct the PROJ string
proj_string = (
    f"+proj=lcc "
    f"+lon_0={proj_attr['longitude_of_central_meridian']} "
    f"+lat_0={proj_attr['latitude_of_projection_origin']} "
    f"+x_0={proj_attr['false_easting']} "
    f"+y_0={proj_attr['false_northing']} "
    f"+lat_1={proj_attr['standard_parallel'][0]} "
    f"+lat_2={proj_attr['standard_parallel'][1]} "
    f"+a={proj_attr['semi_major_axis']} "
    f"+rf={proj_attr['inverse_flattening']}"
)

# Open the input dataset
src_ds = gdal.Open( str(dem_folder/dem_file) )

# Perform the warp using the PROJ string
gdal.Warp(str(dem_folder/lcc_file), src_ds, dstSRS=proj_string)

# Clean up
src_ds = None

# --- Extract elevations for the x,y coordinates in the Daymet file

# Open the GeoTIFF file
with rasterio.open(str(dem_folder/lcc_file)) as src:
    # Read the elevation data
    elevation_data = src.read(1)
    # Get the affine transformation of the GeoTIFF
    transform = src.transform
    # Get the coordinate reference system (CRS) of the GeoTIFF
    tif_crs = src.crs

# Create a transformer object
transformer = rasterio.transform.AffineTransformer(transform)

# Extract the x and y coordinates from the NetCDF file
x_coords = ds['x'].values
y_coords = ds['y'].values

# Create a meshgrid of coordinates if necessary
x_mesh, y_mesh = np.meshgrid(x_coords, y_coords)

# Flatten the meshgrid to have a list of (x, y) pairs
x_flat = x_mesh.flatten()
y_flat = y_mesh.flatten()

# Extract elevation values for each (x, y) pair
elevation_values = []
for x, y in zip(x_flat, y_flat):
    row, col = transformer.rowcol(x, y)
    if 0 <= row < elevation_data.shape[0] and 0 <= col < elevation_data.shape[1]:
        elevation = elevation_data[row, col]
        elevation_values.append(elevation)
    else:
        elevation_values.append(np.nan)  # Handle out-of-bounds coordinates

# Reshape the elevation values to match the shape of the original coordinate meshgrid
elevation_values = np.array(elevation_values).reshape(x_mesh.shape)

# --- Save to file
# Extract the x, y, lat, and lon coordinates and their attributes
x_coords = ds['x'].values
y_coords = ds['y'].values
lat_coords = ds['lat'].values
lon_coords = ds['lon'].values

x_attrs = ds['x'].attrs
y_attrs = ds['y'].attrs
lat_attrs = ds['lat'].attrs
lon_attrs = ds['lon'].attrs

# Assuming elevation_values is your data matrix corresponding to these coordinates
# Create a DataArray with elevation values and copy the coordinate attributes
elevation_da = xr.DataArray(
    elevation_values,
    dims=["y", "x"],
    coords={
        "x": ("x", x_coords, x_attrs),
        "y": ("y", y_coords, y_attrs),
        "lat": (["y", "x"], lat_coords, lat_attrs),
        "lon": (["y", "x"], lon_coords, lon_attrs)
    },
    attrs={
        "units": "m.a.s.l.",
        "long_name": "Elevation"
    }
)

# Create a Dataset and add the DataArray
elevation_ds = xr.Dataset({"elevation": elevation_da, 
                           'lambert_conformal_conic': ds['lambert_conformal_conic']})

# Add global attributes if needed
elevation_ds.attrs = {
    "title": "Elevation Data",
    "source": "Extracted from MERIT Hydro DEM, with coordinates from Daymet v4.1 Lambert Conformal Conic projection",
    "Conventions": "CF-1.6"
}

# To file
elevation_ds.to_netcdf(daymet_folder / elev_file)

# --- Make a plot to check
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# file name
elev_fig = daymet_folder / 'elevation_map.png'

# Extract the elevation data and coordinates
elevation = elevation_ds['elevation']

# Create a plot
fig, ax = plt.subplots()

# Plot the elevation data
elevation.plot(ax=ax, cmap='terrain', vmin=-100, vmax=5000,
               cbar_kwargs={'label': 'Elevation (meters)'})

# Set the title and labels
ax.set_title('Elevation Plot')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

# Save the plot to a file
plt.savefig(elev_fig, bbox_inches='tight')

# Close the plot to free up memory
plt.close(fig)

# --- Remove the LCC dem again, to save space
os.remove(dem_folder/lcc_file)