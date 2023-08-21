# Find shapefile bounding box
def find_shapefile_bounds(path):
    # Modified from: https://github.com/CH-Earth/CWARHM/blob/main/0_tools/ERA5_find_download_coordinates_from_shapefile.ipynb
    shp = gpd.read_file(path)
    return shp.total_bounds

