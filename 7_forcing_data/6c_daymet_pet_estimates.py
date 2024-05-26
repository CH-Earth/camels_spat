#import dask.array as da
from datetime import datetime
import netCDF4 as nc
import numpy as np
from pathlib import Path
import sys
import time
import xarray as xr
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

# Start the timer
start = time.time()

# --- Command line arguments
if len(sys.argv) > 1:
    ix = int(sys.argv[1])
years = np.arange(1980,2024) # convert SLURM array index to year
year = years[ix]

# --- Make an empty PET file for this year, based on the dayl file
# Data location
daymet_folder = Path("/project/gwf/gwf_cmt/wknoben/camels_spat/temp/daymet/")

# Open the source NetCDF file
src_ds = nc.Dataset(daymet_folder / f'daymet_v4_daily_na_dayl_{year}.nc', 'r')
dst_ds = nc.Dataset(daymet_folder / f'daymet_v4_daily_na_pet_{year}.nc', 'w')

# Copy dimensions from source to destination
for dim_name, dim in src_ds.dimensions.items():
    dst_ds.createDimension(dim_name, (len(dim) if not dim.isunlimited() else None))

# Copy everything but the data variable
for var_name in src_ds.variables:
    var = src_ds.variables[var_name]
    if var_name != 'dayl': 
        dst_var = dst_ds.createVariable(var_name, var.dtype, var.dimensions, fill_value=False)
        dst_var[:] = var[:]
        for attr_name in var.ncattrs():
            dst_var.setncattr(attr_name, var.getncattr(attr_name))

# Create the new variable (empty)
new_var = dst_ds.createVariable('pet', 'f4', ('time', 'y', 'x'), 
                                zlib=True, complevel=4, shuffle=True, fill_value=False)
new_var.long_name = 'potential evapotranspiration estimated with Priestley-Taylor method'
new_var.units = 'mm/day'
new_var.grid_mapping = 'lambert_conformal_conic'

# Set the global attributes
dst_ds.start_year = year
dst_ds.Conventions = 'CF-1.6'
dst_ds.source = 'Priestley-Taylor PET derived from Daymet Data Version 4.0'
dst_ds.history = f'Created for CAMELS-SPAT data set on {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'

# Close the source and destination files
src_ds.close()
dst_ds.close()
print(f'Empty PET file created in {time.time()-start:.2f} seconds')

# --- Existing data loading
elev_ds = xr.open_dataset(daymet_folder / 'daymet_v4_na_elev.nc')
dayl_ds = xr.open_dataset(daymet_folder / f'daymet_v4_daily_na_dayl_{year}.nc')
srad_ds = xr.open_dataset(daymet_folder / f'daymet_v4_daily_na_srad_{year}.nc')
tmax_ds = xr.open_dataset(daymet_folder / f'daymet_v4_daily_na_tmax_{year}.nc')
tmin_ds = xr.open_dataset(daymet_folder / f'daymet_v4_daily_na_tmin_{year}.nc')
vp_ds   = xr.open_dataset(daymet_folder / f'daymet_v4_daily_na_vp_{year}.nc')
print(f'Files opened in {time.time()-start:.2f} seconds')

# --- Calculate the PET estimates
# Loop over the days and append the PET estimates to file
for t in range(0,len(dayl_ds['time'])):
    pet = cs.calculate_priestley_taylor_pet(dayl_ds['yearday'].isel(time=t),
                            dayl_ds['dayl'].isel(time=t),
                            srad_ds['srad'].isel(time=t),
                            tmax_ds['tmax'].isel(time=t),
                            tmin_ds['tmin'].isel(time=t),
                            vp_ds['vp'].isel(time=t),
                            elev_ds['elevation'],
                            elev_ds['lat'])
    # Append to file
    with nc.Dataset(daymet_folder / f'daymet_v4_daily_na_pet_{year}.nc', 'a') as ds:
        ds.variables['pet'][t,:,:] = pet
    print(f'Day {t} processed in {time.time()-start:.2f} seconds')
