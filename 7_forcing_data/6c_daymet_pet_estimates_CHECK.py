import gc
import glob
import matplotlib.pyplot as plt
import xarray as xr
from pathlib import Path

# Data location
daymet_folder = Path("/project/gwf/gwf_cmt/wknoben/camels_spat/temp/daymet/")
pet_files = sorted(glob.glob(str(daymet_folder / 'daymet_v4_daily_na_pet_*.nc')))

# 1. Check if we have different values in different years
# Loop over the files and check if Jan-1 and Jul-1 are different
ds_1980 = xr.open_dataset(pet_files[0])
jan1_1980 = ds_1980['pet'].isel(time=0)
jul1_1980 = ds_1980['pet'].isel(time=181)

for pet_file in pet_files[1:]:
    # Open the file
    ds = xr.open_dataset(pet_file)
    
    # Check the values
    jan1 = ds['pet'].isel(time=0)
    jul1 = ds['pet'].isel(time=181)
    if (jan1 == jan1_1980).all():
        print(f'PET values are the same as Jan-1-1980 for Jan-1 in {pet_file}')
    else:
        print(f'PET values are different for Jan-1-1980 and Jan-1 in {pet_file}')
    if (jul1 == jul1_1980).all():
        print(f'PET values are the same as Jul-1-1980 for Jul-1 in {pet_file}')
    else:
        print(f'PET values are different for Jul-1-1980 and Jul-1 in {pet_file}')
    
    # Close the file
    ds.close()

ds_1980.close()

# 2. Create images for random days in each season and year, to see if things look reasonable
images_folder = Path("/project/gwf/gwf_cmt/wknoben/camels_spat/temp/daymet/images_pet_checks/")
images_folder.mkdir(exist_ok=True)

# Loop over the files and create images for random days
for pet_file in pet_files:
    # Open the file
    ds = xr.open_dataset(pet_file)
    
    # Create images for random days
    for t in [0, 90, 181, 273]:
        pet = ds['pet'].isel(time=t)
        pet.plot.imshow()
        plt.title(f'PET on {ds["time"].isel(time=t).values}')
        plt.savefig(images_folder / f'pet_{ds["time"].isel(time=t).values}.png')
        plt.close()
    
    # Close the file
    ds.close()
    gc.collect()