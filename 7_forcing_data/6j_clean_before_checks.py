# Checks if certain old daymet files exist, and removes them if so

import os
from pathlib import Path

# Define the main data folder
main_folder = Path("/project/gwf/gwf_cmt/wknoben/camels_spat/camels-spat-data/basin_data/")

# Find the basin folders in main
basin_folders = [f for f in main_folder.iterdir() if f.is_dir()]

# Sort the basin folders
basin_folders.sort()

# Loop over the basin folders
for basin_folder in basin_folders:
    print(f'Processing {basin_folder.name}')
    
    # Lumped folder
    # Remove any files with the 'Daymet_lumped_remapped_yyyy-01-01-*.nc' pattern
    lumped_folder = basin_folder / "forcing" / "lumped"
    for file in lumped_folder.glob("Daymet_lumped_remapped_*-01-01-*.nc"):
        #print(f'Removing {file}')
        os.remove(file)

    # Rename the 'Daymet_lumped_remapped_daymet_yyyy.nc' files
    for file in lumped_folder.glob("Daymet_lumped_remapped_daymet_*.nc"):
        new_name = file.name.replace("Daymet_lumped_remapped_daymet", "daymet_lumped_remapped")
        #print(f'Renaming {file.name} to {new_name}')
        file.rename(file.with_name(new_name))
    
    # Distributed folder
    dist_folder = basin_folder / "forcing" / "distributed"
    for file in dist_folder.glob("Daymet_dist_remapped_*-01-01-*.nc"):
        #print(f'Removing {file}')
        os.remove(file)
    
    # Rename the 'Daymet_dist_remapped_daymet_yyyy.nc' files
    for file in dist_folder.glob("Daymet_dist_remapped_daymet_*.nc"):
        new_name = file.name.replace("Daymet_dist_remapped_daymet", "daymet_dist_remapped")
        #print(f'Renaming {file.name} to {new_name}')
        file.rename(file.with_name(new_name))
    