import glob
import os
from pathlib import Path

def list_folders(directory):
    folders = []
    directory = str(directory)
    for item in os.listdir(directory):
        if os.path.isdir(os.path.join(directory, item)):
            folders.append(item)
    return sorted(folders)

def list_forcing_files(folder,data):
    search_string = str(folder) + f'/{data}_*.nc'
    files = sorted(glob.glob(search_string))
    return files

def replace_substring_in_file(files,current,new):
    for file in files:
        new_name = file.replace(current,new)
        os.rename(file,new_name)

# main function
def rename_forcing_files(folder,current,new):
    files = list_forcing_files(folder,current)
    replace_substring_in_file(files,current,new)

# processing
main_folder = Path('/project/gwf/gwf_cmt/wknoben/camels_spat/camels-spat-data/basin_data')
sub_folders = list_folders(main_folder)
current_name = 'EM-Earth' # with '-'
new_name = 'EM_Earth' # with '_'
for sub_folder in sub_folders:
    for spatial_res in ['lumped', 'distributed']: # We know files in 'raw' are as 'EM_Earth'
        data_folder = main_folder / sub_folder / 'forcing' / spatial_res
        rename_forcing_files(data_folder, current_name, new_name)

