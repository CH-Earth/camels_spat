import glob
import os
import pandas as pd
import sys
from pathlib import Path
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

# --- Config handling
# Specify where the config file can be found
config_file = '../0_config/config.txt'

# Get the required info from the config file
data_path = cs.read_from_config(config_file,'data_path')

# Basin folder
cs_basin_folder = cs.read_from_config(config_file, 'cs_basin_path')
basins_path = Path(data_path) / cs_basin_folder

# Get the attribute folder
att_folder = cs.read_from_config(config_file, 'att_path')
att_path = basins_path / att_folder

# -- Find files
all_files = sorted(glob.glob(str(att_path / 'attributes*.csv')))

# -- Make one
def open_attribute_csv(file,update_column_name=False):
    '''Opens an attribute CSV and creates a multi-index for it'''
    df = pd.read_csv(file)
    mi = pd.MultiIndex.from_frame(df[['Category','Attribute','Unit','Source']])
    df.index = mi
    df = df.drop(['Category','Attribute','Unit','Source'], axis=1)
    if update_column_name:
        basin_id = os.path.basename(file).split('.')[0].replace('attributes_','')
        df = df.rename(columns={df.columns[0]: basin_id})
    return df

# Iterate over each file
dfs = []
for filename in all_files:
    df = open_attribute_csv(filename, update_column_name=True)
    dfs.append(df)

# Concatenate all DataFrames in the list into a single DataFrame
combined_df = pd.concat(dfs, axis=1)
combined_df.to_csv(Path(data_path) / 'camels-spat-data' / 'camels_spat_attributes.csv')