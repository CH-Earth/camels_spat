# Fixes minor processing errors in several shapefiles
import sys
import geopandas as gpd
import pandas as pd
from pathlib import Path
from shapely.ops import unary_union
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

# --- CONFIG HANDLING
# Specify where the config file can be found
config_file = '../0_config/config.txt'

# Get the required info from the config file
data_path     = cs.read_from_config(config_file,'data_path')

# CAMELS-spat metadata
cs_meta_path  = cs.read_from_config(config_file,'cs_basin_path')
cs_meta_name  = cs.read_from_config(config_file,'cs_meta_name')

# Basin folder
cs_basin_folder = cs.read_from_config(config_file, 'cs_basin_path')

# --- DATA LOADING
# CAMELS-spat metadata file
cs_meta_path = Path(data_path) / cs_meta_path
cs_meta = pd.read_csv(cs_meta_path / cs_meta_name)

# --- LIST FIXES
issue_dict = {'CAN_02GG003': ['hru_area'], 
              'CAN_04FC001': ['hru_area', 'double_comid'], 
              'CAN_07AD002': ['hru_area'], 
              'CAN_07HC002': ['hru_area'], 
              'CAN_08NE077': ['hru_area'], 
              'CAN_08NH007': ['hru_area'], 
              'USA_07142300':['hru_area'], 
              'USA_08198500':['hru_area'],
              'CAN_06AB001': ['double_comid'],
              'CAN_06BB005': ['double_comid'],
              'CAN_08MH076': ['double_comid']}

for i,row in cs_meta.iterrows():
    basin = row['Country'] + '_' + row['Station_id']
    if basin in issue_dict.keys():
        issues = issue_dict[basin]
        print(f'\nFixing {issues} in {basin}')
        gauge_id, shp_lump_path, shp_dist_path, _, _ = cs.prepare_delineation_outputs(cs_meta, i, Path(data_path)/cs_basin_folder)
        assert gauge_id == basin, f'Gauge ID {gauge_id} does not match {basin}'
        for issue in issues:            
            if issue == 'hru_area':
                # Need to load these in here for the basin that needs both fixes
                shp_bas = gpd.read_file(str(shp_dist_path).format('basin'))
                shp_riv = gpd.read_file(str(shp_dist_path).format('river'))
                # Check for HRUs with area below 0.0005 km2
                assert any(shp_bas['unitarea'] < 0.0005), 'No HRUs with area below 0.0005 km2 found'
                assert any(shp_bas['unitarea'] / shp_bas['unitarea'].sum() < 0.001), 'No HRUs with area below 0.1% of total area found'
                mask = shp_bas['unitarea'] < 0.0005
                print(f'HRUs with area below 0.0005 km2:\n {shp_bas.loc[mask]}')
                # Drop the small HRU from both files
                comid_to_remove = shp_bas.loc[mask]['COMID'].values # We know through manual checks that in all cases here this is only 1 COMID
                new_bas = shp_bas[~mask].copy().reset_index(drop=True)
                new_riv = shp_riv[~shp_riv['COMID'].isin(comid_to_remove)].copy().reset_index(drop=True)
                # Check that we lost what we think we should lose
                riv_count = shp_riv['COMID'].isin(comid_to_remove).sum()
                assert len(new_bas) == len(shp_bas) - len(comid_to_remove), f'Expected # HRUs removed = {len(comid_to_remove)}, but old # HRUs = {len(shp_bas)}, and new # HRUs = {len(new_bas)}'
                assert len(new_riv) == len(shp_riv) - riv_count, f'Expected # HRUs removed = {riv_count}, but old # HRUs = {len(shp_riv)}, and new # HRUs = {len(new_riv)}'
                # Save the new files
                new_bas.to_file(str(shp_dist_path).format('basin'))
                new_riv.to_file(str(shp_dist_path).format('river'))
            if issue == 'double_comid':
                # Need to load this in here for the basin that needs both fixes
                shp_riv = gpd.read_file(str(shp_dist_path).format('river'))
                # Check for duplicate COMIDs
                assert any(shp_riv['COMID'].duplicated()), 'No duplicate COMIDs found'
                dupl_mask = shp_riv['COMID'].duplicated()
                print(f'Duplicate COMIDs found:\n {shp_riv.loc[dupl_mask]}')
                # Get all the unique duplicates; we need to merge these on a per-COMID basis
                comids_to_merge = pd.unique(shp_riv.loc[dupl_mask]['COMID'].values) # Get all the unique duplicates; we need to merge these on a per-COMID basis
                new_riv = shp_riv[~dupl_mask].copy().reset_index(drop=True) # retains one of each duplicate
                # Loop over the COMIDs that need to be merged
                for comid in comids_to_merge:
                    print(f'Merging COMID {comid}')
                    merge_mask = shp_riv['COMID'] == comid
                    new_geom = unary_union(shp_riv.loc[merge_mask]['geometry'])
                    update_mask = new_riv['COMID'] == comid
                    new_riv.loc[update_mask,'geometry'] = new_geom
                # Check that we lost what we think we should lose
                assert len(new_riv) == len(shp_riv) - dupl_mask.sum(), f'Expected # COMIDs removed = {dupl_mask.sum()}, but old # COMIDs = {len(shp_riv)}, and new # COMIDs = {len(new_riv)}'
                new_riv.to_file(str(shp_dist_path).format('river'))
