# Takes the individual CSV files we have and merges into a single one
# We can use this to check and rerun forcing creation for basins
import glob
import pandas as pd
import xarray as xr
from pathlib import Path

# Define unusable location
unusable_path = Path('/project/gwf/gwf_cmt/wknoben/camels_spat/camels-spat-data/camels_spat_unusable.csv')
cs_unusable = pd.read_csv(unusable_path)

# Define meta location
meta_path = Path('/project/gwf/gwf_cmt/wknoben/camels_spat/camels-spat-data/camels_spat_metadata.csv')
cs_meta = pd.read_csv(meta_path)

# Define csv location
rerun_path = Path('/project/gwf/gwf_cmt/wknoben/camels_spat/camels-spat-data/forcing_reruns')
csv_files = sorted(glob.glob(str(rerun_path) + '/*.csv'))

# Read in all CSV files
dfs = []
for csv_file in csv_files:
    df = pd.read_csv(csv_file)
    dfs.append(df)

# Concatenate all dataframes
df = pd.concat(dfs, ignore_index=True)

# --- ANALYSIS ---
# Check unique number of basins
basins = df['basin'].unique()
print(f'Unique basins: {len(basins)}') # 1426

# Check how often each reason appears
reasons = df['reason'].value_counts()
print(reasons)

'''
reason
NaNs in variables in raw EM_Earth                1374 # Most of these are cases where the basin is close to open water and EM-Earth contains a few NaN values for open water cells
NaNs in variables in lumped EM_Earth             1372 # These are the result of basins having 1950-01-01 and/or 1979-01-01 as part of there time domain. These files contain NaNs in the first few timesteps
NaNs in variables in distributed EM_Earth        1372 #   This is because EM-Earth relies on ERA5, which initially only covered the time period 1979-01-01 to current, and later extended back to 1950-01-01.
Time zone is NST                                   57 # These are known
Missing invariant file                             13 # These are all lake stations that we exclude from the dataset
Basin had tiny HRUs removed                         8 # These are known: ['CAN_02GG003', 'CAN_04FC001', 'CAN_07AD002', 'CAN_07HC002', 'CAN_08NE077', 'CAN_08NH007', 'USA_07142300', 'USA_08198500']
Non-consecutive dates in raw EM_Earth               7 # v-- These are the same basins --v
Non-consecutive dates in distributed ERA5           7 # ['USA_01638480', 'USA_02384540',
Non-consecutive dates in raw ERA5                   7 #  'USA_02479300', 'USA_04224775',
Non-consecutive dates in lumped EM_Earth            7 #  'USA_08171300', 'USA_08198500',
Non-consecutive dates in lumped ERA5                7 #  'USA_12375900']
Non-consecutive dates in distributed EM_Earth       7 # ^-- These are the same basins --^
Missing variables in distributed ERA5               2 # Same as below: ['CAN_09ED001', 'CAN_10ED002']
NaNs in variables in distributed ERA5               2 # Same as above: ['CAN_09ED001', 'CAN_10ED002']
Missing variables in lumped EM_Earth                1 # ['CAN_10JC003']
'''

# Check which basin is missing variables in lumped EM_Earth
missing_variables = df[df['reason'].str.contains('Missing variables in lumped EM_Earth')]
print(missing_variables['basin']) # ['CAN_10JC003']

# Cech which basins are missing variables in distributed ERA5
missing_variables = df[df['reason'].str.contains('Missing variables in distributed ERA5')]
print(missing_variables['basin']) # ['CAN_09ED001' 'CAN_10ED002']

# Check which basins have NaNs in variables in distributed ERA5
nan_variables = df[df['reason'].str.contains('NaNs in variables in distributed ERA5')]
print(nan_variables['basin']) # ['CAN_09ED001' 'CAN_10ED002']

# Check if all the 'non-consecutive dates' are in the same basins
non_consecutive_dates = df[df['reason'].str.contains('Non-consecutive dates')]
print(non_consecutive_dates['basin'].unique()) # yes: ['USA_01638480' 'USA_02384540' 'USA_02479300' 'USA_04224775' 'USA_08171300' 'USA_08198500' 'USA_12375900']

# Check which basins had tiny HRUs removed
tiny_hrus = df[df['reason'].str.contains('Basin had tiny HRUs removed')]
print(tiny_hrus['basin']) # ['CAN_02GG003', 'CAN_04FC001', 'CAN_07AD002', 'CAN_07HC002', 'CAN_08NE077', 'CAN_08NH007', 'USA_07142300', 'USA_08198500']

# Check which basins are missing invariant files
missing_invariant = df[df['reason'].str.contains('Missing invariant file')]
print(missing_invariant['basin'])

count = 0
for basin in missing_invariant['basin']:
    country = basin[0:3]
    station = basin[4:]
    unusable_mask = (cs_unusable['Country'] == country) & (cs_unusable['Station_id'] == station)
    if len(cs_unusable.loc[unusable_mask]) > 0:
        count += 1
        print(f'Basin {basin} is in the unusable list. Reason: {cs_unusable.loc[unusable_mask,"Reason"].values[0]}')

if count == len(missing_invariant):
    print('All basins are in the unusable list')

# --- MANUAL CHECKS ---
# CAN_10JC003: Missing variables in lumped EM_Earth
# Looks like an interruption to remapping script, where EM-EARTH prcp and tmean 
# are missing from a single file: EM_Earth_lumped_remapped_2005-09-01-00-00-00.nc
print(cs_meta[(cs_meta['Country'] == 'CAN') & (cs_meta['Station_id'] == '10JC003')].index)

# CAN_09ED001: Missing variables and NaNs in distributed ERA5
# Looks like interruptions to processing scripts.
# Missing phi (wind direction) in 1 file (u and v are present): ERA5_dist_remapped_1983-01-01-00-00-00.nc
# Missing 13 variables and NaNs in 1 file: ERA5_dist_remapped_1990-09-01-00-00-00.nc
print(cs_meta[(cs_meta['Country'] == 'CAN') & (cs_meta['Station_id'] == '09ED001')].index)

# CAN_10ED002: Missing variables and NaNs in distributed ERA5
# Looks like interruptions to processing scripts.
# Missing phi in 1 file (u and v are present): ERA5_dist_remapped_1972-05-01-00-00-00.nc
# Missing 3 variables and Nans in 2 variables in 1 file: ERA5_dist_remapped_1981-09-01-00-00-00.nc
print(cs_meta[(cs_meta['Country'] == 'CAN') & (cs_meta['Station_id'] == '10ED002')].index)

# -- NaNs in EM-EARTH
# There appear to be two cases:
# 1. NaNs are present in (almost?) every raw file for the basin
#  - This seems to never occur for lumped and distributed, where we first subset the raw data to the basin outlines
#
#                 raw  lumped  distributed
# basin
# CAN_01AD002   1681       2            2 --> close to open water, where EM-EARTH has NaN data
# USA_11148900  1161       1            1 --> close to open water, where EM-EARTH has NaN data
# CAN_02PG022    903       0            0 --> close to open water, where EM-EARTH has NaN data
# CAN_02QA002   1227       1            1 --> close to open water, where EM-EARTH has NaN data
# CAN_03QC001   1279       1            1 --> close to open water, where EM-EARTH has NaN data
# CAN_08MG013   1649       1            1 --> close to open water, where EM-EARTH has NaN data
# CAN_09CD001   1521       1            1 --> close to open water, where EM-EARTH has NaN data
# CAN_09ED001    841       0            0 --> close to open water, where EM-EARTH has NaN data
# 
# 2. NaNs are present in a few files for the basin:
#  - There are two sub-cases here:
#    - NaNs are present in 1 raw file, 1 lumped file, and 1 distributed file. This is always 1979-01
#    - NaNs are present in 2 raw files, 1 lumped file, and 1 distributed file. This is always 1950-01
#  - This only happens with variable prcp
#
# Using CAN_10ED002 as an example:
# raw file has NaN for first 7 timesteps, and this propages to the lumped and distributed files:
# 
# em_lump_1['prcp'].head(10)
# <xarray.DataArray 'prcp' (time: 10, hru: 1)> Size: 40B
# array([[         nan],
#        [         nan],
#        [         nan],
#        [         nan],
#        [         nan],
#        [         nan],
#        [         nan],
#        [2.276279e-08],
#        [1.975637e-08],
#        [1.633712e-08]], dtype=float32)
# Coordinates:
#   * time     (time) datetime64[ns] 80B 1978-12-31T17:00:00 ... 1979-01-01T02:...
#
# This originates in the source EM-Earth files.
# A plausible explanation is that EM_Earth relies on ERA5 reanalysis. While ERA5 was in production,
# it initially covered the time period 1979-01-01 to current, and later extended back to 1950-01-01.
# I noticed this earlier in ERA5 precipitation data and asked their support (2020-03-20):
#
# The reason is that the accumulated variables for ERA5 come from the forecast runs of the model which 
# are made at 0600 and 1800 each day. The hourly forecast data from  1800 are used to provide data 
# between 1900 ann 0600 the following day. So for 1/1/1979 0000-0600, these data would come from the 
# 1800 forecast run on 31/12/1978. As these 1978 data have not yet been released, there is no accumulated 
# variable data for 0000-0600 on 1/1/1979. Analysis variables such as sp are present for this period, 
# which explains the difference. 

# Define the log file location
log_path = Path('/globalhome/wmk934/HPC/camels_spat/7_forcing_data/forcing_check_logs/')
log_files = sorted(glob.glob(str(log_path) + '/checks_screen_*.out'))

# Count per log file (i.e. basin), how often we encounter the EM-Earth variable is NaN issue
basin_list = []
nan_counts_raw = []
nan_counts_lumped = []
nan_counts_distributed = []
for log_file in log_files:
    with open(log_file, 'r') as f:
        lines = f.readlines()
    nan_count_raw = 0
    nan_count_lumped = 0
    nan_count_distributed = 0
    for line in lines:
        if 'Checking basin forcing data for' in line:
            basin_list.append(line.split(' ')[-1].strip())
        if 'NaNs in variable' in line and 'EM_Earth' in line and not 'Reasons' in line:
            if '/raw/' in line:
                nan_count_raw += 1
            elif '/lumped/' in line:
                nan_count_lumped += 1
            elif '/distributed/' in line:
                nan_count_distributed += 1
    nan_counts_raw.append(nan_count_raw)
    nan_counts_lumped.append(nan_count_lumped)
    nan_counts_distributed.append(nan_count_distributed)

# Store everything in a dataframe
em_earth_nan_counts = pd.DataFrame({'basin': basin_list, 
                                    'raw': nan_counts_raw, 
                                    'lumped': nan_counts_lumped, 
                                    'distributed': nan_counts_distributed})
em_earth_nan_counts.set_index('basin', inplace=True)

# Check how many basins have NaNs in (almost) all raw files
print(len(em_earth_nan_counts[em_earth_nan_counts['raw'] > 2])) # >0: 1374; >1: 363; >2: 8; >3-499: 8
    
# Get the 8 basins where we have a lot of files with NaNs
#  Added to the comments above
print(em_earth_nan_counts[em_earth_nan_counts['raw'] > 2])

# Check if the cases where we have NaN in 1 or raw files
# Hypothesis: this NaN value is also always present in the lumped and distributed files
mask = (em_earth_nan_counts['raw'] < 3) & (em_earth_nan_counts['lumped'] > 0) & (em_earth_nan_counts['distributed'] > 0)
raw_vals = em_earth_nan_counts[mask]['raw']
lump_vals = em_earth_nan_counts[mask]['lumped']
dist_vals = em_earth_nan_counts[mask]['distributed']
assert (raw_vals == lump_vals).all()  # If this happens, we always have NaNs in 2 raw files
assert (raw_vals == dist_vals).all() # This always correspond with NaNs in 1 lumped file

# Check the filenames of the files where we have NaNs in the lumped and distributed files
raw_nan_files = []
lump_nan_files = []
dist_nan_files = []
for log_file in log_files:
    with open(log_file, 'r') as f:
        lines = f.readlines()
    basin = lines[0].split(' ')[-1].strip()
    if basin in em_earth_nan_counts[mask].index:
        for line in lines:
            if 'NaNs in variable' in line and 'EM_Earth' in line:
                file = line.split('/')[-1].strip()
                var = line.split(' ')[4].strip()
                if '/lumped/' in line:
                    lump_nan_files.append((file, var))
                elif '/distributed/' in line:
                    dist_nan_files.append((file, var))
                elif '/raw/' in line:
                    raw_nan_files.append((file, var))

print(f'unique files with NaNs in raw: {set(raw_nan_files)}')
print(f'unique files with NaNs in lumped: {set(lump_nan_files)}')
print(f'unique files with NaNs in distributed: {set(dist_nan_files)}')

# Check a few random files to try and figure out what the problem is
em_raw_file_1 = Path('/project/gwf/gwf_cmt/wknoben/camels_spat/camels-spat-data/basin_data/CAN_10ED002/forcing/raw/EM_Earth_1979-01.nc')
em_lump_file_1 = Path('/project/gwf/gwf_cmt/wknoben/camels_spat/camels-spat-data/basin_data/CAN_10ED002/forcing/lumped/EM_Earth_lumped_remapped_1979-01-01-00-00-00.nc')
em_dist_file_1 = Path('/project/gwf/gwf_cmt/wknoben/camels_spat/camels-spat-data/basin_data/CAN_10ED002/forcing/distributed/EM_Earth_dist_remapped_1979-01-01-00-00-00.nc')
em_raw_1 = xr.open_dataset(em_raw_file_1) 
em_lump_1 = xr.open_dataset(em_lump_file_1) # first 7 timesteps are all NaN
em_dist_1 = xr.open_dataset(em_dist_file_1)

# EM-EARTH checks of 1950-01 and 1979-01 files
source_195001 = Path('/project/gwf/gwf_cmt/wknoben/camels_spat/temp_em_earth/EM_Earth_v1/deterministic_hourly/prcp/NorthAmerica/EM_Earth_deterministic_hourly_NorthAmerica_195001.nc')
source_197901 = Path('/project/gwf/gwf_cmt/wknoben/camels_spat/temp_em_earth/EM_Earth_v1/deterministic_hourly/prcp/NorthAmerica/EM_Earth_deterministic_hourly_NorthAmerica_197901.nc')
native_195001 = xr.open_dataset(source_195001)
native_197901 = xr.open_dataset(source_197901)

# check if lat and long are nan on time=0
print(f'file: {source_195001}')
for time in range(0,10):
    all_nan = native_195001['prcp'].isel(time=0).isnull().all().values
    print(f'time={time}: {all_nan}')

print(f'file: {source_197901}')
for time in range(0,10):
    all_nan = native_197901['prcp'].isel(time=0).isnull().all().values
    print(f'time={time}: {all_nan}')

# --- DEFINE THE FINAL FORCING RERUN LIST
# All reasons other than the NaN in EM-Earth variables should be re-run
reasons_to_keep = ['Time zone is NST', 
                   'Basin had tiny HRUs removed', 
                   'Non-consecutive dates in raw EM_Earth', 
                   'Non-consecutive dates in raw ERA5',
                   'Non-consecutive dates in lumped EM_Earth', 
                   'Non-consecutive dates in lumped ERA5', 
                   'Non-consecutive dates in distributed EM_Earth',
                   'Non-consecutive dates in distributed ERA5', 
                   'Missing variables in distributed ERA5', 
                   'NaNs in variables in distributed ERA5', 
                   'Missing variables in lumped EM_Earth']
final_reruns = df[df['reason'].isin(reasons_to_keep)]['basin'].unique()
print(len(final_reruns)) # 74

# Save to file
rerun_file = Path('/globalhome/wmk934/HPC/camels_spat/7_forcing_data/forcing_check_logs/reruns_20240516.csv')
final_reruns_df = pd.DataFrame({'basin': final_reruns})
final_reruns_df.to_csv(rerun_file, index=False)