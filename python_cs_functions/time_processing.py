'''Contains functions for temporal processing of downloaded data (timezones, averaging, etc.)'''

import netCDF4 as nc4
import numpy as np
import pandas as pd
from scipy import integrate
import time
from datetime import timedelta
from datetime import datetime

# Time zone conversion
# ---------------------------------------
def tz_abbreviation_to_utc(tz_str):
    '''Takes a timezone abbreviation used by USGS and converts to UTC-xx'''
    
    # Source: https://nrc.canada.ca/en/certifications-evaluations-standards/canadas-official-time/time-zones-daylight-saving-time
    # Source: https://www.timeanddate.com/time/zones/akst
    
    # Check inputs
    assert type(tz_str) == str, f'Input tz_str must be {str} but is {type(tz_str)}'
    
    # Define abbreviation-to-UTC dictionary
    dict_abbr_to_utc = {'AST': 'UTC-04',    #            Atlantic Standard Time. 
                        'CDT': 'UTC-05',    # Ambiguous: Central Daylight Time.  Possibly legacy part of pytz: CST6CDT
                        'CST': 'UTC-06',    # Ambiguous: Central Standard Time.  Not in pytz
                        'EDT': 'UTC-04',    #            Eastern Daylight Time.  legacy part of pytz:          EST5EDT
                        'EST': 'UTC-05',    #            Eastern Standard Time.  Part of pytz:                 EST
                        'LST': 'UTC-09',    #            Official abbreviation AKST: Alaska Standard Time
                        'MDT': 'UTC-06',    #            Mountain Daylight Time. legacy part of pytz:          MST7MDT
                        'MST': 'UTC-07',    # Ambiguous: Mountain Standard Time. Possibly part of pytz:        MST
                        'NST': 'UTC-0330',  #            Newfoundland Standard Time.
                        'PDT': 'UTC-07',    #            Pacific Daylight Time.  Legacy part of pytz:          PST8PDT
                        'PST': 'UTC-08'}    #            Pacific Standard Time.  Not in pytz   
    
    return dict_abbr_to_utc[tz_str]

def relative_utc_to_timedelta(utc_str):
    '''Takes a UTC-xx string and converts to Timedelta string'''
    
    # Check inputs
    assert type(utc_str) == str, f'Input tz_str must be {str} but is {type(utc_str)}'
    
    # Define UTC-to-timedelta dictionary
    dict_utc_to_timedelta = {'UTC-0330': '-3:30:00',
                             'UTC-04'  : '-4:00:00',
                             'UTC-05'  : '-5:00:00',
                             'UTC-06'  : '-6:00:00',
                             'UTC-07'  : '-7:00:00',
                             'UTC-08'  : '-8:00:00',
                             'UTC-09'  : '-9:00:00'}
    
    return dict_utc_to_timedelta[utc_str]

def relative_utc_to_float_offset_in_hours(utc_str):
    '''Takes a UTC-xx string and converts to a float in hours'''
    
    # Check inputs
    assert type(utc_str) == str, f'Input tz_str must be {str} but is {type(utc_str)}'
    
    # Define UTC-to-timedelta dictionary
    dict_utc_to_float = {'UTC-0330': -3.5,
                         'UTC-04'  : -4.0,
                         'UTC-05'  : -5.0,
                         'UTC-06'  : -6.0,
                         'UTC-07'  : -7.0,
                         'UTC-08'  : -8.0,
                         'UTC-09'  : -9.0}
    
    return dict_utc_to_float[utc_str]

def datetime_str_to_timeaware_datetime(date_str, offset='+00:00:00', localize_to_UTC=True):
    
    '''Takes a datetime string (yyyy-mm-dd hh:mm) and possible an offset 
    (+/-hh:mm:ss) to generate a timezone-aware datetime in UTC'''
    
    # Check inputs
    assert type(date_str) == str, f'Input tz_str must be {str} but is {type(date_str)}'
    assert type(offset) == str, f'Input tz_str must be {str} but is {type(offset)}'
    
    # Convert the time
    new_time = pd.to_datetime(date_str) + pd.Timedelta(offset)
    if localize_to_UTC: new_time = new_time.tz_localize('UTC')
    
    return new_time

# Resampling to hourly
# --------------------
# Fill gaps with temporary values where needed
def fill_data_gaps_on_integration_bounds(df, data='data', center_window=True):
    
    for ix,(time,row) in enumerate(df.iterrows()):

        # Missing value at a given point in time
        if np.isnan(row[data]):

            # Only proceed if this is a timestamp that we'll need as a boundary value for volume integration
            if ((center_window) & (time.minute == 30)) or \
            ((not center_window)& (time.minute == 0)):

                # Check if we have any valid values in the past and/or next hour
                prev_time_mask = slice(time - timedelta(hours=1),time) # Returns empty series if time slice is out of index range
                prev_data_mask = (df[data].loc[prev_time_mask] >= 0) & \
                                 (df['is_obs'].loc[prev_time_mask]) # Has data and is observation (i.e. not previously filled value)
                next_time_mask = slice(time, time + timedelta(hours=1))
                next_data_mask = (df[data].loc[next_time_mask] >= 0) & \
                                 (df['is_obs'].loc[next_time_mask])

                # Select the correct case
                if any(prev_data_mask) and any(next_data_mask):
                    # do nothing: we'll interpolate the missing value later
                    continue # to next row in the dataframe
                elif any(prev_data_mask):
                    # Copy closest value from the previous hour into this position
                    df.at[time,data] = df[prev_time_mask][prev_data_mask].iloc[-1][data]
                elif any(next_data_mask):
                    # Copy the closest value from the next hour into this position
                    df.at[time,data] = df[next_time_mask][next_data_mask].iloc[-1][data]
                else:
                    # do nothing: we don't have observations from either the previous or the next hour, so we'll leave this empty
                    continue # to next row in the dataframe
    
    return df

def resample_arbitrary_flux_observations_to_hourly(df, data='data', center_window=True, keep_tmp_col=False):
    
    '''Takes a dataframe with a timestamp index and a "data" column and returns old and new dataframe with hourly data values'''
    
    # Inputs:
    # df            - dataframe with a column with name [data]
    # data          - column name; default 'data'
    # center_window - flag to state if hourly resampling should take the whole hour as the mid-point of the window
    #                 (e.g. if the value for 12:00 should be based on [True] 11:30-12:30; or [False] on 12:00-13:00)
    # keep_tmp_col  - flag to keep temporary columns (useful for debugging; default False, i.e. don't keep)
    
    # Copy the input dataframe so we don't change it outside this function
    df = df.copy()
    
    # Ensure we have values at every half hour, so that resampling works properly
    half_hour_start = df.index[0].ceil('30min') # round up to the nearest half hour, so that we can safely interpolate
    half_hour_end   = df.index[-1].floor('30min') # round down to nearest half hour so we only interpolate between observations
    half_hour_times = pd.date_range(half_hour_start,half_hour_end, freq='30min') # Adding 30 minute intervals means this works 
                                                                                 #    with both center_window settings
    
    # Insert the new times into the old dataframe
    df['is_obs'] = 1 # Add a new column that flags the original observations, so we know what was interpolated later on
    tmp_df = pd.DataFrame(index=half_hour_times) # automatically puts in NaNs for data
    df = pd.concat([df,tmp_df], axis=1) # automatically gives us NaNs in places where we don't have 'data' values
    
    # Attempt to fill gaps on the left-hand and right-hand side of the integration windows, where gaps longer than 1h exist
    #  between actual observed values. We don't want to use interpolated values here because they might be too far off reality
    df = fill_data_gaps_on_integration_bounds(df, data=data, center_window=center_window)
    
    # Interpolate the missing values
    df[data] = df[data].interpolate(method='time', limit_area='inside') # only interpolate between valid values
    
    # Get an x-axis for integration in seconds
    #df['time_diff_in_sec'] = (df.index - df.index[0]).astype('timedelta64[s]')
    df['time_diff_in_sec'] = (df.index - df.index[0]).total_seconds()
    
    # Convert all arrays from shape (n,) to (1,n), so we can use np.hstack() later
    # left side of interval
    x1 = df['time_diff_in_sec'].values.reshape(-1,1)
    y1 = df[data].values.reshape(-1,1)
    
    # right side of interval: shift entire array by 1 timestep, to get the consectuive values for everything 
    # in x1 and x2 at the same indices as x1 and y1 use
    y2 = df[data].shift(periods=-1).values.reshape(-1,1) 
    x2 = df['time_diff_in_sec'].shift(periods=-1).values.reshape(-1,1)
    
    # stick paired values side-by-side, so we have each pair of consecutive values as a single entry
    x = np.hstack([x1,x2])
    y = np.hstack([y1,y2])
    
    # Integrate to find volumes between current timestamps
    df['volume_between_now_and_next'] = integrate.trapezoid(y,x=x)
    
    # Resample to hourly using the whole hour as a mid-point
    if center_window:
        #  Source: https://stackoverflow.com/questions/59948078/resample-to-pandas-dataframe-to-hourly-using-hour-as-mid-point
        #  This shifts the entire time index by 30 minutes, meaning what was previously 11.30 is now 12.00
        #  Because pandas.resample() operates on whole hours, it will with this new index therefore be working with the
        #  data from 11.30-12.30, even though the new index would list times between 12.00-13.00. One extra benefit of
        #  this is that the resampled data will automatically be associated with the correct mid-point time (e.g. 12.00).
        df_H     = df.shift(freq='30min').resample('1H', closed='left').sum()
        df_H_RHS = df.shift(freq='30min').resample('1H', closed='right').sum() # This is solely to check if we have an obs at
    else:                                                                      #  the right edge of the summing interval
        df_H     = df.resample('1H').sum()
        df_H_RHS = df.resample('1H', closed='right').sum()
        
    # Convert volume back into flux
    df_H[data] = df_H['volume_between_now_and_next'] / 3600
    
    # CLEAN-UP
    # Remove the interpolated values that are not based on any observations
    df_H_RHS = df_H_RHS.rename(columns={'is_obs':'is_obs_RHS'}) # cannot merge this and DF_H (next line) if names are the same
    df_H = df_H.join(df_H_RHS['is_obs_RHS'], how='outer') # outer matches on index, keeping _all_ unique index values
    df_H['based_on_obs'] = 0 # New column "is this based on at least one observation?"
    df_H.loc[(df_H['is_obs'] > 0) | (df_H['is_obs_RHS'] > 0), 'based_on_obs'] = 1
    df_H.loc[df_H['based_on_obs'] == 0, data] = np.nan
    
    # Remove the values in the dataframe that are for timesteps before we have obs, or after obs have ended
    df_H = df_H[df_H.index >= df.index[0]]
    df_H = df_H[df_H.index <= df.index[-1]]
    
    # Account for ICE data values
    # < TO DO: needs to come earlier than here >

    # Remove the temporary columns in the hourly dataframe
    if not keep_tmp_col:
        df_H = df_H.drop(['time_diff_in_sec','volume_between_now_and_next','is_obs','is_obs_RHS'], axis=1)
    
    return df, df_H

def return_data_quality_flag_meaning(l,country):
    
    '''Loops over a list of data quality flags and returns a list with the appropriate flag meanings'''
    
    # Sources:
    # https://help.waterdata.usgs.gov/codes-and-parameters/instantaneous-value-qualification-code-uv_rmk_cd
    # https://help.waterdata.usgs.gov/code/grade_cd_query?fmt=html
    # https://wateroffice.ec.gc.ca/contactus/faq_e.html#Q15
    # HYDAT database (daily data quality flags)
    
    # Select the correct dictionary
    if country.lower() == 'usa':
        standards = {# Values below are USGS' quality flags that come with the data. Special conditions (backwater effects
                     #  equipment malfunction, ice affected) are listed as strings in the data column itself. We have removed
                     #  these strings in an earlier step, and traced those occurences in separate variables 
                     #  (e.g. is_ice_affected). There is no need to list those abbreviations in this function, because this
                     #  function is only applied to the data_quality columns and in the USGS case that column does not contain
                     #  these extra flags. 
                     'nan'    : 'Unknown',
                     'A:[0]'  : 'Undefined',
                     'A:<'    : 'Approved, but reported value known to be inaccurate (real value is lower)',
                     'A:>'    : 'Approved, but reported value known to be inaccurate (real value is higher)',
                     'A:[4]'  : 'Approved, but with Incomplete or Partial Aggregated Record',
                     'P:[4]'  : 'Provisional, with Incomplete or Partial Aggregated Record',
                     'P:e'    : 'Provisional AND estimated',
                     'P'      : 'Provisional, not approved',
                     'A:R'    : 'Approved, but revised',
                     'A:e'    : 'Approved AND estimated, with unknown data grade code',
                     'A'      : 'Approved, with unknown data grade code',
                     'A:[93]' : 'Approved, with IV verification DV <= 10 percent diff',
                     'A:[92]' : 'Approved, with IV verification DV <= 5 percent diff',
                     'A:[91]' : 'Approved, with IV verification DV <= 1 percent diff',
                     'A:[90]' : 'Approved, with IV verification DV <= 0.01 orig DV = 0'
                    }
    elif country.lower() == 'can':
        standards = {# Instantaneous data quality flags
                     'nan'                          : 'Unknown',
                     'Provisional/Provisoire:0'     : 'Provisional, flag unknown (not described in WSC docs)',
                     'Provisional/Provisoire:40'    : 'Provisional, dry (water level below sensor)',
                     'Provisional/Provisoire:10'    : 'Provisional, ice-affected',
                     'Provisional/Provisoire:20'    : 'Provisional, estimated',
                     'Provisional/Provisoire:30'    : 'Provisional, partial day (relevant only for daily means)',
                     'Provisional/Provisoire:nan'   : 'Approved, no qualifier specified',
                     'Provisional/Provisoire:-1'    : 'Provisional, no special conditions',
                     'Provisional/Provisoire:50'    : 'Provisional, revised',
                     'Final/Finales:0'              : 'Approved, flag unknown (not described in WSC docs)',
                     'Final/Finales:40'             : 'Approved, dry (water level below sensor)',
                     'Final/Finales:10'             : 'Approved, ice-affected',
                     'Final/Finales:50'             : 'Approved, but revised',
                     'Final/Finales:20'             : 'Approved, but estimated',
                     'Final/Finales:30'             : 'Approved, partial day (relevant only for daily means)',
                     'Final/Finales:nan'            : 'Approved, no qualifier specified',
                     'Final/Finales:-1'             : 'Approved, no special conditions',
                     # Daily data quality flags
                     'A'                            : 'Partial day',
                     'B'                            : 'Ice Conditions',
                     'D'                            : 'Dry',
                     'E'                            : 'Estimated',
                     'S'                            : 'Sample(s) collected this day'
                    }
    
    # Map flags onto meanings
    meanings = [standards[item] for item in l]
    
    return meanings

def select_minimal_usgs_data_quality_flag(flags):
    
    '''Compares a list of flags to defined standards and selects the lowest quality flag in the list as a string'''
    
    # https://help.waterdata.usgs.gov/codes-and-parameters/instantaneous-value-qualification-code-uv_rmk_cd
    # https://help.waterdata.usgs.gov/code/grade_cd_query?fmt=html
    
    # quality order
    standards = ['nan',    # Unknown
                 'A:[0]',  # Undefined
                 'A:<',    # Approved, but reported value known to be inaccurate (real value is lower)
                 'A:>',    # Approved, but reported value known to be inaccurate (real value is higher)
                 'A:[4]',  # Approved, but with Incomplete or Partial Aggregated Record
                 'P:e',    # Provisional AND estimated
                 'P',      # Provisional, not approved
                 'A:R',    # Approved, but revised
                 'A:e',    # Approved AND estimated, with unknown data grade code
                 'A',      # Approved, with unknown data grade code
                 'A:[93]', # Approved, with IV verification DV <= 10 percent diff
                 'A:[92]', # Approved, with IV verification DV <= 5 percent diff
                 'A:[91]', # Approved, with IV verification DV <= 1 percent diff
                 'A:[90]', # Approved, with IV verification DV <= 0.01 orig DV = 0
                ]
    
    if len(flags) > 0:
        order = []
        for flag in flags:
            if (type(flag) == float) and (np.isnan(flag)): flag = 'nan'
            order.append(standards.index(flag))
        min_qc = standards[min(order)]
    else:
        min_qc = 'n/a'
    
    return min_qc

def select_minimal_wsc_data_quality_flag(flags):
    
    '''Compares a list of flags to defined standards and selects the lowest quality flag in the list as a string'''
    
    # https://collaboration.cmc.ec.gc.ca/cmc/hydrometrics/www/Document/WebService_Guidelines.pdf
    # https://wateroffice.ec.gc.ca/contactus/faq_e.html#Q15
    
    # quality order
    standards = ['Provisional/Provisoire:0',     # Provisional, flag unknown (not described in WSC docs)
                 'Provisional/Provisoire:40',    # Provisional, dry (water level below sensor)
                 'Provisional/Provisoire:10',    # Provisional, ice-affected
                 'Provisional/Provisoire:20',    # Provisional, estimated
                 'Provisional/Provisoire:30',    # Provisional, partial day (relevant only for daily means)
                 'Provisional/Provisoire:nan',   # Approved, no qualifier specified
                 'Provisional/Provisoire:-1',    # Provisional, no special conditions
                 'Provisional/Provisoire:50',    # Provisional, revised
                 'Final/Finales:0',              # Approved, flag unknown (not described in WSC docs)
                 'Final/Finales:40',             # Approved, dry (water level below sensor)
                 'Final/Finales:10',             # Approved, ice-affected
                 'Final/Finales:50',             # Approved, but revised
                 'Final/Finales:20',             # Approved, but estimated
                 'Final/Finales:30',             # Approved, partial day (relevant only for daily means)
                 'Final/Finales:nan',            # Approved, no qualifier specified
                 'Final/Finales:-1'              # Approved, no special conditions
                 ]
    
    if len(flags) > 0:
        order = []
        for flag in flags:
            if (type(flag) == float) and (np.isnan(flag)): flag = 'nan'
            order.append(standards.index(flag))
        min_qc = standards[min(order)]
    else:
        min_qc = 'n/a'
    
    return min_qc

def assign_hourly_quality_flag(df, df_H, country, center_window=True):
    
    '''Checks the quality of data our hourly estimates are based on and assigns ice and qc flags accordingly'''
    
    # 0. Prep new columns
    df_H['is_ice_affected'] = 0
    df_H['is_malfunction_affected'] = 0
    df_H['is_backwater_affected'] = 0
    df_H['is_below_sensor_level'] = 0
    df_H['minimum_data_quality'] = 'n/a'
    
    # Check for each value - maybe this can be more elegant but it's a one-off thing so looping is fine
    for ix,row in df_H.iterrows():
        
        # 1. Find the time window this value is based on
        if center_window:
            time_mask = slice(row.name - timedelta(minutes=30), row.name + timedelta(minutes=30))
        else:
            time_mask = slice(row.name, row.name + timedelta(hours=1))
        
        # 2. Select the subset of the original dataframe for this time window
        subset = df.loc[time_mask].copy()
        
        # 3. Assign quality flags to lists (data, ice)
        if country == 'USA':
            
            # USGS stores certain conditions as strings in the flow observation field
            # https://help.waterdata.usgs.gov/codes-and-parameters/instantaneous-and-daily-value-status-codes
            if subset['obs_00060'].dtype == 'O': # Object, we get this if mixing float and str
                if (subset['obs_00060'].str.lower() == 'ice').any():
                    df_H.at[row.name,'is_ice_affected'] = 1
                if (subset['obs_00060'].str.lower() == 'eqp').any():
                    df_H.at[row.name,'is_malfunction_affected'] = 1
                if (subset['obs_00060'].str.lower() == 'bkw').any():
                    df_H.at[row.name,'is_backwater_affected'] = 1
        
            # Use a dedicated function to select the lowest quality flag from the 'data quality' field
            flags = subset['obs_00060_cd'].to_list()
            df_H.at[row.name, 'minimum_data_quality'] = select_minimal_usgs_data_quality_flag(flags)
            
        elif country == 'CAN':
            
            # WSC stores certain conditions as integer values in a dedicated column
            if (subset['Qualifier/Qualificatif'] == 10).any():
                df_H.at[row.name,'is_ice_affected'] = 1
            if (subset['Qualifier/Qualificatif'] == 40).any():
                df_H.at[row.name,'is_below_sensor_level'] = 1
            
            # Use a dedicated function to select the lowest quality flag from the 'data quality' field
            flags = subset['Approval/Approbation'].to_list() # Get provisional/approved tag
            dataq = subset['Qualifier/Qualificatif'].to_list() # Get data quality tag
            merge = ["{}:{}".format(a_, b_) for a_, b_ in zip(flags, dataq)] # combine into a single list 
            df_H.at[row.name, 'minimum_data_quality'] = select_minimal_wsc_data_quality_flag(merge)
            
        elif country == 'MEX':
            
            # < TO DO >
            a = 0
    
    return df_H

# To netcdf
# -----------
def prep_subdaily_country_csv_for_netcdf(csv_path,country):
    
    '''Loads a .csv with observed flow data and processes according to the country the data originates from'''
    
    # Load the data
    csv = pd.read_csv(csv_path, index_col=0, parse_dates=True)    
    
    # General processing
    csv.index.name = 'time'
    data_name = 'q_obs'
    flag_name = 'q_obs_data_quality'
    data_conversion = 0.0283168466 # m^3 ft^-3
    
    # Ensure the index is not timezone-aware, because this trips up the conversion to netcdf dimension later
    csv.index = csv.index.tz_localize(None)
    
    # Country-specific processing
    if country.lower() == 'usa':
        csv = csv.rename(columns={'obs_00060': data_name})
        
        # Unit conversion
        # We know all the USGS data is in cubic feet per second, because we checked this in 1b_usa_flow_obs_to_utc.ipynb
        print('Warning: converting data from units feet^3 s^-1 to m^3 s^-1')
        csv[data_name] = csv[data_name] * data_conversion
        
        # Drop the columns we added to the csv's for the Canadian data
        if 'is_below_sensor_level' in csv.columns:
            csv = csv.drop(columns=['is_below_sensor_level'])
        
    elif country.lower() == 'can':
        csv = csv.rename(columns={'Value/Valeur': data_name})
        
        # Unit conversion
        # Data comes in m3/s (p.3): https://collaboration.cmc.ec.gc.ca/cmc/hydrometrics/www/Document/WebService_Guidelines.pdf
        
        # Drop the columns we added to the csv's for the US data
        if 'is_malfunction_affected' in csv.columns:
            csv = csv.drop(columns=['is_malfunction_affected'])
        if 'is_backwater_affected' in csv.columns:
            csv = csv.drop(columns=['is_backwater_affected'])
        
    # Check that all data values are derived from observations, and proceed if so
    assert all(csv['based_on_obs'][csv['q_obs'].notna()] == 1), f'Not all data values in {csv_path} are derived from observations. Aborting.'
    csv = csv.drop(columns=['based_on_obs'])
    
    # Find and rename all auxilliary variables
    csv = csv.rename(columns={'minimum_data_quality': flag_name})
    for column in csv.columns:
        if not 'q_obs' in column:
            csv = csv.rename(columns={column: 'q_obs_'+column}) # prepends 'q_obs' to any remaining variables
    
    return csv

def prep_daily_country_csv_for_netcdf(csv_path,country,tz):
    
    # 0. General processing settings
    indx_name = 'time'
    data_name = 'q_obs'
    flag_name = 'q_obs_data_quality'
    data_conversion = 0.0283168466 # m^3 ft^-3
    
    # USA-specific codes
    usa_streamflow_codes = ['***', 'Bkw', 'Dis', 'Eqp', 'Ice', 'Rat', 'Ssn'] # https://help.waterdata.usgs.gov/codes-and-parameters/instantaneous-and-daily-value-status-codes
    usa_required_columns = ['q_obs', 'q_obs_data_quality','tz_cd',
                            'is_ice_affected','is_malfunction_affected','is_backwater_affected']
    
    # CAN-specific codes
    can_streamflow_codes = ['A','B','D','E','S'] # HYDAT; see code blocks below
    can_required_columns = ['q_obs', 'q_obs_data_quality','tz_cd',
                            'is_ice_affected','is_partial_day','is_dry_day','is_estimated_value']
    
    # 1. Load the data
    csv = pd.read_csv(csv_path, index_col=0, parse_dates=True)
    
    # 2. Ensure the index is not timezone-aware, because this trips up the conversion to netcdf dimension later
    csv.index = csv.index.tz_localize(None)
    csv.index.name = indx_name
    
    # 3. Prep new columns for status code processing
    csv['is_ice_affected'] = 0
    csv['is_malfunction_affected'] = 0
    csv['is_backwater_affected'] = 0
    csv['is_partial_day'] = 0
    csv['is_dry_day'] = 0
    csv['is_estimated_value'] = 0

    # 4. Country-specific processing
    if country.lower() == 'usa':
        
        # 3a. Replace column names
        csv = csv.rename(columns={'obs_00060_00003': data_name,
                                  'obs_00060_00003_cd': flag_name})
        
        # 3b. Move condition strings from the observation column into dedicated columns
        for code in usa_streamflow_codes:
            
            # Handle specific flags
            if code == 'Ice': 
                csv.loc[csv[data_name] == code, 'is_ice_affected'] = 1
            elif code == 'Eqp': 
                csv.loc[csv[data_name] == code, 'is_malfunction_affected'] = 1
            elif code == 'Bkw': 
                csv.loc[csv[data_name] == code, 'is_backwater_affected'] = 1
            
            # Remove all flags so we can transform to float
            csv[data_name] = csv[data_name].replace(code,np.nan)
        
        # 3c. Transform to float and remove negative values
        csv[data_name] = csv[data_name].astype('float')
        csv.loc[csv[data_name] < 0, data_name] = np.nan
        
        # 3d. Unit conversion
        # We know all the USGS data is in cubic feet per second, because we checked this in 1d_usa_daily_flow_obs_to_csv.ipynb
        print('Warning: converting data from units feet^3 s^-1 to m^3 s^-1')
        csv[data_name] = csv[data_name] * data_conversion
        
        # 3e. Drop any columns we don't need later
        for column in csv.columns:
            if column not in usa_required_columns: 
                csv = csv.drop(columns=[column])
        
        # 3f. Extra check: ensure that metadata timezone matches csv timezone
        assert csv['tz_cd'].unique()[0] == tz, "metadata and csv timezone specification don't match"
    
    elif country.lower() == 'can':
        
        # 3a. Replace column names
        csv = csv.rename(columns={'FLOW': data_name,
                                  'SYMBOL': flag_name})
        
        # 3b. Move condition strings from the observation column into dedicated columns
        for code in can_streamflow_codes:
            
            # Handle specific flags
            if code == 'A': 
                csv.loc[csv[data_name] == code, 'is_partial_day'] = 1
            elif code == 'B': 
                csv.loc[csv[data_name] == code, 'is_ice_affected'] = 1
            elif code == 'D': 
                csv.loc[csv[data_name] == code, 'is_dry_day'] = 1
            elif code == 'E': 
                csv.loc[csv[data_name] == code, 'is_estimated_value'] = 1
        
        # 3c. Remove negative values
        csv.loc[csv[data_name] < 0, data_name] = np.nan
        
        # 3d. Add the timezone locale
        # Get station location and match this to the "Standard timezone" shapefile
        csv['tz_cd'] = tz
        
        # 3e. Drop any columns we don't need later
        for column in csv.columns:
            if column not in can_required_columns: 
                csv = csv.drop(columns=[column])
    
    # 5. Find and rename all auxilliary variables
    csv = csv.rename(columns={'minimum_data_quality': flag_name})
    for column in csv.columns:
        if not 'q_obs' in column:
            csv = csv.rename(columns={column: 'q_obs_'+column}) # prepends 'q_obs' to any remaining variables
    
    return csv

   
def subdaily_flow_csv_to_netcdf(csv, nc_path, country, station, tz, center_window):
    
    '''Converts a standardized csv file with flow observations to xarray data set and saves as netcdf'''
    
    # 1. Define standard values
    # -------------------------
    
    # Auxiliary
    global_att_countries = ['USA', 'CAN', 'MEX']
    global_att_i = global_att_countries.index(country)
    global_att_now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Global attributes
    global_att_ttl = 'CAMELS-spat streamflow data'
    global_att_con = 'CF-1.10'
    global_att_src = 'Streamflow derived from observed water levels'
    global_att_ins = ['United States Geological Survey',
                      'Water Survey of Canada']
    global_att_ref = [('U.S. Geological Survey, 2016, National Water Information System data available ' +
                       'on the World Wide Web (USGS Water Data for the Nation), accessed 2023-03-23, at '+
                       'URL [http://waterdata.usgs.gov/nwis/]'),
                      ('Original data extracted from the Environment and Climate Change Canada Real-time ' +
                       'Hydrometric Data web site (https://wateroffice.ec.gc.ca/mainmenu/real_time_data_index_e.html) ' + 
                       'on 2023-04-05')]
    global_att_his = (f'{global_att_now} | File prepared using CAMELS-spat scripts. See:' + 
                       'https://github.com/CH-Earth/camels-spat')
    global_att_com = 'n/a'
    
    # Data variables
    q_obs_unit = 'm3 s-1'
    q_obs_long = 'observed streamflow values'
    q_obs_anc = [column for column in csv.columns if '_is_' in column] # Get names of all ancillary variables in .csv 
    q_obs_anc.append('q_obs_data_quality') # add the 'q_obs_data_quality' variable that's not captured by the above
    q_obs_anc = ' '.join([f"'{anc}'" for anc in q_obs_anc]) # convert full list into single string
    
    # Time settings
    time_unit = 'minutes since 1950-01-01 00:00:00'
    time_cal = 'proleptic_gregorian'
    
    # 2. Create a basic data set to build from
    ds = csv.to_xarray()
    
    # 3. Global attributes
    ds.attrs['title'] = global_att_ttl
    ds.attrs['conventions'] = global_att_con
    ds.attrs['source'] = global_att_src
    ds.attrs['country'] = country
    ds.attrs['station'] = station
    ds.attrs['institution'] = global_att_ins[global_att_i]
    ds.attrs['references'] = global_att_ref[global_att_i]
    ds.attrs['history'] = global_att_his
    #ds.attrs['comment'] = global_att_com

    # 4a. Time attributes (coordinate already exists)
    # NOTE: attributes 'units' and 'calendar' are automatically specified when writing to netcdf
    #       This can be checked by saving to netcdf, and then loading as follows: xr.open_dataset(nc_path, decode_times=False)
    ds.time.attrs['standard_name'] = 'time'
    ds.time.attrs['bounds'] = 'time_bnds'
    ds.time.encoding['units'] = time_unit
    ds.time.encoding['calendar'] = time_cal
        
    # 4b. Time bounds variable
    ds = ds.assign_coords(nbnds=[1,2])
    if center_window:    
        ds = ds.assign(time_bnds=(['nbnds','time'], [csv.index - pd.Timedelta('30min'), csv.index + pd.Timedelta('30min')]))
    else:
        ds = ds.assign(time_bnds=(['nbnds','time'], [csv.index, csv.index + pd.Timedelta('1h')]))
    ds.nbnds.attrs['standard_name'] = 'bounds for timestep intervals'
    ds.time_bnds.attrs['long_name'] = 'start and end points of each time step'
    ds.time_bnds.attrs['time_zone'] = tz
    
    # 5. Observed streamflow
    ds.q_obs.attrs['units'] = q_obs_unit
    ds.q_obs.attrs['long_name'] = q_obs_long
    ds.q_obs.attrs['cell_methods'] = 'time:mean' # indicating that values are average values over the timestep
    ds.q_obs.attrs['ancillary_variables'] = q_obs_anc
    ## TO DO: add other variables to ancillary_variables list
    
    # 6. Data quality flags
    flags = [str(s) for s in csv['q_obs_data_quality'].unique()]
    flags.sort()
    meanings = return_data_quality_flag_meaning(flags,country)
    ds.q_obs_data_quality.attrs['standard_name'] = 'quality_flag'
    ds.q_obs_data_quality.attrs['long_name'] = 'lowest data quality flag listed in the values used to generate an average flow value for each timestep'
    ds.q_obs_data_quality.attrs['flag_values'] = ' '.join([f"'{flag}'" for flag in flags])
    ds.q_obs_data_quality.attrs['flag_meanings'] = ' '.join([f"'{meaning}'" for meaning in meanings])
    
    # 7. Other status variables
    for variable in ds.variables:
        if '_is_' in variable:
            ds[variable].attrs['standard_name'] = 'quality_flag'
            ds[variable].attrs['long_name'] = 'flag indicating if main variable is affected by process in variable name'
            ds[variable].attrs['flag_values'] = "'0' '1'"
            ds[variable].attrs['flag_meanings'] = "'no' 'yes'"
    
    # Save to file
    ds = ds.drop_indexes(['time','nbnds'])
    ds.to_netcdf(nc_path)
    
    return ds

def daily_flow_csv_to_netcdf(csv, nc_path, country, station):
    
    '''Converts a standardized csv file with flow observations to xarray data set and saves as netcdf'''
    
    # 1. Define standard values
    # -------------------------
    
    # Auxiliary
    global_att_countries = ['USA', 'CAN', 'MEX']
    global_att_i = global_att_countries.index(country)
    global_att_now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Global attributes
    global_att_ttl = 'CAMELS-spat streamflow data'
    global_att_con = 'CF-1.10'
    global_att_src = 'Streamflow derived from observed water levels'
    global_att_ins = ['United States Geological Survey',
                      'Water Survey of Canada']
    global_att_ref = [('U.S. Geological Survey, 2016, National Water Information System data available ' +
                       'on the World Wide Web (USGS Water Data for the Nation), accessed 2023-03-23, at '+
                       'URL [http://waterdata.usgs.gov/nwis/]'),
                      ('Original data extracted from the Environment and Climate Change Canada Real-time ' +
                       'Hydrometric Data web site (https://wateroffice.ec.gc.ca/mainmenu/real_time_data_index_e.html) ' + 
                       'on 2023-04-05')]
    global_att_his = (f'{global_att_now} | File prepared using CAMELS-spat scripts. See:' + 
                       'https://github.com/CH-Earth/camels-spat')
    global_att_com = (f'{global_att_ins[global_att_i]} calculates daily average flow values' +
                      ' in the station\'s standard time (i.e., not UTC). See: variable time_bnds.')
    
    # Data variables
    q_obs_unit = 'm3 s-1'
    q_obs_long = 'observed streamflow values'
    q_obs_anc = [column for column in csv.columns if '_is_' in column] # Get names of all ancillary variables in .csv 
    q_obs_anc.append('q_obs_data_quality') # add the 'q_obs_data_quality' variable that's not captured by the above
    q_obs_anc = ' '.join([f"'{anc}'" for anc in q_obs_anc]) # convert full list into single string
    
    # Time settings
    time_unit = 'minutes since 1950-01-01 00:00:00'
    time_cal = 'proleptic_gregorian'
    
    # Time settings - ensure the data only specifies a single time zone that was used for calculating averages
    #                 Communications with USGS and WSC state that this should be the case
    assert len(csv['q_obs_tz_cd'].unique()) == 1, "Multiple timezones specified in csv; there should be only one"
    time_original_tz = csv['q_obs_tz_cd'].unique()[0]
    
    # 2. Create a basic data set to build from
    ds = csv.to_xarray()
        
    # 3. Global attributes
    ds.attrs['title'] = global_att_ttl
    ds.attrs['conventions'] = global_att_con
    ds.attrs['source'] = global_att_src
    ds.attrs['country'] = country
    ds.attrs['station'] = station
    ds.attrs['institution'] = global_att_ins[global_att_i]
    ds.attrs['references'] = global_att_ref[global_att_i]
    ds.attrs['history'] = global_att_his
    ds.attrs['comment'] = global_att_com
    
    # 4a. Time attributes (coordinate already exists)
    # NOTE: attributes 'units' and 'calendar' are automatically specified when writing to netcdf
    #       This can be checked by saving to netcdf, and then loading as follows: xr.open_dataset(nc_path, decode_times=False)
    ds.time.attrs['standard_name'] = 'time'
    ds.time.attrs['bounds'] = 'time_bnds'
    ds.time.encoding['units'] = time_unit
    ds.time.encoding['calendar'] = time_cal
        
    # 4b. Time bounds variable
    ds = ds.assign_coords(nbnds=[1,2])
    ds = ds.assign(time_bnds=(['nbnds','time'],
                              [csv['time_bnds_l'], csv['time_bnds_r']]))
    ds.nbnds.attrs['standard_name'] = 'bounds for timestep intervals'
    ds.time_bnds.attrs['long_name'] = 'start and end points of each time interval'
    #ds.time_bnds.attrs['time_zone'] = 'UTC'
    #ds.time_bnds.attrs['station_standard_time'] = time_original_tz
    ds.time_bnds.attrs['time_zone'] = time_original_tz

    # 5. Observed streamflow
    ds.q_obs.attrs['units'] = q_obs_unit
    ds.q_obs.attrs['long_name'] = q_obs_long
    ds.q_obs.attrs['cell_methods'] = 'time:mean' # indicating that values are average values over the timestep
    ds.q_obs.attrs['ancillary_variables'] = q_obs_anc
    ## TO DO: add other variables to ancillary_variables list
    
    # 6. Data quality flags
    flags = [str(s) for s in csv['q_obs_data_quality'].unique()]
    flags.sort()
    while ' ' in flags: flags.remove(' ')  # Sometimes we have empty spaces with no specific meaning in the data quality column: take those out
    meanings = return_data_quality_flag_meaning(flags,country)
    ds.q_obs_data_quality.attrs['standard_name'] = 'quality_flag'
    ds.q_obs_data_quality.attrs['long_name'] = 'data quality flag'
    ds.q_obs_data_quality.attrs['flag_values'] = ' '.join([f"'{flag}'" for flag in flags])
    ds.q_obs_data_quality.attrs['flag_meanings'] = ' '.join([f"'{meaning}'" for meaning in meanings])
    
    # 7. Other status variables
    for variable in ds.variables:
        if '_is_' in variable:
            ds[variable].attrs['standard_name'] = 'quality_flag'
            ds[variable].attrs['long_name'] = 'flag indicating if main variable is affected by process in variable name'
            ds[variable].attrs['flag_values'] = "'0' '1'"
            ds[variable].attrs['flag_meanings'] = "'no' 'yes'"
    
    # 8. Remove the timezone variables we added to get the time_bnds
    ds = ds.drop_vars(['q_obs_tz_cd', 'time_bnds_l', 'time_bnds_r'])
    
    # Save to file
    ds = ds.drop_indexes(['time','nbnds'])
    ds.to_netcdf(nc_path)
        
    return ds

# Update netcdf
# -------------

def add_time_bnds(file, dataset=[], timezone=[]):

    '''Adds a time_bnds variable to a netcdf'''

    with nc4.Dataset(file, 'a') as f: # (a)ppend

        # Check that the approach within this function is valid
        assert 'hours' in f['time'].getncattr('units'), f'ERROR: time units in {file} not in hours'

        # Check that a source dataset is specified
        assert (dataset.lower() == 'era5') | (dataset.lower() == 'em-earth'), \
            f'add_time_bnds() contains no settings for dataset = {dataset}'

        # Connect variable 'time_bounds' to variable 'time' through time attribute 'bounds'
        f['time'].setncattr('bounds','time_bnds')

        # Add nbnds dimension
        f.createDimension('nbnds', 2)
        f.createVariable('nbnds', 'i', 'nbnds')
        f.variables['nbnds'][:] = [1,2]
        f['nbnds'].setncattr('standard_name','bounds for timestep intervals')

        # Add time_bnds variable
        f.createVariable('time_bnds', f.variables['time'].datatype, ('nbnds','time'), fill_value = False, zlib=True, shuffle=True)
        f['time_bnds'].setncattr('long_name', 'start and end points of each time step')
        f['time_bnds'].setncattr('time_zone', timezone)
        f['time_bnds'].setncattr('calendar', f['time'].getncattr('calendar'))
        f['time_bnds'].setncattr('units', f['time'].getncattr('units'))

        # Add the actual data
        if dataset.lower() == 'era5':
            f['time_bnds'][:] = np.array([f['time'][:]-1, f['time'][:]]) # Period-ending timestamps: t(n) is valid over t(n-1) to t(n)
        if dataset.lower() == 'em-earth':
            f['time_bnds'][:] = np.array([f['time'][:], f['time'][:]+1]) # Period-starting timestamps: t(n) is valid over t(n) to t(n+1)

        # Update the file history
        new_history = f' On {time.ctime(time.time())}: add_time_bnds().'
        old_history = f.History
        hist = f'{old_history} {new_history}'
        f.setncattr('History',hist)