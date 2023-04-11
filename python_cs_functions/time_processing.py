'''Contains functions for temporal processing of downloaded data (timezones, averaging, etc.)'''

# Time zone conversion
# ---------------------------------------
def tz_abbreviation_to_utc(tz_str):
    '''Takes a timezone abbreviation used by USGS and converts to UTC-xx'''
    
    # Check inputs
    assert type(tz_str) == str, f'Input tz_str must be {str} but is {type(tz_str)}'
    
    # Define abbreviation-to-UTC dictionary
    dict_abbr_to_utc = {'CDT': 'UTC-05', # Ambiguous: Central Daylight Time.  Possibly legacy part of pytz: CST6CDT
                        'CST': 'UTC-06', # Ambiguous: Central Standard Time.  Not in pytz
                        'EDT': 'UTC-04', #            Eastern Daylight Time.  legacy part of pytz:          EST5EDT
                        'EST': 'UTC-05', #            Eastern Standard Time.  Part of pytz:                 EST
                        'MDT': 'UTC-06', #            Mountain Daylight Time. legacy part of pytz:          MST7MDT
                        'MST': 'UTC-07', # Ambiguous: Mountain Standard Time. Possibly art of pytz:         MST
                        'PDT': 'UTC-07', #            Pacific Daylight Time.  Legacy part of pytz:          PST8PDT
                        'PST': 'UTC-08'} #            Pacific Standard Time.  Not in pytz
    
    return dict_abbr_to_utc[tz_str]

def relative_utc_to_timedelta(utc_str):
    '''Takes a UTC-xx string and converts to Timedelta string'''
    
    # Check inputs
    assert type(utc_str) == str, f'Input tz_str must be {str} but is {type(utc_str)}'
    
    # Define UTC-to-timedelta dictionary
    dict_utc_to_timedelta = {'UTC-04': '+4:00:00',
                             'UTC-05': '+5:00:00',
                             'UTC-06': '+6:00:00',
                             'UTC-07': '+7:00:00',
                             'UTC-08': '+8:00:00'}
    
    return dict_utc_to_timedelta[utc_str]

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
import numpy as np
import pandas as pd
from scipy import integrate
from datetime import timedelta

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
    df['time_diff_in_sec'] = (df.index - df.index[0]).astype('timedelta64[s]')
    
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

def select_minimal_usgs_data_quality_flag(flags):
    
    '''Compares a list of flags to defined standards and selects the lowest quality flag in the list as a string'''
    
    # https://help.waterdata.usgs.gov/codes-and-parameters/instantaneous-value-qualification-code-uv_rmk_cd
       
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