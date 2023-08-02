'''Contains functions to interact with the HYDAT sqlite3 database.'''

import sqlite3
import numpy as np
import pandas as pd

# Source: https://www.sqlitetutorial.net/sqlite-python/sqlite-python-select/
def connect_to_sqlite_database(db_file):
    
    ''' create a database connection to the SQLite database
        specified by the db_file
    :param db_file: database file
    :return: Connection object or None
    '''
    
    conn = None
    try:
        conn = sqlite3.connect(db_file)
    except Error as e:
        print(e)
    
    return conn


def find_all_table_names(conn, to_screen=True):
    
    ''' Reads all available table names in the database.
        Prints to screen by default.'''
    
    # Source: https://stackoverflow.com/a/10746045
    cur = conn.cursor()
    cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = cur.fetchall()
    cur.close()
    
    if to_screen:
        for table in tables:
            print(table)
        
    return tables


def find_table_contents(conn, table_name, to_screen=True):
    
    ''' Reads all rows in a given database table. 
        Prints to screen by default.'''
    
    # Adapted from: https://www.sqlitetutorial.net/sqlite-python/sqlite-python-select/
    cur = conn.cursor()
    cur.execute("SELECT * FROM {};".format(table_name))
    rows = cur.fetchall()
    cur.close()
    
    if to_screen:
        for row in rows:
            print(row)
    
    return rows


def find_table_headers(conn, table_name, to_screen=True):
    
    ''' Reads headers from a row in a given database table. 
        Prints to screen by default.'''
    
    # Source: https://stackoverflow.com/a/7831685
    cur = conn.execute("SELECT * FROM {};".format(table_name))
    headers = [description[0] for description in cur.description]
    cur.close()
    
    if to_screen:
        print(headers)
    
    return headers


def sql_query_to_dataframe(conn, query):
    
    ''' Returns the result of a databse query as a Pandas dataframe'''
    
    # Source: https://stackoverflow.com/a/36029761
    df = pd.read_sql_query(query, conn)
    
    return df


def construct_query_from_dataframe_column(table,field,values):
    
    ''' Constructs an SQL query to obtain an entire table 
        for a given list of values in the table''s field.
        Query format returned:
        "SELECT * FROM [table] WHERE [field] in (list[0],..,list[n])"
        '''
    
    values_str = "('" + "','".join(values) + "')"
    query = "SELECT * FROM {} WHERE {} in {}".format(table,field,values_str)
    
    return query

def location_abbreviation_to_full_name(abbr):
    
    '''Converts two-letter state/territory/province abbreviations in HYDAT to full name'''
    
    # Source 1: https://www.canada.ca/en/revenue-agency/services/tax/businesses/topics/completing-slips-summaries/financial-slips-summaries/return-investment-income-t5/provincial-territorial-codes.html
    # Source 2: https://en.wikipedia.org/wiki/List_of_U.S._state_and_territory_abbreviations
    
    short_to_full = {# Canadian lcoations
            'AB': 'Alberta',
            'BC': 'British Columbia',
            'MB': 'Manitoba',
            'NB': 'New Brunswick',
            'NL': 'Newfoundland and Labrador',
            'NS': 'Nova Scotia',
            'NT': 'Northwest Territories',
            'NU': 'Nunavut',
            'ON': 'Ontario',
            'PE': 'Prince Edward Island',
            'QC': 'Quebec',
            'SK': 'Saskatchewan',
            'YT': 'Yukon',
            # Cross-border basins: USA abbreviations
            'AK': 'Alaska',
            'ID': 'Idaho',
            'ME': 'Maine',
            'MT': 'Montana',
            'WA': 'Washington'}
    
    return short_to_full[abbr]
    
def restructure_daily_flows_in_hydat_table(table, time_s, time_e, site):
    
    '''Takes a dataframe containing a HYDAT table for a given station and extracts daily flow and QC values'''
    
    # Define years and months to cycle over
    years = table['YEAR'].unique()
    months= np.arange(1,13)
    
    # Make a list of monthly dataframes
    dfs = []
    for year in years:
        for month in months:
            
            # Check if this year and month fall within the period of interest
            if ( (year < time_s.year) or 
                 ((year == time_s.year) and (month < time_s.month)) or
                 ((year == time_e.year) and (month > time_e.month)) or
                 (year > time_e.year) ):
                continue
            
            # Select the year and month and continue if this month is in fact part of the data
            mask = (table['YEAR'] == year) & (table['MONTH'] == month) 
            if mask.any():
                
                # Find number of days in the month
                nday = table[mask]['NO_DAYS'].values[0]
                days = np.arange(1,nday+1)
                
                # Make a temporary dataframe for this month
                tmp = pd.DataFrame(data={'site': site, 
                                         'YEAR': year, 'MONTH': month, 'DAY': days, 
                                         'FLOW': np.nan, 'SYMBOL': np.nan})
                
                # Move the daily values into the proper columns
                for day in days:
                    col1 = f'FLOW{day}'
                    col2 = f'FLOW_SYMBOL{day}'
                    flow = table[mask][col1].values[0]
                    symb = table[mask][col2].values[0]
                    tmp.loc[day-1,['FLOW']] = flow
                    tmp.loc[day-1,['SYMBOL']] = symb
                
                # Add this month's dataframe to the list
                dfs.append(tmp)

    # Make a single dataframe with a timeseries, if possible
    if not dfs: # if monthly dataframe list is empty (empty list = False)
        df = pd.DataFrame()
    else:
        df = pd.concat(dfs)
    
        # Convert YEAR, MONTH, DAY columns into a single Date column
        df['Date'] = pd.to_datetime(dict(year=df['YEAR'], month=df['MONTH'], day=df['DAY']))
        df = df.set_index('Date')
        df = df.drop(columns=['YEAR','MONTH','DAY'])
    
    return df