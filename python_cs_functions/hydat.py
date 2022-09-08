'''Contains functions to interact with the HYDAT sqlite3 database.'''

# Source: https://www.sqlitetutorial.net/sqlite-python/sqlite-python-select/
def connect_to_sqlite_database(db_file):
    
    ''' create a database connection to the SQLite database
        specified by the db_file
    :param db_file: database file
    :return: Connection object or None
    '''
    import sqlite3
    
    conn = None
    try:
        conn = sqlite3.connect(db_file)
    except Error as e:
        print(e)
    
    return conn


def find_all_table_names(conn, to_screen=True):
    
    ''' Reads all available table names in the database.
        Prints to screen by default.'''
    import sqlite3
    
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
    import sqlite3
    
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
    import sqlite3
    
    # Source: https://stackoverflow.com/a/7831685
    cur = conn.execute("SELECT * FROM {};".format(table_name))
    headers = [description[0] for description in cur.description]
    cur.close()
    
    if to_screen:
        print(headers)
    
    return headers


def sql_query_to_dataframe(conn, query):
    
    ''' Returns the result of a databse query as a Pandas dataframe'''
    import sqlite3
    import pandas as pd
    
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