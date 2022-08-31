'''Contains functions to interact with the configuration file.'''

def my_tester():
    print('Successfully executed test function.')

def enforce_type( value,type ):
    
    '''Converts [value] to [type], as specified in the config file.'''
    if type == 'string':
        value = str(value)
    elif type == 'float':
        value = float(value)
    elif type == 'integer':
        value = int(value)
    else:
        print('WARNING in enforce_type(): type not specified, defaulting to returning value as string.')
    return value
    
def read_from_config( file, setting ):
    
    '''Reads [setting] from configuration [file].'''
    # Open 'config.txt' and ...
    with open(file) as contents:
        for line in contents:
            
            # ... find the line with the requested setting
            if setting in line:
                break
    
    # Line comes out as: [setting] | [value] | [type] # [comment]
    # Harmonize the delimiters
    line = line.replace('#','|')
    
    # Extract the setting's value and type
    value = line.split('|')[1].strip()
    type = line.split('|')[2].strip()
    
    # Type conversion
    value = enforce_type(value,type)
    
    return value