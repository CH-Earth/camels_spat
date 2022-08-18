'''Contains functions to interact with the configuration file.'''

def my_tester():
    print('Successfully executed test function.')
    
def read_from_config( file, setting ):
    
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
           
    # Return this value    
    return value,type