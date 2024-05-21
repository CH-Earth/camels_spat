# Function sot process RDRS and Daymet

import netCDF4 as nc

def remove_vars_from_rdrs_download(input_file,output_file,variables_to_keep,
                                   compression=True, complevel=4):
    # Open the NetCDF files
    with nc.Dataset(input_file, 'r') as src, nc.Dataset(output_file, 'w') as dst:
    
        # Copy global attributes
        dst.setncatts({attr: src.getncattr(attr) for attr in src.ncattrs()})
        
        # Copy dimensions
        for name, dimension in src.dimensions.items():
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))
        
        # Copy variables except the ones to be removed
        for name, variable in src.variables.items():
            if name in variables_to_keep:
                # Determine the compression settings
                zlib = compression
                complevel = complevel if compression else None
                
                # Create variable in the destination file with compression
                new_var = dst.createVariable(
                    name, 
                    variable.datatype, 
                    variable.dimensions,
                    zlib=zlib,
                    complevel=complevel
                )
                
                # Copy variable attributes
                new_var.setncatts({attr: variable.getncattr(attr) for attr in variable.ncattrs()})
                # Copy variable data
                new_var[:] = variable[:]
    return