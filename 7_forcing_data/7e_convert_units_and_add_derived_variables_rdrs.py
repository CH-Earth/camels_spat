# Convert units and add derived variables
import glob
import netCDF4 as nc4
import pandas as pd
import shutil
import sys
import xarray as xr
from pathlib import Path
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

# --- Config Handling
# Specify where the config file can be found
config_file = '../0_config/config.txt'

# Get the required info from the config file
data_path = cs.read_from_config(config_file,'data_path')

# CAMELS-spat metadata
cs_meta_path = cs.read_from_config(config_file,'cs_basin_path')
cs_meta_name = cs.read_from_config(config_file,'cs_meta_name')
cs_unusable_name = cs.read_from_config(config_file,'cs_unusable_name')

# Basin folder
cs_basin_folder = cs.read_from_config(config_file, 'cs_basin_path')
basins_path = Path(data_path) / cs_basin_folder

# --- Data loading
# CAMELS-spat metadata file
cs_meta_path = Path(data_path) / cs_meta_path
cs_meta = pd.read_csv(cs_meta_path / cs_meta_name)

# Open list of unusable stations; Enforce reading IDs as string to keep leading 0's
cs_unusable = pd.read_csv(cs_meta_path / cs_unusable_name, dtype={'Station_id': object})

# --- Processing
debug_message = f'\n!!! CHECK DEBUGGING STATUS: \n- Development ongoing\n'

# Get the array index from the command-line argument
if len(sys.argv) != 2:
    print("Usage: python 1e_add_derived_variables.py <array_index>")
    sys.exit(1)
else:
    ix = int(sys.argv[1])

# Define cs_meta row
row = cs_meta.iloc[ix] # needs to be between 0  and 1697
basin_id = row.Country + '_' + row.Station_id

# Check if we need to run downloads for this station at all
missing = cs.flow_obs_unavailable(cs_unusable, row.Country, row.Station_id)
if 'iv' in missing and 'dv' in missing: 
    sys.exit(0) # gracefully exit run, because we have no observations at all for this station

# --- Processing
print('--- Now running basin {}. {}'.format(ix, basin_id))

# Get forcing destination path
raw_fold, _, _ = cs.prepare_forcing_outputs(cs_meta, ix, Path(data_path)/cs_basin_folder) # Returns folders only, not file names

# Find the files
rdrs_files = glob.glob(str(raw_fold / 'rdrs_month' / 'RDRS_*.nc'))
rdrs_files.sort()

# Loop over the files and create a backup
backup_folder = raw_fold / 'rdrs_backup_before_units_pet'
backup_folder.mkdir(parents=True, exist_ok=True)
for rdrs_file in rdrs_files:
    filename = Path(rdrs_file).name
    rdrs_file_backup = filename.replace('RDRS_', 'RDRS_beforeUnitsPET_')
    shutil.copyfile(rdrs_file, backup_folder / rdrs_file_backup)

# Loop over the files and update the units
for rdrs_file in rdrs_files:
    filename = Path(rdrs_file).name
    backup_file = filename.replace('RDRS_', 'RDRS_beforeUnitsPET_')
    ds = xr.open_dataset(backup_folder/backup_file)
    ds = cs.convert_rdrs_variables(ds, cs.create_rdrs_unit_conversion_dict())
    ds.to_netcdf(rdrs_file)

# Loop over the files and add the vapor pressure and pet variables
for rdrs_file in rdrs_files:
    with nc4.Dataset(rdrs_file, 'r+') as f:
        f = cs.derive_vapor_pressure(f, 
                                        air_pressure='RDRS_v2.1_P_P0_SFC', # in [Pa] after conversion above
                                        specific_humidity='RDRS_v2.1_P_HU_1.5m',
                                        Pa_to_kPa=True, # need to convert air pressure units
                                        new_name='e', # same as in ERA5 files
                                        dims='rlat/rlon')
        f = cs.derive_penman_monteith_pet(f,
                                            air_pressure='RDRS_v2.1_P_P0_SFC',
                                            air_temperature='RDRS_v2.1_P_TT_1.5m',
                                            relative_humidity='RDRS_v2.1_P_HR_1.5m',
                                            shortwave_radiation='RDRS_v2.1_P_FB_SFC',
                                            longwave_radiation='RDRS_v2.1_P_FI_SFC',
                                            wind_speed='RDRS_v2.1_P_UVC_10m',
                                            wind_height=10,
                                            ground_heat=0,
                                            new_name='PET',
                                            to_kg_m2_s=True)

# --- Functions
'''
def create_rdrs_unit_conversion_dict():
    #
    # List the variables, their units, and what we want them to be
    # Precipitation quantity:           RDRS_v2.1_A_PR0_SFC,    [m],        [kg m-2 s-1]
    # Downward solar flux:              RDRS_v2.1_P_FB_SFC,     [W m-2],    [W m-2]
    # Surface incoming infrared flux:   RDRS_v2.1_P_FI_SFC,     [W m-2],    [W m-2]
    # Geopotential height:              RDRS_v2.1_P_GZ_SFC,     [dam],      < leave as is >
    # Relative humidity:                RDRS_v2.1_P_HR_1.5m,    [1],        [kPa kPa-1]
    # Specific humidity:                RDRS_v2.1_P_HU_1.5m,    [kg kg-1],  [kg kg-1]
    # Surface pressure:                 RDRS_v2.1_P_P0_SFC,     [mb],       [Pa]
    # Air temperature:                  RDRS_v2.1_P_TT_1.5m,    [C],        [K]
    # Corrected U-component wind:       RDRS_v2.1_P_UUC_10m,    [kts],      [m s-1]
    # Corrected V-component wind:       RDRS_v2.1_P_VVC_10m,    [kts],      [m s-1]
    # Wind modulus from UU and VV:      RDRS_v2.1_P_UVC_10m,    [kts],      [m s-1]
    #
    # Constants
    is_same           = 1       # No conversion needed
    rho_water         = 1000    # [kg m-3]
    seconds_per_hour  = 3600    # [s]
    celsius_to_kelvin = 273.15  # [degrees]
    milibar_to_pascal = 100     # [Pa]
    kts_to_ms         = 0.514   # [kts] * 0.514 = [m s-1]
    #
    # Find conversion factors SPECIFIC to RDRSv2.1
    new_dict = {'RDRS_v2.1_A_PR0_SFC':  rho_water / seconds_per_hour, # [m] * [kg m-3] / [s] = [kg m-2 s-1]
                'RDRS_v2.1_P_FB_SFC':   is_same, # [W m-2]
                'RDRS_v2.1_P_FI_SFC':   is_same, # [W m-2]
                'RDRS_v2.1_P_GZ_SFC':   is_same, # [dam]
                'RDRS_v2.1_P_HR_1.5m':  is_same, # No value conversion needed, but units will need to change: [1] to [kPa kPa-1]
                'RDRS_v2.1_P_HU_1.5m':  is_same, # [kg kg-1]
                'RDRS_v2.1_P_P0_SFC':   milibar_to_pascal, # [mb] * 100 = [Pa]
                'RDRS_v2.1_P_TT_1.5m':  celsius_to_kelvin, # [C] + 273.15 = [K]
                'RDRS_v2.1_P_UUC_10m':  kts_to_ms, # [kts] * 0.514 = [m s-1]
                'RDRS_v2.1_P_VVC_10m':  kts_to_ms, # [kts] * 0.514 = [m s-1]
                'RDRS_v2.1_P_UVC_10m':  kts_to_ms} # [kts] * 0.514 = [m s-1]
    return new_dict

def convert_rdrs_variables(data_in, rdrs_unit_conversion_dict):
    data_out = data_in.copy()
    
    # Convert values
    for variable, factor in rdrs_unit_conversion_dict.items():
        if variable == 'RDRS_v2.1_P_TT_1.5m':
            data_out[variable].values = data_in[variable].values + factor
        else:
            data_out[variable].values = data_in[variable].values * factor
    
    # Set the units
    data_out['RDRS_v2.1_A_PR0_SFC'].attrs['units'] = 'kg m**-2 s**-1'
    data_out['RDRS_v2.1_P_HR_1.5m'].attrs['units'] = 'kPa kPa**-1'
    data_out['RDRS_v2.1_P_P0_SFC'].attrs['units']  = 'Pa'
    data_out['RDRS_v2.1_P_TT_1.5m'].attrs['units'] = 'K'
    data_out['RDRS_v2.1_P_UUC_10m'].attrs['units'] = 'm s**-1'
    data_out['RDRS_v2.1_P_VVC_10m'].attrs['units'] = 'm s**-1'
    data_out['RDRS_v2.1_P_UVC_10m'].attrs['units'] = 'm s**-1'
    
    # Update the history
    now = time.strftime("%c")
    old_history = data_out.attrs['history']
    new_history = f'{old_history}. {now}: Converted units.'
    data_out.attrs['history'] = new_history
    
    return data_out

# Functions now in cs.geospatial
def derive_penman_monteith_pet(file, air_pressure='RDRS_v2.1_P_P0_SFC',
                                     air_temperature='RDRS_v2.1_P_TT_1.5m',
                                     relative_humidity='RDRS_v2.1_P_HR_1.5m',
                                     shortwave_radiation='RDRS_v2.1_P_FB_SFC',
                                     longwave_radiation='RDRS_v2.1_P_FI_SFC',
                                     wind_speed='RDRS_v2.1_P_UVC_10m',
                                     wind_height=10,
                                     ground_heat=0,
                                     new_name='PET',
                                     dims='rlat/rlon'):
    
    # 1. Get the values of this variable from the source (this automatically applies scaling and offset)
    airpres = file.variables[air_pressure][:]
    airtemp = file.variables[air_temperature][:]
    relhum = file.variables[relative_humidity][:]
    swdown = file.variables[shortwave_radiation][:]
    lwdown = file.variables[longwave_radiation][:]
    windspd = file.variables[wind_speed][:]
    
    # 2. Compute the new values
    values_new = calculate_penman_monteith_pet(airpres, airtemp, relhum, swdown, lwdown, windspd, windht=wind_height, groundheat=ground_heat)
    
    # 3. Define the other inputs
    unit_new = 'mm hr**-1'
    long_name = 'Potential evapotranspiration calculated with the Penman-Monteith method'
    new_history = f' On {time.ctime(time.time())}: Derived PET using Penman-Monteith method.'
    
    # 4. Create the new variable
    file = cs.make_nc_variable(file, new_name, unit_new, values_new, long_name=long_name, history=new_history, dims=dims)
    
    return file
        
def calculate_penman_monteith_pet(airpres, airtemp, relhum, swdown, lwdown, windspd, windht=10, groundheat=0):
    #
    # Calculation steps to find Penman-Monteith PET. The following is
    # mostly based on Zotarelli et al, 2009. Step by Step Calculation 
    # of the Penman-Monteith Evapotranspiration (FAO-56 Method).
    #
    # Zotarelli, L., Dukes, M. D., Romero, C. C., Migliaccio, K. W., and
    # Morgan, K. T.: Step by step calculation of the Penman-Monteith
    # Evapotranspiration (FAO-56 Method). University of Florida Extension,
    # AE459, available at: http://edis.ifas.ufl.edu (last access:
    # 21 January 2018), 10 pp., 2009.
    #
    # Main equation, modified for hourly data:
    # -----------------------------------------
    # PET    = (0.408 * delta * (Rn - G) + gamma * ((900/24) / (T + 273) * u2 * (es - ea))) / (delta + gamma * (1 + 0.34 * u2))
    #
    # PET    = potential evapotranspiration [mm/h]
    # delta  = slope of the saturation vapor pressure-temperature relationship [kPa/oC]
    # Rn     = net radiation at the crop surface [MJ/m2/h]
    # G      = soil heat flux density [MJ/m2/h]
    # gamma  = psychrometric constant [kPa/oC]
    # T      = air temperature at 2 m height [oC]
    # u2     = wind speed at 2 m height [m/s]
    # es     = saturation vapor pressure [kPa]
    # ea     = actual vapor pressure [kPa]
    #
    # Inputs
    # ------
    # airpres  = air pressure           (RDRS: RDRS_v2.1_P_P0_SFC)    [Pa]
    # airtemp  = 2-m air temperature    (RDRS: RDRS_v2.1_P_TT_1.5m)   [K]
    # relhum   = relative humidity      (RDRS: RDRS_v2.1_P_HR_1.5m)   [-]
    # swdown   = shortwave radiation    (RDRS: RDRS_v2.1_P_FB_SFC)    [W/m2]
    # lwdown   = longwave radiation     (RDRS: RDRS_v2.1_P_FI_SFC)    [W/m2]
    # windspd  = wind speed             (RDRS: RDRS_v2.1_P_UVC_10m)   [m/s]
    #
    # Optional inputs
    # ------------
    # windht   = wind height [m].           RDRS: 10 m
    # groundheat = soil heat flux [MJ/m2/h] Assumed: 0
    #
    # 0. Define constants
    con_Cp = 1.013 * 10**(-3)  # Specific heat of air at constant pressure [MJ kg-1 oC-1]
    con_eps = 0.622  # Ratio of the molecular weight of water vapor to dry air [-]
    con_lam = 2.45  # Latent heat of vaporization of water [MJ kg-1]
    con_albedo = 0.23  # Albedo
    con_sigmaSB = 5.67 * 10**(-8)  # Stefan-Boltzmann constant [W m-2 K-4]
    con_emissivity = 0.95  # Emissivity of the crop surface [-], still grass (Stull, Table 2-4)
    #
    # Convert units
    airpres = airpres / 1000  # [Pa] to [kPa]
    airtemp = airtemp - 273.15  # [K] to [oC]
    swdown = swdown * 3600 * 0.000001 # [W/m2] * [s hr-1] * [MJ J-1] = [MJ m-2 h-1]
    lwdown = lwdown * 3600 * 0.000001 # [W/m2] to [MJ m-2 h-1]
    #
    # 2. Convert wind speed to 2m
    u2 = windspd*(4.87 / np.log(67.8 * windht - 5.42)) # [m/s]
    #
    # 3. Calculate vapor pressure
    # saturation vapor pressure
    es = 0.6108 * np.exp((17.27 * airtemp) / (airtemp + 237.3)) # [kPa]
    #
    # We need actual vapor pressure. We can get this from dewpoint
    # temperature (Stull, 2018; Table 4-2b, Eq. 4.15), but calculating
    # dewpoint temperature is a bit cumbersome. We can also get it from
    # RDRS directly but this would mean redoing a lot of data processing,
    # because we initially didn't include dewpoint in the data subset.
    # Therefore, we will calculate it from relative humidity and saturation
    # vapor pressure. Other equations left here for reference.
    #
    # actual vapor pressure from dewpoint temperature
    #con_Rv_L = 1.844 * 10**(-4) # K-1
    #con_e0 = 0.6113  # kPa
    #con_T0 = 273.15  # K
    #dewtemp = (1/con_T0 - con_Rv_L * np.log(e/con_e0))**(-1) # Alternative: obtain directly from RDRS
    #ea = 0.6108 * np.exp((17.27 * dewtemp) / (dewtemp + 237.3)) # [kPa]
    #
    # actual vapor pressure from relative humidity
    ea = relhum * es # (Stull, 2018; Table 4-2b, Eq. 4.14a)
    #
    # 4. Delta and gamma
    # Psychometric constant
    gamma = con_Cp * airpres / (con_eps * con_lam) # [kPa/oC]
    #
    # Slope of the saturation vapor pressure-temperature relationship
    delta = (4098 * es) / (airtemp + 237.3)**2 # [kPa/oC]
    #
    # 5. Net radiation
    # Net short net shortwave radiation, [MJ m-2 hr-1] - Zea2009 eq 29
    Rns = (1 - con_albedo) * swdown
    #
    # Estimated outgoing longwave radiation
    # Original equation gives [W m-2], i.e. [J m-2 s-1]
    lwup = con_emissivity * con_sigmaSB * (airtemp + 273.15)**4
    lwup = lwup * 3600 * 0.000001 # [W/m2] to [MJ m-2 h-1]
    #
    # Net longwave radiation, [MJ m-2 hr-1]
    Rnl = lwdown - lwup
    #
    # Net radiation, [MJ m-2 hr-1]
    Rn = Rns - Rnl
    #
    # 6. Penman-Monteith PET [mm/hr]
    pet = (0.408*delta*(Rn-groundheat) + gamma*((900/24)/(airtemp+273))*u2*(es-ea)) / (delta + gamma*(1+0.34*u2))  
    #
    return pet

def make_nc_variable(file, name, units, values, long_name='', standard_name='', history='', dims='lat/lon'):
    #
    # Makes a netcdf4 variable with dimensions (time,latitude,longitude) in file
    #
    # Assumptions
    # - Dimension in input file are 'time', 'latitude', 'longitude'
    # - No fill value specified
    # - Compression options zlib=True, shuffle=True are active
    #
    # 1. Create the .nc variable
    if dims.lower() == 'lat/lon':
        file.createVariable(name,'f4',('time','latitude','longitude'), fill_value = False, zlib=True, shuffle=True)
    elif dims.lower() == 'hru':
        file.createVariable(name,'f4',('time','hru'), fill_value = False, zlib=True, shuffle=True)
    elif dims.lower() == 'rlat/rlon':
        file.createVariable(name,'f4',('time','rlat','rlon'), fill_value = False, zlib=True, shuffle=True)
    else:
        print(f'!!! Warning: make_nc_variable(): dimensions option {dims} not recognized. Exiting.')
        return file
    #
    # 2. Set the attributes FIRST, so we don't run into any scaling/offset issues
    file[name].setncattr('missing_value',-999)
    file[name].setncattr('units',units)
    if long_name: file[name].setncattr('long_name',long_name)
    if standard_name: file[name].setncattr('standard_name', standard_name)
    if dims.lower() == 'rlat/rlon': 
        file[name].setncattr('grid_mapping', 'rotated_pole') # matches RDRS existing var attributes
        file[name].setncattr('coordinates', 'lon lat') # add the secondary RDRS coordinate
    #
    # 3. Copy the data SECOND
    file[name][:] = values
    #
    # 4. Update history
    if history:
        # check if file has a History or history attribute
        if 'History' in file.ncattrs():
            current_history = file.getncattr('History')
            new_history = f'{current_history}{history}'
            file.setncattr('History', new_history)
        elif 'history' in file.ncattrs():
            current_history = file.getncattr('history')
            new_history = f'{current_history}{history}'
            file.setncattr('history', new_history)
    #
    return file
'''