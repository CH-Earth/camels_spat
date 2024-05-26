# Functions to process RDRS and Daymet

import netCDF4 as nc
import numpy as np

def calculate_priestley_taylor_pet(jdoy, dayl, srad, tmax, tmin, wtvp, elev, lati):
    
    # Calculation steps to find Priestley-Taylor PET. The following is
    # based on Zotarelli et al, 2009. Step by Step Calculation of the
    # Penman-Monteith Evapotranspiration (FAO-56 Method).
    #
    # Zotarelli, L., Dukes, M. D., Romero, C. C., Migliaccio, K. W., and
    # Morgan, K. T.: Step by step calculation of the Penman-Monteith
    # Evapotranspiration (FAO-56 Method). University of Florida Extension,
    # AE459, available at: http://edis.ifas.ufl.edu (last access:
    # 21 January 2018), 10 pp., 2009.

    # PET    = aPT / lambda * delta / (delta + gamma) * (Rn-G)
    #
    # aPT    = Priestley-Taylor coefficient [-]
    # lambda = latent heat of vaporization of water [MJ/kg]
    # delta  = slope of the saturation vapor pressure-temperature
    #           relationship [kPa/oC]
    # gamma  = psychrometric constant [kPa/oC]
    # Rn     = net radiation [MJ/m2/d]
    # G      = soil heat flux [MJ/m2/d], assumed 0

    # Constants
    con_Gsc = 0.0820  # Solar constant [MJ/m2/min]
    con_albedo = 0.23  # Albedo
    con_sigma = 4.903e-9  # Stefan-Boltzmann constant [MJ/K4/m2/d]
    con_aPT = 1.26  # Priestley-Taylor coefficient
    con_lambda = 2.45  # Latent heat of vaporization [MJ/kg]

    # Convert units
    wtvp = wtvp / 1000  # [Pa] to [kPa]
    lati = np.deg2rad(lati)  # [degree] to [rad]

    # 1. delta: slope of the saturation vapor pressure-temperature relationship, [kPa/oC]
    # ---------------------------------------------------------------------

    # 1a. daily mean temperature [oC] - Zea2009 eq 5
    meanT = (tmax + tmin) / 2

    # 1b. delta [kPa/oC] - Zea2009 eq 9
    delta = (4098 * (0.6108 * np.exp((17.27 * meanT) / (meanT + 237.3)))) / ((meanT + 237.3) ** 2)

    # 2. gamma: psychrometric constant, [kPa/oC]
    # ---------------------------------------------------------------------

    # 2a. atmospheric pressure [kPa] - Zea2009 eq 10
    atmoP = 101.3 * (((293 - 0.0065 * elev) / 293) ** 5.26)

    # 2b. gamma, [kPa/oC] - Zea2009 eq 11
    gamma = 0.000665 * atmoP

    # 3. Rn: net radiation, [MJ/m2/d]
    # ---------------------------------------------------------------------

    # 3a. inverse relative distance Sun-Earth, [-] - Zea2009 eq 23
    dr = 1 + 0.033 * np.cos(2 * np.pi / 365 * jdoy)

    # 3b. solar declination [rad] - Zea2009 eq 24
    delta_sol = 0.409 * np.sin(2 * np.pi / 365 * jdoy - 1.39)

    # 3c. sunset hour angle, [rad] - Zea2009 eq 26
    ws_in = -np.tan(lati) * np.tan(delta_sol) # Input to np.arccos: this MUST be between -1 and 1 or np.arccos will return nan
    ws_in = ws_in.where(ws_in > -1, -1) # This is the case for high latitudes in summer (polar day, 24h light)
    ws_in = ws_in.where(ws_in <  1,  1) # High latitude winter (polar night, 24h dark)
    ws = np.arccos(ws_in)

    # 3d. extraterrestrial radiation, [MJ/kg/day] - Zea2009 eq 27
    Ra = (24 * 60 / np.pi) * con_Gsc * dr * (
        (ws * np.sin(lati) * np.sin(delta_sol)) + (np.cos(lati) * np.cos(delta_sol) * np.sin(ws))
    )

    # 3e. clear sky solar radiation, [MJ/m2/day] - Zea2009 eq 28
    Rs0 = (0.75 + 2e-5 * elev) * Ra

    # 3f. average daily total radiation, [MJ/m2/day]
    # daylight [s] * incident shortwave radiation [W/m2] = [J/m2/day]
    Rs = srad * dayl / 1e6

    # 3g. net solar or net shortwave radiation, [MJ/kg/day] - Zea2009 eq 29
    Rns = (1 - con_albedo) * Rs

    # 3h. net outgoing long wave solar radiation, [MJ/m2/day] - Zea2009 eq 30
    RsRs0Ratio = Rs / Rs0 # this can lead to division by zero, which is why we take it out of the equation
    RsRs0Ratio = RsRs0Ratio.where((Rs0 > 0) * (Rs >= 0), 0)
    Rnl = con_sigma * (((tmax + 273.16) ** 4 + (tmin + 273.16) ** 4) / 2) * (
        0.34 - 0.14 * np.sqrt(wtvp)
    ) * (1.35 * RsRs0Ratio - 0.35)

    # 3i. Rn, [MJ/m2/day] - Zea2009 eq 31
    Rn = Rns - Rnl

    # 4. Priestley-Taylor PET, [kg/m2/d]
    # assuming density water = 1000 kg/m3, this simplifies to [mm/d]
    # ---------------------------------------------------------------------
    pet = con_aPT / con_lambda * delta * Rn / (delta + gamma)
    pet = np.maximum(pet, 0)  # correct PET to 0 if P-T estimates negative PET due to Rnl > Rs (happens occasionally in high-latitude winters)

    return pet


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