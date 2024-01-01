# Forcing data downloads
We use ERA5 (all variables) and deterministic EM-Earth (precipitation and temperature).

## ERA5
Code adapted from CWARHM: https://github.com/CH-Earth/CWARHM/

###  Time-invariant variables
#### Table
- Time-invariant variables (Table 1): https://confluence.ecmwf.int/pages/viewpage.action?pageId=82870405#ERA5:datadocumentation-Table1

#### Summary
- Geopotential [m^2 s^-2]
    - #9 in Table 1 
    - Description: https://codes.ecmwf.int/grib/param-db/?id=129
    - Purpose: needed for temperature lapse rate application
- ECMWF provides a larger set of time-invariant data. We download those as well because the data load is minimal, and these values might be useful for comparison purposes in certain studies. We thus download all of Table 1 linked above.

###  Time-variant variables
#### Tables
- Instantaneous surface level variables (Table 2): https://confluence.ecmwf.int/pages/viewpage.action?pageId=82870405#ERA5:datadocumentation-Table2
- Time-averaged surface level variables (Table 4): https://confluence.ecmwf.int/pages/viewpage.action?pageId=82870405#ERA5:datadocumentation-Table4
- Instantaneous pressure level variables (Table 9): https://confluence.ecmwf.int/pages/viewpage.action?pageId=82870405#ERA5:datadocumentation-Table9

#### Summary surface level variables
- Precipitation
    - #32 in Table 4
    - Variable name: `mean_total_precipitation_rate`
    - Description: https://apps.ecmwf.int/codes/grib/param-db?id=235055
- Downward shortwave radiation at the surface
    - #13 in Table 4
    - Variable name: `mean_surface_downward_short_wave_radiation_flux`
    - Description: https://apps.ecmwf.int/codes/grib/param-db?id=235035
- Net shortwave radiation at the surface
    - #15 in Table 4
    - Variable name: `mean_surface_net_short_wave_radiation_flux`
    - Description: https://apps.ecmwf.int/codes/grib/param-db?id=235037
- Downward longwave radiation at the surface
    - #14 in Table 4
    - Variable name: `mean_surface_downward_long_wave_radiation_flux`
    - Description: https://apps.ecmwf.int/codes/grib/param-db?id=235036
- Net longwave radiation at the surface
    - #16 in Table 4
    - Variable name: `mean_surface_net_long_wave_radiation_flux`
    - Description: https://apps.ecmwf.int/codes/grib/param-db?id=235038
- Air pressure
    - #39 in Table 2
    - Variable name: `surface_pressure`
    - Description: https://apps.ecmwf.int/codes/grib/param-db?id=134
- Potential evaporation (PET, see notes in paper)
    - #39 in Table 4
    - Variable name: `mean_potential_evaporation_rate`
    - Description: https://codes.ecmwf.int/grib/param-db/?id=235070
    - Additional details: https://codes.ecmwf.int/grib/param-db/?id=228251

#### Summary pressure level variables
-  Temperature
    - #5 in Table 9
    - Variable name: `temperature`
    - Description: https://apps.ecmwf.int/codes/grib/param-db?id=130
- U-component wind
    - #6 in Table 9
    - Variable name: `u_component_of_wind`
    - Description: https://apps.ecmwf.int/codes/grib/param-db?id=131
- V-component wind
    - #7 in Table 9
    - Variable name: `v_component_of_wind`
    - Description: https://apps.ecmwf.int/codes/grib/param-db?id=132
- Specific humidity
    - #8 in Table 9
    - Variable name: `specific_humidity`
    - Description: https://apps.ecmwf.int/codes/grib/param-db?id=133
- Relative humidity
    - #12 in Table 9
    - Variable name: `relative_humidity`
    - Description: https://apps.ecmwf.int/codes/grib/param-db?id=157

## EM-Earth
EM-Earth downloads rely on the Globus Transfer functionality. Initial downloads must be performed manually. See steps below.

Time-variant variables needed:
- Precipitation
    - [add details]
- Air temperature
    - [add details]
 
### Steps
1. Install the Globus Transfer client.
	- Install instructions (Linux): https://docs.globus.org/how-to/globus-connect-personal-linux/#installation
	- Install instructions (MacOS): https://docs.globus.org/how-to/globus-connect-personal-mac/#installation
	- Install instructions (Windows): https://docs.globus.org/how-to/globus-connect-personal-windows/#installation
2. Set up a local Globus Endpoint.
	- Setup instructions (Linux): https://docs.globus.org/how-to/globus-connect-personal-linux/#configuration
	- Setup instructions (MacOS): https://docs.globus.org/how-to/globus-connect-personal-mac/#configuration
	- Setup instructions (Windows): https://docs.globus.org/how-to/globus-connect-personal-windows/#configuration
3. Launch Globus Connect Personal
4. Navigate to the landing page of the EM-Earth data
	- Link: https://www.frdr-dfdr.ca/repo/dataset/8d30ab02-f2bd-4d05-ae43-11f4a387e5ad
5. Scroll down this page and click the `Download dataset` button. 
6. You will be prompted to log in to the Globus system. Do so.
7. You will be forwarded to a new webpage, where on the left side the EM-Earth folder structure is visible.
8. For precipitation downloads (~850GB):
	- Navigate to the folder with continental precipitation data, by clicking `EM_Earth_v1` > `deterministic_hourly` > `prcp`
	- Select, but do not open, the `NorthAmerica` folder
	- On the taskbar on the right, click `Transfer or Sync to..`
	- In the newly appeared search bar in the top right, find your personal Globus Endpoint (step 2)
	- In the file explorer that shows your personal machine, navigate to the folder where you want to keep the EM-Earth downloads. Note that you will need to specify this folder's location in the CAMELS-spat config file, so that the EM-Earth processing code knows where to look for the raw EM-Earth data downloads.
	- On the left side of the Globus file manager, click the `Start` button to initiate the transfer.
9. For the temperature downloads ():
	- Navigate to the folder with continental precipitation data, by clicking `EM_Earth_v1` > `deterministic_hourly` > `tmean`
	- Follow the instructions as given under (8)
