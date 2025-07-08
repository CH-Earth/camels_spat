# CAMELS-spat
Repository that contains processing code used to generate the CAMELS-spat (Catchment Attributes and MEteorology for Large-Sample studies - SPATially distributed) data set.

For a description of the data set, processing steps and other hopefully helpful information, see here: https://egusphere.copernicus.org/preprints/2025/egusphere-2025-893/

To access the data, see here (and check the README on the data repository for download help, particularly about the sizes of various data set components and how to include/exclude these from downloads): https://dx.doi.org/10.20383/103.01216

## Issues and questions
If you are looking to use the CAMELS-SPAT data and want to know about potential pitfalls, please see the list of currently known issues here: https://github.com/CH-Earth/camels_spat/issues/42

If you are here to report a potential issue please create a new, separate post here: https://github.com/CH-Earth/camels_spat/issues

If you are here to ask a question about data set usage, please go to the Q&A section: https://github.com/CH-Earth/camels_spat/discussions/categories/q-a


## Repository description

This repository contains the following sub-folders:
- **0_config** - contains a configuration file with all high-level decisions. Prevents needing to hard-code paths in the remainder of the code and this in turns leads to more efficient reproducibility.
- **1_Python_setup** - contains: compiled `fiona` and `gdal` libraries for Windows 64-bit, Python package requirements and code to configure a Python virtual environment.
- **2_reference_shapefiles** - contains code to obtain and process reference shapefiles for the CAMELS-spat basins, obtained from the CAMELS-US data set and Water Survey of Canada data sets.
- **3_merit_hydro** - contains code to obtain MERIT Hydro flow accumulation and flow direction grids, and MERIT Hydro-derived shapefiles.
- **4_data_structure_prep** - contains code to generate a meta data file and prepare the necessary folder structures (makes folders and copies reference shapes if available).
- **5_basin_delineation** - contains code to subset downloaded MERIT Hydro-derived shapefiles to the basins of interest, and ensure that the shapefile extent aligns with the USGS and RHBN gauge locations.
- **6_flow_data** - contains code to obtain and process flow data from United States Geological Survey (USGS) and Water Survey of Canada (WSC) at daily and sub-daily time resolutions.
- **7_forcing_data** - contains code to obtain and process four forcing data sets: ERA5, EM-Earth, RDRS, and Daymet.
- **8_geospatial_data** - contains code to obtain and process 12 geospatial data sets, covering: climate (WorldClim), topography (MERIT Hydro), land cover (GLCLU2019 land cover and forest height, MODIS land cover and LAI, LGRIP30) surface water (HydroLAKES), soil (Pelletier, SOILGRIDS), geology (GLHYMPS).
- **9_attributes** - contains code to calculate basin attributes from the streamflow, forcing and geospatial data, at lumped and sub-basin resolutions.
- **10_analysis** - contains code for paper plots and various checks.
- **11_updates** - contains code for data set updates after initial publication.
