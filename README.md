# CAMELS-spat
Repository that contains processing code used to generate the CAMELS-spat (Catchment Attributes and MEteorology for Large-Sample studies - SPATially distributed) data set.

## Repository description

This repository contains the following sub-folders:
- **0_config** - contains a configuration file with all high-level decisions. Prevents needing to hard-code paths in the remainder of the code and this in turns leads to more efficient reproducibility.
- **1_Python_setup** - contains:
	- **0_tools**: folder with shared functions.
	- **1_make_venv**: requirements of and code to configure a Python virtual environment.
- **2_reference_shapefiles** - contains code to obtain and process reference shapefiles for the CAMELS-spat basins, obtained from the CAMELS-US data set and Water Survey of Canada data sets.

	
## Reproducibility

To reproduce the data processing steps, execute scripts in order, starting at folder `1_Python_setup`, before moving on to main folder `2_`, main folder `3_`, etc. Folder names starting `0_` do not contain anything that needs to be executed manually. 

Further instructions and descriptions are found in the Readme's contain in sub-folders.