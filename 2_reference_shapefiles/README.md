# Reference shapefiles
For this project, we'll extend the existing CAMELS-US data (Newman et al., 2014, 2015) with basins from the Canadian Reference Hydrometric Basin Network (RHBN; Environment and Climate Change Canada, 2020). 

CAMELS-US shapefiles are available as a direct download: https://gdex.ucar.edu/dataset/camels/file.html

RHBN shapefiles need to be constructed as follows:
1. Download the list of RHBN stations through ftp: https://collaboration.cmc.ec.gc.ca/cmc/hydrometrics/www/RHBN/
2. Download two separate set of shapefiles:
	- Water Survey of Canada 2016: National hydrometric network basin polygons (Government of Canada, n.d.; Environment and Climate Change Canada, 2016)
	- Water Survey of Canada 2018: National Water Data Archive: HYDAT (Environment and Climate Change Canada, 2018; Water Survey of Canada, 2022)
3. Match the RHBN station metadata to the downloaded shapfiles and retain only the RHBN shapes.


## References
Environment and Climate Change Canada (2016). National hydrometric network basin polygons—Open Government Portal. https://open.canada.ca/data/en/dataset/0c121878-ac23-46f5-95df-eb9960753375. Last access: 2022-08-23

Environment and Climate Change Canada (2018, July 5). National Water Data Archive: HYDAT [Service description]. https://www.canada.ca/en/environment-climate-change/services/water-overview/quantity/monitoring/survey/data-products-services/national-archive-hydat.html. Last access: 2022-08-23

Environment and Climate Change Canada (2020, October 20). Reference Hydrometric Basin Network. https://www.canada.ca/en/environment-climate-change/services/water-overview/quantity/monitoring/survey/data-products-services/reference-hydrometric-basin-network.html. Last access: 2022-08-18

Government of Canada. (n.d.). National hydrometric network basin polygons—Environment and Climate Change Canada Data. Retrieved August 23, 2022, from https://donnees.ec.gc.ca/data/water/products/national-hydrometric-network-basin-polygons/?lang=en

Newman, A., Sampson, K., Clark, M., Bock, A., Viger, R., & Blodgett, D. (2014). A large-sample watershed-scale hydrometeorological dataset for the contiguous USA (p. approximately 2.5 GB) [Text/plain, text/tab-separated-values, png, shp]. UCAR/NCAR. https://doi.org/10.5065/D6MW2F4D

Newman, A. J., Clark, M. P., Sampson, K., Wood, A., Hay, L. E., Bock, A., Viger, R. J., Blodgett, D., Brekke, L., Arnold, J. R., Hopson, T., & Duan, Q. (2015). Development of a large-sample watershed-scale hydrometeorological data set for the contiguous USA: Data set characteristics and assessment of regional variability in hydrologic model performance. Hydrology and Earth System Sciences, 19(1), 209–223. https://doi.org/10.5194/hess-19-209-2015

Water Survey of Canada. (2022). National Hydrometric Network Basin Polygons [README]. https://collaboration.cmc.ec.gc.ca/cmc/hydrometrics/www/HydrometricNetworkBasinPolygons/Read_me.pdf