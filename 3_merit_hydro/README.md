# MERIT Hydro
We base the CAMELS-spat basins on the MERIT Hydro Digital Elevation Model (DEM; Yamazaki et al., 2019), and existing basin delineations using this DEM (Lin et al., 2019). Code is largely based on existing code (Knoben et al., 2022a, 2022b).

The code here downloadsand preprocesses where needed:
1. MERIT Hydro flow accumulation grid;
2. MERIT Hydro flow direction grid;
3. MERIT Hydro shapefiles.


## Note on catchment shapes
Original MERIT Hydro shapefiles were modified to more accurately discretize coastal hillslopes (Clark, unpublished) and are hence downloaded from a different source then the one mentioned in Lin et al. (2019).


## References
Knoben, W. J. M., Clark, M. P., Bales, J., Bennett, A., Gharari, S., Marsh, C. B., Nijssen, B., Pietroniro, A., Spiteri, R. J., Tang, G., Tarboton, D. G., & Wood, A. W. (2022a). Community Workflows to Advance Reproducibility in Hydrologic Modeling: Separating model-agnostic and model-specific configuration steps in applications of large-domain hydrologic models [Preprint]. Hydrology. https://doi.org/10.1002/essoar.10509195.2

Knoben, W. J. M., Marsh, C. B., & Tang, G. (2022b). CH-Earth/CWARHM: Initial release (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.6968609

Lin, P., Pan, M., Beck, H. E., Yang, Y., Yamazaki, D., Frasson, R., . . .  Wood, E. F.1662(2019).Global Reconstruction of Naturalized River Flows at 2.941663Million Reaches.Water Resources Research,55(8), 6499–6516. Retrieved 2021-166402-03, fromhttps://onlinelibrary.wiley.com/doi/10.1029/2019WR0252871665doi:  10.1029/2019WR025287

Yamazaki, D., Ikeshima, D., Sosa, J., Bates, P. D., Allen, G. H., & Pavelsky, T. M.1867(2019).MERIT Hydro:  A High-Resolution Global Hydrography Map1868Based on Latest Topography Dataset.Water Resources Research,55(6), 5053–18695073.Retrieved 2021-02-03, fromhttps://onlinelibrary.wiley.com/doi/187010.1029/2019WR024873doi:  10.1029/2019WR024873