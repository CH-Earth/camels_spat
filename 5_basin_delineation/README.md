# Basin delineation
Here we subset the downloaded MERIT Hydro shapefiles to the RHNB and CAMELS-US basins of interest. 

## Prepare metadata file
We'll track high-level info during the data preparation process in a single metadata file, stored as `.csv`. See table for descriptions.

The station (or basin outlet) locations provided by governing bodies don't always align with where the MERIT Hydro DEM thinks the river
network is. In most cases it is sufficient to map the provided station (`station_lat/lon`) or outlet (`outlet_lat/lon`) location onto the 
MERIT Hydro river network (`maaped_lat/lon`). In some cases this is not sufficient and we need to provide better coordinates for basin
delineation (`manual_lat/lon`). These are the result of an iterative procedure where we first delineated the basins using automated outlet
mapping and manually provided better locations for basins where this was needed.

| Field                                | Value(s)/example              | Description                                                                                          |
|--------------------------------------|-------------------------------|------------------------------------------------------------------------------------------------------|
| Country                              | CAN                           | Country the station is in                                                                            |
| Station_id                           | 01AD002                       | Station ID issued by responsible governing body; used as persistent identifier in CAMELS-spat        |
| Station_name                         | SAINT JOHN RIVER AT FORT KENT | Station name                                                                                         |
| Station_lat                          | 47.25806                      | Station latitude as defined by governing body                                                        |
| Station_lon                          | -68.59583                     | Station longitude as defined by governing body                                                       |
| Station_source                       | WSC 20222 data set            | Source that provided station lat/lon                                                                 |
| Outlet_lat                           | 47.25787655                   | Basin outlet latitde if different from station latitude, if specified by a data source               |
| Outlet_lon                           | -68.59491873                  | Basin outlet longitude if different from station latitude, if specified by a data source             |
| Outlet_source                        | WSC 20222 data set            | Source that provided basin outlet lat/lon                                                            |
| Mapped_lat                           | 47.25731224                   | Outlet latitude after automatic mapping of outlet or station lat onto the MERIT Hydro river network  |
| Mapped_lon                           | -68.59435442                  | Outlet longitude after automatic mapping of outlet or station lon onto the MERIT Hydro river network |
| Manual_lat                           | -999                          | Outlet latitude if specified in `2_manually_define_outlets.ipynb`                                    |
| Manual_lon                           | -999                          | Outlet longitude if specified in `2_manually_define_outlets.ipynb`                                   |
| Manual_outlet_location               |                               | Flag to indicate if manually specified the outlet (yes or [blank])                                   |
| Basin_area_km2                       | 14691.61894                   | Area of delineated basin [km^2]                                                                      |
| Ref_area_1_src                       | HYDAT gross drainage area     | First source that provides a reference basin area                                                    |
| Ref_area_1_km2                       | 14700                         | First reference basin area [km^2]; blank if none available                                           |
| Ref_area_2_src                       | HYDAT effective drainage area | Second source that provides a reference basin area                                                   |
| Ref_area_2_km2                       |                               | Second reference basin area [km^2]; blank if none available                                          |
| Ref_shape                            | yes                           | Flag to indicate if reference shape is available (yes/no)                                            |
| Ref_shape_source                     | WSC 2022 data set             | Source that provided the reference shape                                                             |
| Ref_shape_area_km2                   | 14677.4                       | Area of the reference shape in [km^2], if available                                                  |
| Ref_and_new_shape_fractional_overlap | 0.993791572                   | Fractional overlap of delineated basin shape and reference basin shape, if available                 |
| Manual_delineation_notes		       |                               | Reason for needing to specify a manual outlet location                                               |
| Delineation_confidence               | high                          | High: close match between delineation and reference data. Medium: Unsure if CAMELS-spat or reference data is more accurate. Low: reason to believe CAMELS-spat delineation is inaccurate |
| Delineation_notes                    |                               | Reason for delineation confidence statement                                                          |