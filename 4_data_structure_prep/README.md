# Data structure preparation
Prepares the metadata file and the destination file structure where the CAMELS-spat data will go.

## Prepare metadata file
We'll track high-level info during the data preparation process in a single metadata file, stored as `.csv`. See table for descriptions.

Basin delineation from gauge locations does not work in certain cases (e.g. when the gauge is located on the river bank of a wide river)
and for these we manually need to define a more DEM-appropriate outlet location. These outlet locations are the result of a combination
of automated gauge-to-river mapping procedures (part of the delineation code) and manual iterative procedures for a small number of basins.
Manually found outlet locations are hard-coded when the metadata file is generated and used by the delineation code. Automatically found
outlet locations are added to the metadata file during basin delineation steps.

| Field            | Value(s)/example              | Description                                                                                  |
|------------------|-------------------------------|----------------------------------------------------------------------------------------------|
| Country          | CAN/US                        | Country the gauge is in                                                                      |
| Station_id       | 01AD002                       | Gauge ID issued by responsible governing body; used as persistent identifier in CAMELS-spat  |
| Station_name     |                               |                                                                                              |
| Station_lat      | <>                            | Gauge latitude as defined by governing body                                                  |
| Station_lon      | <>                            | Gauge longitude as defined by governing body                                                 |
| Outlet_lat       | <>                            | Longitude/latitude location used for basin delineation; result of iterative manual procedure |
| Outlet_lon       | <>                            | Longitude/latitude location used for basin delineation; result of iterative manual procedure |
| Basin_area_km2   | <>                            | Area of the delineated CAMELS-spat basin [km^2]                                              |
| Ref_area_1_src   | HYDAT gross drainage area     | Source of the first reference basin area                                                     |
| Reference area   | 14700                         | Area of the first reference shape in [km^2], if available                                    |
| Ref_area_2_src   | HYDAT effective drainage area | Source of the second reference basin area                                                    |
| Reference area   |                               | Area of the second reference shape in [km^2], if available                                   |
| Reference shape  |                               | Flag to indicate if reference shape is available                                             |
| Reference area   |                               | Area of the reference shape in [km^2], if available                                          |

## Data structure prep
Creates an output folder structure. Copies reference shape to destination folder and stores reference shape area in the metadata file if a reference shape is available.