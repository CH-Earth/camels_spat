# Flow data downloads

## USGS metadata
We strip the header from USGS data downloads to ease the downloading process: instead of downloading a full file and modifying it later, we can now directly load the data into a Pandas dataframe. An example header is stored here for future reference:

URL request: 
```
https://nwis.waterservices.usgs.gov/nwis/iv/?format=rdb&sites=02395120&startDT=1950-01-01&endDT=2023-01-02&parameterCd=00060&siteStatus=all
```

Header + first two data lines:
```
# ---------------------------------- WARNING ----------------------------------------
# Some of the data that you have obtained from this U.S. Geological Survey database may not 
# have received Director's approval.  Any such data values are qualified as provisional and 
# are subject to revision.  Provisional data are released on the condition that neither the 
# USGS nor the United States Government may be held liable for any damages resulting from its use.
#  Go to http://help.waterdata.usgs.gov/policies/provisional-data-statement for more information.
#
# File-format description:  http://help.waterdata.usgs.gov/faq/about-tab-delimited-output
# Automated-retrieval info: http://help.waterdata.usgs.gov/faq/automated-retrievals
#
# Contact:   gs-w_support_nwisweb@usgs.gov
# retrieved: 2023-02-22 12:58:07 -05:00	(nadww01)
#
# Data-value grade codes included in this output:
#    91    IV verification DV <= 1 percent diff
#    92    IV verification DV <= 5 percent diff
#    93    IV verification DV <= 10 percent diff
#
# Data for the following 1 site(s) are contained in this file
#    USGS 02395120 TWO RUN CREEK NEAR KINGSTON, GA
# -----------------------------------------------------------------------------------
#
# TS_ID - An internal number representing a time series.
#
# Data provided for site 02395120
#    TS_ID       Parameter Description
#    40346       00060     Discharge, cubic feet per second
#
# Data-value qualification codes included in this output:
#     A  Approved for publication -- Processing and review completed.
#     e  Value has been estimated.
#
agency_cd	site_no	datetime	tz_cd	40346_00060	40346_00060_cd
5s	15s	20d	6s	14n	10s
USGS	02395120	1989-10-01 01:30	EDT	864	A:[91]
```