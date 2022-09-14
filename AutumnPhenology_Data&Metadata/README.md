# This folder contains all post-processed data of phenological observations, climatic layers and environmental drivers for the estimation of past and future autumn phenology, and related scripts

### The subfolder contains the following documents and scripts:
- Codes:
	- DataMeta_1_PEPdownload_cleaning.r: data cleaning of phenological observations (output: DataMeta_2_PhenologyObs_PEP725_CleanData.csv).
	- DataMeta_2_GLDAS_Download.js: downloading Global Land Data Assimilation System (GLDAS) climatic layers for 1948-2015 from Google Earth Engine (GEE). 
	- DataMeta_3_Add_drivers.r: extracting climatic variables and calculating drivers of autumn phenology of all timeseries for 1948-2015 (output: DataMeta_3_Drivers.csv).
	- DataMeta_4_FutureSpringProjections_3Models.r: estimating future leaf flushing dates during 2016-2100 using three spring phenology models: AT, M1 and PM1 (output: DataMeta_5_FutureSpringProjections_3Models.csv).
	- DataMeta_5_FutureClimate&Drivers_Calculation.r: extracting future climatic variables and calculating drivers of autumn phenology for 2016-2100.
