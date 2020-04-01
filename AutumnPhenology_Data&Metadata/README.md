# This folder contains all post-processed data of phenological observations, climatic layers and environmental drivers for the estimation of past and future autumn phenology, and related scripts

### The two subfolders contain the following documents and scripts:
- Data:
	- DataMeta_1_PhenologyObs_PEP725_Stations.csv: locations of all phenological records included in the study. 
	- DataMeta_2_PhenologyObs_PEP725_CleanData.csv: leaf phenology records and summary statistics of all timeseries after filtering. 
	- DataMeta_3_Drivers.csv: drivers of autumn leaf phenology of all timeseries for 1948-2015. 
	- DataMeta_4_co2_soilParameters.csv: soil parameters and atmospheric CO2 concentration of all timeseries for 1948-2015. 
	- DataMeta_5_FutureSpringProjections_3Models.csv: future projections of spring leaf phenology for 2016-2100 using AT, M1 and PM1 models.
	- DataMeta_6_FutureDrivers.csv: drivers of autumn leaf phenology for 2016-2100.

- Codes:
	- DataMeta_1_PhenoData_Cleaning.r: data cleaning of phenological observations (output: DataMeta_2_PhenologyObs_PEP725_CleanData.csv).
	- DataMeta_2_GLDAS_Download.js: downloading Global Land Data Assimilation System (GLDAS) climatic layers for 1948-2015 from Google Earth Engine (GEE). 
	- DataMeta_3_Climate&Drivers_Calculation.r: extracting climatic variables and calculating drivers of autumn phenology of all timeseries for 1948-2015 (output: DataMeta_3_Drivers.csv and DataMeta_4_co2_soilParameters.csv).
	- DataMeta_4_FutureSpringProjections_3Models.r: estimating future leaf flushing dates during 2016-2100 using three spring phenology models: AT, M1 and PM1 (output: DataMeta_5_FutureSpringProjections_3Models.csv).
	- DataMeta_5_FutureClimate&Drivers_Calculation.r: extracting future climatic variables and calculating drivers of autumn phenology for 2016-2100 (output: DataMeta_6_FutureDrivers.csv).
