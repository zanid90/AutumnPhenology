###########################################################
## Future Climate and Soil Extraction & Drivers Calculation

# See Materials and Methods:
# Climate and soil data sets
# Summary of seasonal photosynthesis calculation

# Define directory paths
setwd(".../AutumnPhenology/AutumnPhenology_Data&Metadata/Data/")

# Load libraries
library(data.table)
library(readxl)
library(tidyverse)

# Which model for future spring phenology?
model <- "M1" # AT, M1, PM1


##----------------------------------------
# Import data
##----------------------------------------

##----------------------------------------
## Future spring phenology

# Import dataset and select model of interest
wide_data <- fread("DataMeta_5_FutureSpringProjections_3Models.csv")
colnames(wide_data)
long_data <- wide_data %>% 
  filter(Model == model) %>% 
  select(-Scenario,-Model) %>% 
  pivot_longer(-c(timeseries,Species,PEP_ID,LAT,LON), values_to="DOY_out",names_to="Year")
long_data$Year <- as.numeric(gsub("[^0-9]", "", long_data$Year))


##----------------------------------------
## Soil texture

# Get coordinates
stations <- long_data %>% 
  select(timeseries,LON,LAT) 
stations <- stations[!duplicated(stations$timeseries),]
xy <- stations %>% 
  select(LON,LAT)

# Import raster files of soil texture
library(sp)
library(raster)
library(soiltexture)
library(rgdal)

# Soil images provided by ISRIC (World Soil Information)
# SoilGrids
# https://maps.isric.org/
Soils_coarse <- raster("SoilTexture_0cm.tif") # Coarse Fragments Volumetric in % at surface
Soils_fine <- raster("ClayContent_0cm.tif") # Fine "Clay" content Mass Fraction in % at surface

# Create a spatial object
spdf <- SpatialPointsDataFrame(coords = xy, data = stations)
head(coordinates(spdf))
class(spdf)

# Extract fragment values from rasters using coordinates 
proj4string(spdf) <- CRS("+init=epsg:4326")
Clay <- extract(Soils_coarse, spdf, cellnumbers = T)[,2]
Silt <- extract(Soils_fine, spdf, cellnumbers = T)[,2]
Sand <- 100-(Clay+Silt)
fragment.df <- data.frame(
  "CLAY" = Clay,
  "SILT" = Silt,
  "SAND" = Sand
)

# Get soil texture from fragment percentages
texture.df <- TT.points.in.classes(
  tri.data = fragment.df,
  class.sys = "HYPRES.TT",
  PiC.type = "l"
)

# Create dataframe of soil parameters and texture-to-parameter conversion
# Ref. Sitch, S. et al. Evaluation of ecosystem dynamics, plant geography and terrestrial carbon cycling in the LPJ dynamic global vegetation model. Glob. Chang. Biol. 9, 161-185 (2003).
soil_pars.df <- data.frame(
  "w_max" = c(13,13,21,14.8,10),
  "k_perc" = c(2,3,4,3.5,5),
  "c_soil" = c(0.1,0.1,0.02,0.02,0.2),
  row.names = c("VF","F","M","MF","C")
)

# Find soil parameters based on soil texture for each timeseries
soil_pars_ts.df <- data.frame()
for(row in 1:nrow(stations)) {
  sub.df <- texture.df[row,]
  index <- min(which(sub.df == TRUE))
  pars_ts <- soil_pars.df[index,]
  pars_ts$timeseries <- stations[row,]$timeseries
  soil_pars_ts.df <- rbind(soil_pars_ts.df,pars_ts)
  print(paste0("Soil parameters for ",stations[row,]$timeseries," executed!"))
}

# Add soil parameters to the phenology dataset
pheno_soil.df <- merge(long_data,soil_pars_ts.df,by="timeseries")
rm(list=setdiff(ls(), c("pheno_soil.df", "model","main_folder")))

# Unload libraries
detach("package:rgdal", unload=TRUE)
detach("package:raster", unload=TRUE)
detach("package:soiltexture", unload=TRUE)
detach("package:sp", unload=TRUE)


##----------------------------------------
## Future CO2 concentrations

# IPCC dataset (1900-2100)
# Fossil fuel and land use CO2 emissions and other well mixed GHG emissions follow IPCC's RCP8.5 (BAU)
# Global and hemishperic means w/ yearly resolution
# Unit: Gtons CO^2 per year
CO2_fut.df <- read_excel("2100-projections_climate-scoreboard_2015-1214.xlsx")
pheno_soil_co2.df <- CO2_fut.df %>% 
  filter(Year > 2015) %>% 
  select(Year,BAU) %>% 
  left_join(pheno_soil.df,., by="Year")
rm(CO2_fut.df,pheno_soil.df)


##----------------------------------------
## Get starting value for photosynthesis calculation

# w1 = soil moisture for the upper layer (0-10cm)
# w2 = soil moisture for the bottom layer (10-40cm)
# ETA = evapotranspiration
# For last day of the year before future projections:
# DOY = 364 & Year = 2015

# Get soil moisture for the upper (0-10cm) and lower (10-40cm) layers
vn <- c("SoilMoi0_10cm_inst","SoilMoi10_40cm_inst")
DataList <- replicate(2,data.frame())

for(i in 1:2){
  
  # Import the datasets per variable
  data <- fread(paste0(vn[i],"_2015_Extracted.csv"))
  
  # Pick the value for DOY = 364 and id (whihc corresponds to the PEP_ID)
  if(i == 1) {
    data <- data %>% 
      select(id,SoilMoi0_10cm_inst_364)
  }
  if(i == 2) {
    data <- data %>% 
      select(id,SoilMoi10_40cm_inst_364)
  }

  # Filter for the future sites
  data <- data %>% 
    filter(id %in% unique(pheno_soil_co2.df$PEP_ID))
  
  # Store in the datalist
  DataList[[i]] <- data
}

library(rlist)
start_val <- list.cbind(DataList)
start_val <- start_val[,-3] #eliminate replicate cols for id
colnames(start_val) <- c("PEP_ID","w1","w2")

# Add starting moisture values to the main dataset
pheno_soil_co2.df <- merge(pheno_soil_co2.df,start_val,by="PEP_ID")
pheno_soil_co2.df$ETA <- NA
colnames(pheno_soil_co2.df)

# Data wrangling 
pheno_soil_co2.df <- pheno_soil_co2.df %>% 
  rename(CO2 = BAU, YEAR = Year) %>% 
  select(timeseries, everything()) %>% 
  arrange(timeseries,YEAR)
rm(list=setdiff(ls(), c("pheno_soil_co2.df", "model","main_folder")))

# Export
setwd(main_folder)
write.table(pheno_soil_co2.df,paste0("Future_SpringPhenology_soil_CO2_",model,".csv"),sep=";",row.names=FALSE)


##----------------------------------------
## Extract future climatic variables

# Load libraries
library(raster)

# Set directory 
pheno_soil_co2.df <- fread("Future_SpringPhenology_soil_CO2_M1.csv")
future_data_folder <- "D:/FutureBeerData"
extract_folder <- paste0(future_data_folder,"/Extract")
setwd(extract_folder)

# Extract climatic values to dataframes
for (yr in 2016:2100) {
  
  # Subset the records data frame by year
  SubRecords <- subset(pheno_soil_co2.df,pheno_soil_co2.df$YEAR==yr)
  
  # Get the aggregated records' coordinate
  UniqueSites <- SubRecords[!duplicated(SubRecords$PEP_ID),]
  AggreCoordinates <- UniqueSites[,c("LON","LAT")]
  
  # Load the data
  TminStack <- brick(paste(future_data_folder,"/TmindegC/RCP85_0d50_",yr,"_stack.tif",sep=""))
  TmaxStack <- brick(paste(future_data_folder,"/TmaxdegC/RCP85_0d50_",yr,"_stack.tif",sep=""))
  TmeanStack <- brick(paste(future_data_folder,"/Mean_max_min/RCP85_Mean_max_min_",yr,".tif",sep=""))
  precStack <- brick(paste(future_data_folder,"/precip/RCP85_0d50_",yr,"_stack.tif",sep=""))
  SWradStack <- brick(paste(future_data_folder,"/SWdown/RCP85_0d50_",yr,"_stack.tif",sep=""))
  LWradStack <- brick(paste(future_data_folder,"/LWdown/RCP85_0d50_",yr,"_stack.tif",sep=""))
  
  # Extract the values by points' coordinate
  ExtractedTmin <- extract(TminStack,AggreCoordinates,method="simple")
  ExtractedTmax <- extract(TmaxStack,AggreCoordinates,method="simple")
  ExtractedTmean <- extract(TmeanStack,AggreCoordinates,method="simple")
  Extractedprec <- extract(precStack,AggreCoordinates,method="simple")
  ExtractedSWrad <- extract(SWradStack,AggreCoordinates,method="simple")
  ExtractedLWrad <- extract(LWradStack,AggreCoordinates,method="simple")
  
  # Transfer the output matrix to data frame 
  ExtractedTminData <- as.data.frame(ExtractedTmin)
  ExtractedTmaxData <- as.data.frame(ExtractedTmax)
  ExtractedTmeanData <- as.data.frame(ExtractedTmean)
  ExtractedprecData <- as.data.frame(Extractedprec)
  ExtractedSWradData <- as.data.frame(ExtractedSWrad)
  ExtractedLWradData <- as.data.frame(ExtractedLWrad)
  
  # cbind the coordinates to the extracted data
  ExtractedTminDataFrame <- cbind(UniqueSites$PEP_ID,AggreCoordinates,ExtractedTminData)
  ExtractedTmaxDataFrame <- cbind(UniqueSites$PEP_ID,AggreCoordinates,ExtractedTmaxData)
  ExtractedTmeanDataFrame <- cbind(UniqueSites$PEP_ID,AggreCoordinates,ExtractedTmeanData)
  ExtractedprecDataFrame <- cbind(UniqueSites$PEP_ID,AggreCoordinates,ExtractedprecData)
  ExtractedSWradDataFrame <- cbind(UniqueSites$PEP_ID,AggreCoordinates,ExtractedSWradData)
  ExtractedLWradDataFrame <- cbind(UniqueSites$PEP_ID,AggreCoordinates,ExtractedLWradData)
  
  # Write the data frame into the local directories
  write.csv(ExtractedTminDataFrame,paste("BeerDataFuture_TmindegC_", yr, "_Extracted.csv",sep=""))
  write.csv(ExtractedTmaxDataFrame,paste("BeerDataFuture_TmaxdegC_", yr, "_Extracted.csv",sep=""))
  write.csv(ExtractedTmeanDataFrame,paste("BeerDataFuture_TmeandegC_", yr, "_Extracted.csv",sep=""))
  write.csv(ExtractedprecDataFrame,paste("BeerDataFuture_prec_", yr, "_Extracted.csv",sep=""))
  write.csv(ExtractedSWradDataFrame,paste("BeerDataFuture_SWrad_", yr, "_Extracted.csv",sep=""))
  write.csv(ExtractedLWradDataFrame,paste("BeerDataFuture_LWrad_", yr, "_Extracted.csv",sep=""))
  
  # Show the progress of the extraction
  print(paste("------Climatic data for year ",yr," has already been extarcted and saved------",sep=""))
}

# Import datasets per variable for different years and bind them

# Set directory
future_data_folder <- paste0(main_folder,"/FutureBeerData")
setwd(future_data_folder)

# Variable names
var <- c("TmaxdegC","TmindegC","TmeandegC","prec","SWrad","LWrad")
vn <- c("Maximum Temperature","Minimum Temperature","Mean Temperature","Precipitation","Short-wave Radiation","Long-wave Radiation")

# Generate filenames
FileNames_var   <- data.frame()
for (i in 1:6) { 
  fn_var <- paste0("BeerDataFuture_",var[i],"_",as.character(2016:2100),"_Extracted.csv")
  FileNames_var[i,1:85]      <- fn_var
  rownames(FileNames_var)[i] <- vn[i] 
}

detach("package:raster", unload=TRUE)
for (x in vn) {
  allData_var     <- lapply(FileNames_var[x==vn,], function(.file) {
    
    # Load files
    dat           <- fread(paste("D:/FutureBeerData/Extract/",.file,sep=""))  
    dat           <- as.data.frame(dat) 
    dat           <- dat[,-1]
    
    # Account for different number of days per year
    if (ncol(dat)==368) {
      names(dat)  <- c('PEP_ID','LAT','LON',as.character(1:365)) 
    }
    if (ncol(dat)==369) {
      names(dat)  <- c('PEP_ID','LAT','LON',as.character(1:366))   
    }
    
    # add identifier column as YEAR
    dat$YEAR <- as.numeric(strsplit(.file,"_")[[1]][3])
    dat <- dat %>% select(YEAR, everything()) 
  })
  
  # Bind the yearly subsets in one single dataset for each variable & export them
  finalData_var     <- do.call(rbind.fill, allData_var)
  
  # Export datasets
  write.table(finalData_var, paste("Future ",x,".csv",sep=""),sep=";",row.names=FALSE)
  print(paste("---the",x,"dataset has been exported---"))
}
rm(list=setdiff(ls(), c("pheno_soil.df", "model","main_folder")))


##----------------------------------------
## Calculate predictors of future senescence
##----------------------------------------

##----------------------------------------
## Calculate photoperiod

# Set directory
setwd(main_folder)
setwd(extract_folder)

# Load libraries
library(geosphere)
library(zoo)
library(lubridate)

# Get all time-points
pheno_soil_co2.df$lat_yr <- paste0(pheno_soil_co2.df$LAT,"_",pheno_soil_co2.df$YEAR)
latitude_year <- unique(pheno_soil_co2.df$lat_yr)

# Initialize data frame to store results
photo.df <- data.frame()

for(ly in latitude_year) {
  
  # Subset table according to ty, i.e. location x year
  photo.sub <- pheno_soil_co2.df %>% 
    filter(lat_yr==ly) %>% 
    select(lat_yr,PEP_ID,LAT,YEAR)
  photo.sub <- photo.sub[!duplicated(photo.sub$lat_yr),]
  
  # Get the start and end date for the calculation 
  start_doy <- paste(photo.sub$YEAR,"-01-01",sep="") 
  end_doy <- paste(photo.sub$YEAR,"-12-31",sep="")
  
  # Define the time interval as the whole year
  doys <- seq(as.Date(start_doy), as.Date(end_doy), by="days")
  
  # Calculate daily photoperiod for the whole year
  photo <- daylength(photo.sub$LAT,doys)
  
  # Add daily photoperiods to the subset table
  photo.sub[as.character(1:366)] <- 0
  if(length(photo)==365) {
    photo.sub[,5:369] <- photo
  }
  if(length(photo)==366) {
    photo.sub[,5:370] <- photo
  }
  
  photo.df <- rbind(photo.df,photo.sub)
  print(paste0(ly," photoperiod calculated!"))
}

# Export dataset
write.table(photo.df,"FuturePhotoperiod.csv",sep=";",row.names=FALSE)
print("---the Future photoperiod dataset has been exported---")


##----------------------------------------
## Calculate future climatic predictors

# temp_GS: average temperature during the growing season
# RD_summer: rainy days during summer (i.e. water supply during the driest season)

# Set working directory
future_data_folder <- paste0(main_folder,"/FutureBeerData")
setwd(future_data_folder)
photo.df <- fread("FuturePhotoperiod.csv")

# Import datasets 
# ts_yr = location + year identifier
# DataList[[1]] = mean temperature [�C] --> daily values
# DataList[[2]] = precipiation [mm] --> daily values
# DataList[[3]] = future leaf.out dates, soil parameters and CO2 concentration
# DataList[[4]] = photoperiod [hours] --> daily values
vn <- c('Mean Temperature','Precipitation')
DataList <- replicate(2,data.frame())
for(i in 1:2) {
  data <- fread(paste0("Future ",vn[i],".csv"))
  data$ts_yr <- paste0(data$PEP_ID,"_",data$YEAR)
  DataList[[i]] <- data
}
names(DataList) <- vn
DataList[[3]] <- pheno_soil_co2.df
DataList[[4]] <- photo.df

# Initialize dataset 
Factors.df <- data.frame()

# Get all time-points
pheno_soil_co2.df$ts_yr <- paste0(pheno_soil_co2.df$PEP_ID,"_",pheno_soil_co2.df$YEAR)
timeseries_year <- unique(pheno_soil_co2.df$ts_yr)

# Loop through all time-points
for(ty in timeseries_year) {
  
  # Subset by time-point
  pheno.sub <- DataList[[3]][which(DataList[[3]]$ts_yr==ty),]
  T_mean.sub <- DataList[[1]][which(DataList[[1]]$ts_yr==ty),]
  prec.sub <- DataList[[2]][which(DataList[[2]]$ts_yr==ty),]
  photo.sub <- DataList[[6]][which(DataList[[6]]$lat_yr==pheno.sub$lat_yr),]
  
  # Generate sub-dataframe to store results
  factors.sub <- pheno.sub %>% 
    select(timeseries,ts_yr)
  
  # Define the current year in calendar units
  year <- strsplit(ty,"_")[[1]][2]
  start_doy <- paste(year,"-01-01", sep="") 
  end_doy <- paste(year,"-12-31", sep="")
  days <- seq(as.Date(start_doy), as.Date(end_doy), by="days")
  
  # Calculate the monthly averages
  TMEAN <- T_mean.sub %>% 
    select(as.character(1:366))
  PRCP <- prec.sub %>% 
    select(as.character(1:366))
  daily_vals <- data.frame(Tmean=as.numeric(TMEAN), Prec=as.numeric(PRCP), MONTH=0)
  daily_vals <- daily_vals[complete.cases(daily_vals),]
  daily_vals$MONTH <- lubridate::month(as.Date(days,origin=days[1]))
  
  # Convert precipitation to mm/day
  daily_vals$Prec <- daily_vals$Prec*864000 
  
  # Order time-series
  z  <- zoo(daily_vals, days)
  
  # Go from daily values to mean monthly values 
  month <- function(x)format(x, '%Y-%m')
  monthly_vals <- as.data.frame(aggregate(z, by=month, FUN=mean))

  # Calculate the growing season GS
  # Starting of GS is DoY_off (future projection)
  DoY_out <- pheno.sub$DOY_out
  
  # End of the GS is the first day below 11 hours after the beginning of the growing season
  endGS_site <- photo.sub %>% 
    select(as.character(1:366)) 
  endGS_site <- which(endGS_site<11)
  endGS_site <- endGS_site[which(endGS_site>DoY_out)][1]
 
  # Growing season
  GS_interval <- DoY_out:endGS_site
  
  # Calculate the mean temperature of the growing season
  temp_GS <- T_mean.sub %>% 
    select(as.character(1:366)) %>% 
    select(as.character(GS_interval))
  factors.sub$temp_GS <- mean(as.numeric(temp_GS))

  # Calculate RD_summer: number of rainy days (i.e. >2mm of precipitation)
  # during the driest season
  
  # Get the months of predicted leaf-out and leaf-off
  month_out <- lubridate::month(as.Date(DoY_out,origin=days[1]))
  month_off <- lubridate::month(as.Date(endGS_site,origin=days[1]))

  # Get the top-3 driest months during the growing season
  monthDrys <- sort(monthly_vals$Prec, index.return=TRUE, decreasing=FALSE)$ix
  monthDrys <- monthDrys[which(monthDrys<=month_off & monthDrys>=month_out)]
  monthDry3 <- monthDrys[1:3]
  
  # Get the rainy days of the driest months
  RD_summer <- daily_vals %>% 
    filter(MONTH %in% monthDry3) %>% 
    filter(Prec>=2)
  factors.sub$RD_summer <- nrow(RD_summer)
  
  print(paste0("RUN: ",pheno.sub$timeseries," => ",which(timeseries_year==ty)," OF ",length(timeseries_year)))
  Factors.df <- rbind(Factors.df,factors.sub)
}


##----------------------------------------
## Calculate cGSI

# Import datasets
# DataList[[1]] = maximum temperature [�C] --> daily values
# DataList[[2]] = minimum temperature [�C] --> daily values
# DataList[[3]] = mean temperature [�C] --> daily values
# DataList[[4]] = phenological data on leaf.out [Julian DAY]
# DataList[[5]] = photoperiod [hours] --> daily values

vn <- c('Maximum Temperature','Minimum Temperature','Mean Temperature')
DataList <- replicate(3,data.frame())
for(i in 1:3) {
  data <- fread(paste0("Future ",vn[i],".csv"))
  data <- as.data.frame(data)
  DataList[[i]] <- data
}
names(DataList) <- vn
DataList[[4]] <- pheno_soil_co2.df
DataList[[5]] <- fread("FuturePhotoperiod.csv")

# Add the time-point identifiers
for(i in c(1:4)) {
  DataList[[i]]$ts_yr  <- paste0(DataList[[i]]$PEP_ID,"_",DataList[[i]]$YEAR)
}

# Add unique id for phenological observations
DataList[[4]]$id  <- paste0(DataList[[4]]$PEP_ID,"_",DataList[[4]]$Species,"_",DataList[[4]]$YEAR)

# Add plant functional type label 
# T-BL-SG: Temperate broad-leaved summergreen tree
# B-NL-SG: Boreal needle-leaved summergreen tree
DataList[[4]]$PFT <- "TBL" 
DataList[[4]][which(DataList[[4]]$Species=="Larix decidua"),]$PFT <- "BNL"

## Helper functions
# This function converts degC into kPa
degC_to_kPa.fun <- function(temp) {
  out       <- 0.6108*exp((17.27*temp)/(temp+237.3)) 
  return(out)
}

# Temperature inhibition function from LPJ-GUESS 
temp_opt.fun <- function(temp) {
  x1        <- 1
  x2        <- 18
  x3        <- 25
  x4        <- 45
  k1        <- 2.*log((1/0.99)-1.)/(x1-x2)
  k2        <- (x1+x2)/2
  low       <- 1/(1+exp(k1*(k2-temp)))
  k3        <- log(0.99/0.01)/(x4-x3)
  high      <- 1-0.01*exp(k3*(temp-x3))
  tstress   <- low*high 
  if(tstress>=0) {
    tstress <- tstress
  } else {
    tstress <- 0
  }
  return(tstress)
}

# Photoperiod function
# photo = photoperiod 
# photo_min = minimum value during the growing season --> limited canopy development
# photo_max = maxmum value during the growing season --> allowed canopies to develop unconstrained
photoperiod.fun <- function(photo, photo_min, photo_max) {
  if(photo<=photo_min) {
    photo_resp    <- 0
  }
  if(photo<photo_max & photo>photo_min) {
    photo_resp    <- (photo-photo_min)/(photo_max-photo_min)
  }
  if(photo>=photo_max) {
    photo_resp   <- 1
  }
  return(photo_resp)
}

# Vapour Pressudre Deficit (VPD) function
# VPD = vapour pressure deficit [kPa]
# T_min & T_max = minimum and maximum daily temperature [�C]
# VPD_min --> at low values, latent heat losses are unlikely to exceed available water 
# little effect on stomata
# VPD_max --> at high values, particularly if sustained, photosynthesis and growth are likely to be significantly limited
# complete stomatal closure
VPD.fun <- function(VPD, VPD_min, VPD_max) {
  if(VPD>=VPD_max) {
    y             <- 0
  }
  if(VPD<VPD_max & VPD>VPD_min) {
    y             <- 1-((VPD-VPD_min)/(VPD_max-VPD_min))
  }
  if(VPD<=VPD_min) {
    y             <- 1
  }
  return(y)
}

## Cumulative Growing Season Index (cGSI) function
# (modified by Jolly et al. 2005)

# Initialize dataframes
GSI.df <- data.frame()

# Get all time-points per species
ids = unique(DataList[[4]]$id)

for(id_sub in ids) {
  
  # Subset input data by time-point
  pheno_sub.df = DataList[[4]] %>% 
    filter(id==id_sub)
  T_max_sub.df = DataList[[1]] %>% 
    filter(ts_yr==pheno_sub.df$ts_yr)
  T_min_sub.df = DataList[[2]] %>% 
    filter(Tts_yr==pheno_sub.df$ts_yr)
  T_mean_sub.df = DataList[[3]] %>% 
    filter(ts_yr==pheno_sub.df$ts_yr)
  photoperiod_sub.df = photoperiod.df %>% 
    filter(lat_yr==pheno_sub.df$lat_yr)
  
  # Generate sub-dataframe to store results
  GSI.sub <- pheno_sub.df %>% 
    select(timeseries,ts_yr,Species,PEP_ID,YEAR)
  
  # Estimate VPD parameters based on plant-functional type (PFT)
  if(pheno_sub.df$PFT=="BNL") {
    VPD_min       <- 0.61
    VPD_max       <- 3.1
  } else {
    VPD_min       <- 1.1
    VPD_max       <- 3.6
  }
  
  # Calculate the growing season GS
  # Starting of GS is DoY_off (future projection)
  DoY_out <- pheno_sub.df$DOY_out
  
  # End of the GS is the first day below 11 hours after the beginning of the growing season
  endGS_site <- photoperiod_sub.df %>% 
    select(as.character(1:366)) 
  endGS_site <- which(endGS_site<11)
  endGS_site <- endGS_site[which(endGS_site>DoY_out)][1]
  
  # Growing season
  GS_interval <- DoY_out:endGS_site
  
  # Subset climatic variables from the starting GS
  Tmax <- T_max_sub.df %>% 
    select(as.character(1:366)) %>% 
    select(as.character(GS_interval)) 
  Tmin <- T_min_sub.df %>% 
    select(as.character(1:366)) %>% 
    select(as.character(GS_interval)) 
  Tmean <- T_mean_sub.df %>% 
    select(as.character(1:366)) %>% 
    select(as.character(GS_interval)) 
  photoperiod <- photo_sub.df %>% 
    select(as.character(1:366)) %>% 
    select(as.character(GS_interval)) 
  
  # Estimate phoperiod thresholds based on the maximum and minimum values of the growing season
  photo_min <- min(photoperiod) 
  photo_max <- max(photoperiod) 
  
  # Initialize vector to store daily GSI
  iGSI_year <- 0
  
  # Loop through days in the growing season
  for(day in GS_interval) {
    day <- as.character(day)
    
    # e_s: saturation vapour pressure [kPa] 
    e_s <- (degC_to_kPa.fun(temp=Tmax[,day])+degC_to_kPa.fun(temp=Tmin[,day]))/2
    
    # e_a: derived from dewpoint temperature [kPa]
    e_a <- degC_to_kPa.fun(temp=Tmin[,day])
    
    # VPD: Vapour pressure deficit [kPa]
    VPD <- e_s-e_a
    iVPD <- VPD.fun(VPD, VPD_min, VPD_max)
    
    # iOpt_temp: response to optimal temperature (Gompertz function)
    iOpt <- temp_opt.fun(Tmean[,day])
    
    # iPhoto: photoperiod response
    iPhoto <- photoperiod.fun(photoperiod[,day], photo_min, photo_max)
    
    # Calculate daily GSI
    iGSI <- as.numeric(iVPD*iOpt*iPhoto)
    
    # Add to the cumulative cGSI
    iGSI_year <- c(iGSI_year,iGSI)
    cGSI <- sum(iGSI_year)
  }
  
  # Store results
  GSI.sub$cGSI <- cGSI
  
  # Bind final datasets
  GSI.df <- rbind(GSI.df,GSI.sub)
  print(paste0("COMPLETED: ",id_sub," => ",which(ids==id_sub)," OF ",length(ids)))
}


##----------------------------------------
## Calculate cA_tot

# Import datasets
# DataList[[1]] = net short-wave radiation [W/m^2]
# DataList[[2]] = net long-wave radiation [W/m^2]
# DataList[[3]] = precipitation [mm]
# DataList[[4]] = mean temperature [�C]
# DataList[[5]] = photoperiod [hours] 
# DataList[[6]] = phenology_soil_CO2 data [DAY for leaf.out, pCO2 and soil parameters] --> pCO2 are monthly values

vn  <- c('Short-wave Radiation','Long-wave Radiation','Precipitation','Mean Temperature')
DataList <- replicate(4,data.frame())
for(i in 1:4) {
  data <- fread(paste0("Future ",vn[i],".csv"))
  data <- as.data.frame(data)
  DataList[[i]] <- data
}
DataList[[5]] <- fread("FuturePhotoperiod.csv")
DataList[[6]] <- pheno_soil_co2.df

# Add the time-point identifiers
for(i in c(1:4,6)) {
  DataList[[i]]$ts_yr  <- paste0(DataList[[i]]$PEP_ID,"_",DataList[[i]]$YEAR)
}

# Add unique id for phenological predictions
DataList[[6]]$id  <- paste0(DataList[[4]]$PEP_ID,"_",DataList[[6]]$Species,"_",DataList[[6]]$YEAR)

# Add plant functional type label 
# T-BL-SG: Temperate broad-leaved summergreen tree
# B-NL-SG: Boreal needle-leaved summergreen tree
DataList[[6]]$PFT <- "TBL" 
DataList[[6]][which(DataList[[6]]$Species=="Larix decidua"),]$PFT <- "BNL"

## Helper functions

# Temperate inhibition function from LPJ-GUESS 
temp_opt.fun <- function(temp) {
  x1        <- 1
  x2        <- 18
  x3        <- 25
  x4        <- 45
  k1        <- 2.*log((1/0.99)-1.)/(x1-x2)
  k2        <- (x1+x2)/2
  low       <- 1/(1+exp(k1*(k2-temp)))
  k3        <- log(0.99/0.01)/(x4-x3)
  high      <- 1-0.01*exp(k3*(temp-x3))
  tstress   <- low*high 
  if(tstress>=0) {
    tstress <- tstress
  } else {
    tstress <- 0
  }
  return(tstress)
}

## Parameters for the Photosynthesis module
po2                 <- 20.9e3 #O2 partial pressure in Pa
p                   <- 1.0e5 # atmospheric pressure in Pa
bc3                 <- 0.015 # leaf respiration as fraction of Vmax for C3 plants
theta               <- 0.7 # colimitation (shape) parameter
q10ko               <- 1.2 #q10 for temperature-sensitive parameter ko
q10kc               <- 2.1 # q10 for temperature-sensitive parameter kc
q10tau              <- 0.57 # q10 for temperature-sensitive parameter tau
ko25                <- 3.0e4 # value of ko at 25 deg C
kc25                <- 30.0 # value of kc at 25 deg C
tau25               <- 2600.0 # value of tau at 25 deg C
alphaa              <- 0.5 # fraction of PAR assimilated at ecosystem level relative to leaf level
alphac3             <- 0.08 # intrinsic quantum efficiency of CO2 uptake in C3 plants
lambdamc3           <- 0.8 # optimal (maximum) lambda in C3 plants
cmass               <- 12.0 # atomic mass of carbon
cq                  <- 2.04e-6 # conversion factor for solar radiation from J m-2 to mol m-2
n0                  <- 7.15 # leaf N concentration (mg/g) not involved in photosynthesis
m                   <- 25.0 # corresponds to parameter p in Eqn 28, Haxeltine & Prentice 1996
t0c3                <- 250.0 # base temperature (K) in Arrhenius temperature response function for C3 plants
e0                  <- 308.56 # parameter in Arrhenius temp response function
tk25                <- 298.15 # 25 deg C in Kelvin
tmc3                <- 45.0 # maximum temperature for C3 photosynthesis
## Parameters for the Water balance module
gamma               <- 65 # psychrometer constant gamma [Pa/K]
L                   <- 2.5*10^6 # latent heat of vaporization of water L [J/kg]
emissivity          <- 0.6 # emissivity for coniferous and deciduous surface type
k_sb                <- 5.670367*10^-8 # Stefan-Boltzman constant [W/m^2 K^4]
d1                  <- 0.5 # thickness of upper soil layer [m]
d2                  <- 1 # thickness of lower soil layer [m]    
a_m                 <- 1.391 # maximum Priestley-Taylor coefficient a_m 
g_m                 <- 3.26 # scaling conductance g_m [mm/s]
k_melt              <- 3 # rate of snowmelt [mm/�C d]

## Soil parameters depending one texture [Phenologoy_CO2_soil dataset]
E_max               <- 5 # maximum transpiration rate that can be sustained under well-watered conditions E_max [mm/d] --> depends on plant functional type (same for T-BD-SG and B-NL-SG)
# w_max = soil texture-dependent difference between field capacity and wilting point w_max [%]
# c_soil = soil texture-dependent maximum rate of ETA from the bare soil [mm/h]
# k_perc = soil texture-dependent conductivity cond_soil or percolation rate field capacity [mm/d]

#PHOTOSYNTHESIS-CONDUCTANCE MODEL (modified by Sitch. et al. 2003)

# Initialize dataframe
photosynthesis.df  <- data.frame()

# Get all time-points per species
ids = unique(DataList[[6]]$id)

for(id_sub in ids) {
  
  # Subset input data by time-point
  pheno_sub.df = DataList[[6]] %>% 
    filter(id==id_sub)
  shortw_sub.df = DataList[[1]] %>% 
    filter(ts_yr==pheno_sub.df$ts_yr)
  longw_sub.df = DataList[[2]] %>% 
    filter(Tts_yr==pheno_sub.df$ts_yr)
  prec_sub.df = DataList[[3]] %>% 
    filter(Tts_yr==pheno_sub.df$ts_yr)
  T_mean_sub.df = DataList[[4]] %>% 
    filter(ts_yr==pheno_sub.df$ts_yr)
  photoperiod_sub.df = DataList[[5]] %>% 
    filter(lat_yr==pheno_sub.df$lat_yr)
  
  # Generate sub-dataframe to store results
  photosynthesis_sub.df <- pheno_sub.df %>% 
    select(timeseries,Species,PEP_ID,YEAR)
  
  # Calculate the growing season GS
  # Starting of GS is DoY_off (future projection)
  DoY_out <- pheno_sub.df$DOY_out
  
  # End of the GS is the first day below 11 hours after the beginning of the growing season
  endGS_site <- photoperiod_sub.df %>% 
    select(as.character(1:366)) 
  endGS_site <- which(endGS_site<11)
  endGS_site <- endGS_site[which(endGS_site>DoY_out)][1]
  
  # Growing season
  GS_interval <- DoY_out:endGS_site
  
  # Subset climatic variables from the starting GS
  ShortW <- shortw_sub.df %>% 
    select(as.character(1:366)) %>% 
    select(as.character(GS_interval)) 
  LongW <- longw_sub.df %>% 
    select(as.character(1:366)) %>% 
    select(as.character(GS_interval)) 
  Prec <- prec_sub.df %>% 
    select(as.character(1:366)) %>% 
    select(as.character(GS_interval)) 
  Tmean <- T_mean_sub.df %>% 
    select(as.character(1:366)) %>% 
    select(as.character(GS_interval)) 
  photoperiod <- photoperiod_sub.df %>% 
    select(as.character(1:366)) %>% 
    select(as.character(GS_interval)) 
  
  ## Daily Net Photosynthesis rate (dA_n) and water stress factor (dw) are calculated daily and then accumulated by summation
  dA_tot_year <- vector()
  dA_gross_year <- vector()
  dR_year <- vector()
  dA_totw_year <- vector()
  
  # Loop through days of the growing season
  for(day in GS_interval) {
    day <- as.character(day)
    
    ## Net photosynthesis rate (ref. Sitch et al. 2003)
    
    # apar: daily integral of absorbed photosynthetically active radiation (PAR), J m-2 d-1
    # Eqn 4, Haxeltine & Prentice 1996
    # alphaa: scaling factor for absorbed PAR at ecosystem, versus leaf, scale
    # nearly half of short-wave radiation is PAR --> mean annual value of 0.473 observed for the irradiance ratio in the PAR (ref. Papaioannou et al. 1993) plus 8% reflected and transmitted
    # convert in J/m^-2 day: the power in watts (W) is equal to the energy in joules (J), divided by the time period in seconds (s): 
    # --> 1 Watt = 1 Joule/second, therefore j = W*86400
    apar <- alphaa*ShortW[,day] * 60*60*24 
    
    # Calculate temperature inhibition function limiting photosynthesis at low and high temperatures (ref. Sitch et al. 2002)
    tstress <- temp_opt.fun(Tmean[,day])
    
    # Calculate catalytic capacity of rubisco, Vm, assuming optimal (non-water-stressed) value for lambda, i.e. lambdamc3
    # adjust kinetic parameters for their dependency on temperature 
    # i.e. relative change in the parameter for a 10 degC change in temperature
    # Eqn 22, Haxeltine & Prentice 1996a
    
    # ko: Michaelis constant of rubisco for O2
    ko <- ko25*q10ko**((Tmean[,day]-25.0)/10.0)
    
    # kc: Michaelis constant for CO2
    kc <- kc25*q10kc**((Tmean[,day]-25.0)/10.0) 
    
    # tau: CO2/O2 specificity ratio
    tau <- tau25*q10tau**((Tmean[,day]-25.0)/10.0) 
    
    # gammastar: CO_2 compensation point [CO2 partial pressure, Pa]   
    # Eqn 8, Haxeltine & Prentice 1996
    gammastar <- po2/(2.0*tau)
    
    # Get daily CO2, averaged by year [ppm]
    co2 <- pheno_sub.df %>% 
      select(CO2)
    
    # Convert ambient CO2 level from mole fraction to partial pressure, Pa
    pa <- co2*p
    
    # p_i: non-water-stressed intercellular CO2 partial pressure, Pa
    # Eqn 7, Haxeltine & Prentice 1996
    p_i <- pa*lambdamc3 
    
    # Calculate coefficients
    # Eqn 4, Haxeltine & Prentice 1996
    c1 <- tstress*alphac3*((p_i-gammastar)/(p_i+2.0*gammastar))
    
    # Eqn 6, Haxeltine & Prentice 1996
    c2 <- (p_i-gammastar)/(p_i+kc*(1.0+po2/ko)) 
    b <- bc3 # choose C3 value of b for Eqn 10, Haxeltine & Prentice 1996
    t0 <- t0c3 # base temperature for temperature response of rubisco
    
    # Eqn 13, Haxeltine & Prentice 1996
    s <- (24.0/photoperiod[,day])*b
    
    # Eqn 12, Haxeltine & Prentice 1996
    sigma <- sqrt(max(0.0,1.0-(c2-s)/(c2-theta*s)))
    
    # vm: optimal rubisco capacity, gC m-2 d-1 
    # Eqn 11, Haxeltine & Prentice 1996
    # cmass: the atomic weight of carbon, used in unit conversion from molC to g 
    # cq: conversion factor from apar [J m-2] to photosynthetic photon flux density [mol m-2]
    vm <- (1.0/b)*(c1/c2)*((2.0*theta-1.0)*s-(2.0*theta*s-c2)*sigma)*apar*cmass*cq
    
    # je: PAR-limited photosynthesis rate, gC m-2 h-1
    # Eqn 3, Haxeltine & Prentice 1996
    # Convert je from daytime to hourly basis
    je <- c1*apar*cmass*cq/photoperiod[,day]
    
    # jc: rubisco-activity-limited photosynthesis rate, gC m-2 h-1
    # Eqn 5, Haxeltine & Prentice 1996
    jc <- c2*vm/24.0
    
    # agd: daily gross photosynthesis, gC m-2 d-1
    # Eqn 2, modified with k_shape (theta)
    if(je<1e-10 | jc<=1e-10) {
      agd <- 0
    } else {
      agd <- (je+jc-sqrt((je+jc)**2.0-4.0*theta*je*jc))/(2.0*theta)*photoperiod[,day]
    }
    
    # rd: daily leaf respiration, gC m-2 d-1
    # Eqn 10, Haxeltine & Prentice 1996
    rd <- b*vm
    
    # and: daily net photosynthesis (at leaf level), gC m-2 d-1
    and <- agd-rd
    
    # adt: total daytime net photosynthesis, gC m-2 d-1
    # Eqn 19, Haxeltine & Prentice 1996
    adt <- and+(1.0-photoperiod[,day]/24.0)*rd
    
    # Convert adt from gC m-2 d-1 to mm m-2 d-1 using ideal gas equation
    adtmm <- adt/cmass*8.314*(Tmean[,day]+273.3)/p*1000.0
    
    # Store the daily result in the yearly vector
    dA_tot_year <- c(dA_tot_year,adt)
    dA_gross_year <- c(dA_gross_year,agd)
    dR_year <- c(dR_year,rd)
    
    
    ## Water Stress Factor (ref. Gerten et al. 2004)
    # soil is treated as a simple bucket consisting of two layers with fixed thickness
    
    # Calculate of potential evapotranspiration (ETA) rate, E_pot, mm d-1
    
    # delta: rate of increase of the saturation vapour pressure with temperature
    delta <- (2.503*10^6 * exp((17.269*Tmean[,day])/(237.3+Tmean[,day])))/(237.3+Tmean[,day])^2
    
    # R_n: istantaneous net radiation, W m-2 = R_s net short-wave radiation flux + R_l net long-wave flux 
    R_n <- ShortW[,day] + LongW[,day]
    
    # E_eq: equilibrium EvapoTranspiration
    # from seconds to day
    E_eq <- 24*3600*(delta/(delta+gamma))*(R_n/L) 
    
    # E_pot: potential EvapoTranspiration = equilibrium ETA * Priestley-Taylor coefficient 
    E_pot <- E_eq*a_m
    
    # ratio: stomata-controlled ratio between intercellular and ambient CO_2 partial pressure in the absence of water limitation
    ratio <- p_i/pa # ca. 0.8
    
    # g_min: minimum canopy conductance, mm s-1
    # depends on PFT
    if(pheno_sub.df$PFT=="TBL") { 
      g_min <- 0.5*3600*24 # from seconds to day
    } else {
      g_min <- 0.3*3600*24
    }
    
    # g_pot: nonwater-stressed potential canopy conductance, mm s-1
    g_pot <- g_min + ((1.6*adt)/((pa/p)*(1-ratio)))
    
    # E_demand: atmoshperic demand 
    # unstressed transpiration which occurs when stomatal opening is not limited by reduced water potential in the plant
    E_demand <- E_pot/(1+(g_m/g_pot))
    
    # root1/2: fraction of roots present in the respective layers
    # depends on PFT
    if (pheno_sub.df$PFT=="TBL") { 
      root1 <- 0.7  
      root2 <- 0.3
    } else {
      root1 <- 0.9
      root2 <- 0.1
    }
    
    #######
    # N.B.: we do not have predictions of soil moisture for the future
    # hence we calculate relative soil moisture wr by explicit simulation
    # of soil water processes for each day and update for the next
    #######
    
    # melt: snowmelt
    # above -2 degC the snowpack begins to melt at a maximum rate of melt
    if(Tmean[,day] > -2){
      melt <- (Tmean[,day] + 2)*k_melt 
    }
    
    # w_max: soil texture-dependent difference between field capacity and wilting point, %
    w_max <- pheno_sub.df %>% 
      select(w_max)
    w_max <- as.numeric(w_max)
    
    # E_bare: evapotranspiration from bare soil
    E_bare <- E_pot*(root1*(pheno_sub.df$w1/w_max)) 
    
    # perc1: daily percolation from the upper to the lower soil layer
    w1.rel <- pheno_sub.df$w1/(pheno_sub.df$w1 + pheno_sub.df$w2)
    w2.rel <- pheno_sub.df$w2/(pheno_sub.df$w1 + pheno_sub.df$w2)
    perc1 <- pheno_sub.df$k_perc*w1.rel^2
    
    # perc2: daily percolation from the lower soil layer to the groundwater
    perc2 <- pheno_sub.df$k_perc*w2.rel^2
    
    # betas: rates of water extraction for the upper and the lower soil layers
    beta1 <- (root1*pheno_sub.df$w1)/(root1*pheno_sub.df$w1 + (1-root1)*pheno_sub.df$w2)
    beta2 <- 1-beta1
    
    # deltaw: daily changes in the water content of both layers; take into account:
    # INs: snowmelt, precipitation
    # OUTs: evapotranspiration, percolation, runoff
    if(is.na(pheno_sub.df$ETA)) {
      pheno.sub$ETA <- E_pot
    } else {
      pheno_sub.df$ETA <- pheno_sub.df$ETA
    }
    deltaw1 <- Prec[,day] + melt - E_bare - pheno_sub.df$ETA*beta1 - perc1
    deltaw2 <- perc1 + melt - pheno_sub.df$ETA*beta2 - perc2    
    
    # w1 = soil moisture for the upper layer (0-10cm)
    # w2 = soil moisture for the bottom layer (10-40cm)
    w1 <- pheno_sub.df$w1 + deltaw1
    w2 <- pheno_sub.df$w2 + deltaw2
    
    # Store values for the next day
    pheno_sub.df$w1 <- w1
    pheno_sub.df$w2 <- w2
    
    # wr: relative soil moisture wr
    # ratio between current soil water content and plant-available water capacity
    # wr is computed for both soil layers by weighting their relative soil water contents (w1, w2) with the fraction of roots present in the respective layer
    wr <- root1*(pheno_sub.df$w1/w_max) + root2*(pheno_sub.df$w2/w_max)
    
    # E_supply: plant- and soil-limited supply function 
    E_supply <- as.numeric(E_max*wr)
    
    # ETA: evapotranspiration 
    ETA <- min(E_supply,E_demand)
    pheno_sub.df$ETA <- ETA                 
    
    # dw: daily water stress factor
    dw <- min(1,(E_supply/E_demand))
    
    # dA_totw: daily net photosynthesis modified by water stress factor
    dA_totw <- adt*dw
    
    # Add daily result to the yearly vector
    dA_totw_year <- c(dA_totw_year,dA_totw)
    
  } # END loop through days of the growing season
  
  # Calculate the cumulative values for the growing season
  cA_tot <- sum(dA_tot_year)
  cA_gross <- sum(dA_gross_year)
  cR_d <- sum(dR_year)
  cA_totw <- sum(dA_totw_year)
  
  # Store results of cumulative values
  photosynthesis_sub.df$cA_n <- cA_tot
  photosynthesis_sub.df$cA_nw <- cA_totw
  photosynthesis_sub.df$cA_gross <- cA_gross
  photosynthesis_sub.df$cR_d <- cR_d
  
  # Bind final datasets
  photosynthesis.df  <- rbind(photosynthesis.df,photosynthesis_sub.df)
  print(paste0("COMPLETED: ",id_sub," => ",which(ids==id_sub)," OF ",length(ids)))
}

# Export dataset
write.table(photosynthesis.df,"FuturePhotosynthesis.csv",sep=";",row.names=FALSE)
print("---the Future Photosynthesis rate has been exported---")


##----------------------------------------
## Dataset of Future Autumn phenology drivers

# Add time-spatial labels
pheno.df <- fread(paste0("Future_SpringPhenology_soil_CO2_",model,".csv"))
drivers.df <- pheno.df %>% 
  select(timeseries,Species,PEP_ID,YEAR)

# Add drivers of autumn phenology
drivers.df$temp_GS <- Factors.df$temp_GS
drivers.df$RD_summer <- Factors.df$RD_summer
drivers.df$cGSI <- GSI.df$cGSI
drivers.df$cA_tot <- photosynthesis.df$cA_totw

# Data wrangling
drivers.df <- drivers.df %>% 
  select(timeseries, everything()) %>% 
  arrange(-desc(timeseries))

# Export dataset
write.table(drivers.df,"DataMeta_6_FutureDrivers.csv",sep=";",row.names=FALSE)


##----------------------------------------
## References

# Haxeltine, A., & Prentice, I. C. BIOME3: An equilibrium terrestrial biosphere model based on ecophysiological constraints, resource availability, and competition among plant functional types. Glob Biogeochem Cy. 10, 693-709 (1996).
# Jolly, W. M., Nemani, R. & Running, S. W. A generalized, bioclimatic index to predict foliar phenology in response to climate. Glob. Chang. Biol. 11, 619-632 (2005).
# Smith B., Prentice I. C., & Sykes M. T. Representation of vegetation dynamics in the modelling of terrestrial ecosystems: comparing two contrasting approaches within European climate space. Glob Ecol Biogeogr. 10, 621-37 (2001).