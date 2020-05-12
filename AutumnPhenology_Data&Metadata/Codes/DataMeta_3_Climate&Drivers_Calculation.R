####################################################
## Climate and Soil Extraction & Drivers Calculation

# See Materials and Methods:
# Climate and soil data sets
# Summary of seasonal photosynthesis calculation

# Define directory paths
setwd(".../AutumnPhenology/AutumnPhenology_Data&Metadata/Data/")

# Load libraries
library(data.table)
library(readxl)
library(tidyverse)


##----------------------------------------
# Import data
##----------------------------------------

##----------------------------------------
## Spring phenology

# Import locations
stations <- fread("DataMeta_1_PhenologyObs_PEP725_Stations.csv")
stations <- stations %>% 
    select(PEP_ID,LON,LAT)

# Import dataset and select spring phenology
pheno.df <- fread("DataMeta_2_PhenologyObs_PEP725_CleanData.csv")
pheno_out.df <- pheno.df %>% 
    filter(phenology == "leaf.out") %>% 
    select(timeseries,PEP_ID,Species,YEAR,DoY)
pheno_out.df$ts_yr <- paste0(pheno_out.df$PEP_ID,"_",pheno_out.df$YEAR)

# Add coordinates
pheno_out.df <- merge(pheno_out.df,stations,by="PEP_ID")

# Data wrangling
pheno_out.df <- pheno_out.df %>% 
    rename(DoY_out = DoY) %>%
    select(timeseries, everything()) %>% 
    arrange(-desc(timeseries))


##----------------------------------------
## Soil texture

# Get coordinates
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
pheno_soil.df <- merge(pheno_out.df,soil_pars_ts.df,by="timeseries")
rm(list=setdiff(ls(), c("pheno_soil.df","main_folder")))

# Unload libraries
detach("package:rgdal", unload=TRUE)
detach("package:raster", unload=TRUE)
detach("package:soiltexture", unload=TRUE)
detach("package:sp", unload=TRUE)



##----------------------------------------
## CO2 concentrations

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

# AIC ETH dataset (1948-2014)
# Global and hemishperic means w/ monthly resolution
# Unit: Surface mole fractions (ppm)
CO2_ETH <- fread("mole_fraction_of_carbon_dioxide_in_air_input.csv")

# Mauna Loa dataset (2015)
# Unit: monthly atmospheric CO2 data 1958-present
CO2_mauna <- fread("monthly_in_situ_co2_mlo.csv") 

# Fill missing values
CO2_mauna[which(CO2_mauna$data_mean_nh==-99.99),]$data_mean_n <- CO2_mauna[which(CO2_mauna$data_mean_n==-99.99),]$CO2.1 

# Filter for years of interest (1948-2015)
years <- unique((pheno_soil_co2.df$YEAR)) 
years <- sort(years, decreasing=F)
years <- years[which(years>=1948 & years<=2014)]
CO2_ETH <- CO2_ETH %>% 
    filter(year %in% years)
CO2_mauna <- CO2_mauna %>% 
    filter(year == 2015)

# Select values for nothern hemisphere
CO2_ETH <- CO2_ETH %>% 
    select(year,month,data_mean_nh)
CO2_mauna <- CO2_mauna %>% 
    select(Yr,Mn,CO2)
CO2.df <- rbind(CO2_ETH,CO2_mauna)
colnames(CO2.df) <- c("YEAR","MONTH","CO2")

# Add monthly CO2 values to the phenology dataset
pheno_soil_co2.df <- data.frame()
for(yr in 1948:2015) {
  # Subset by year
  pheno.sub  <- pheno_soil.df %>% 
    filter(YEAR == yr)
  CO2.sub <- CO2.df %>% 
    filter(YEAR == yr)
  
  # Transpose CO2 monthly values to the phenological subset
  pheno.sub[as.character(1:12)] <- 0
  for(r in 1:nrow(pheno.sub)) {
    pheno.sub[r,as.character(1:12)] <- t(CO2.sub$CO2)
  }

  # Bind final dataset
  pheno_soil_co2.df <- rbind(pheno_soil_co2.df,pheno.sub)
  print(yr)
}

# Data wrangling 
pheno_soil_co2.df <- pheno_soil_co2.df %>% 
  arrange(timeseries,YEAR)
rm(list=setdiff(ls(), c("pheno_soil_co2.df", "model","main_folder")))


##----------------------------------------
## Extract climatic variables
# 1. NASA/GLDAS/V020/NOAH/G025/T3H --> Data availability (time): Jan 1, 1948 - Dec 31, 2010
# 2. NASA/GLDAS/V021/NOAH/G025/T3H --> Data availability (time): Jan 1, 2000 - Feb 14, 2019 

## Step 1: extract time-point of interest and bind all years for each climatic layer

# Sort years for each dataset
years <- unique((pheno_soil_co2.df$YEAR)) 
years <- sort(years, decreasing=F)
years_nas1 <- years[which(years>=1948 & years<=2010)] # years of NASA dataset 1
years_nas2 <- years[which(years>=2011 & years<=2015)] # years of NASA dataset 2

# Variable names for each datasets
var_nas1 <- c('SWnetsfc','LWnetsfc','Rainfsfc','SoilM0_10cm','SoilM10_40cm','TSoil0_10cm','Qairsfc','Tairsfc')
var_nas2 <- c('Swnet_tavg','Lwnet_tavg','Rainf_tavg','SoilMoi0_10cm_inst','SoilMoi10_40cm_inst','SoilTMP0_10cm_inst','Qair_f_inst','Tair_f_inst')

# Generate filenames
vn <- c('Net Short-wave Radiation','Net Long-wave Radiation','Precipitation','Soil Moisture 0-10 cm', 'Soil Moisture 10-40 cm', 'Soil temperature','Specific Humidity','Mean Temperature','Maximum Temperature','Minimum Temperature')
FileNames_var <- data.frame()
for (i in 1:length(vn)) { 
  fn_nas1 <- paste(var_nas1[i],"_",years_nas1,"_Extracted.csv",sep="")
  fn_nas2 <- paste(var_nas2[i],"_",years_nas2,"_Extracted.csv",sep="")
  fn_var <- c(fn_nas1,fn_nas2)
  FileNames_var[i,1:68] <- fn_var
  rownames(FileNames_var)[i] <- vn[i] 
}
FileNames_var["Maximum Temperature",] <- c(paste(var1[8],"_",years_nas1,"_Max_Extracted.csv",sep=""),paste(var_nas2[8],"_",years_nas2,"_Max_Extracted.csv",sep=""))
FileNames_var["Minimum Temperature",] <- c(paste(var1[8],"_",years_nas1,"_Min_Extracted.csv",sep=""),paste(var_nas2[8],"_",years_nas2,"_Min_Extracted.csv",sep=""))

# This function imports datasets of one variable for different years and binds them
for (x in vn) {
    # Apply the function for each variable/row 
    allData_var       <- lapply(FileNames_var[x==vn,], function(.file) { 
        dat           <- fread(paste(main_folder,.file,sep="")) 
        dat           <- as.data.frame(dat) 
        dat           <- dat[,-c(1,as.numeric(ncol(dat)))] # remove trash-info
        
        # Order columns according to dates
        dat_order     <- dat[,5:ncol(dat)][,order(parse_number(colnames(dat)[5:ncol(dat)]))]
        dat           <- cbind(dat[,1:4],dat_order)
        
        # Account for different number of days per year
        if (ncol(dat)==368) {
            names(dat) <- c('PEP_ID','LAT','LON',as.character(1:365)) 
        }
        if (ncol(dat)==369) {
            names(dat) <- c('PEP_ID','LAT','LON',as.character(1:366))  
        }

        # Add identifier column as YEAR
        if (x=="Maximum Temperature" | x=="Minimum Temperature") {
            ln        <- length(strsplit(.file,"_")[[1]])-2
        } 
        else {
            ln        <- length(strsplit(.file,"_")[[1]])-1
        }
        dat$YEAR      <- as.numeric(strsplit(.file,"_")[[1]][ln])
        dat           <- dat %>% select(YEAR, everything()) 

        # Add time-point identifier
        dat$ts_yr     <- paste0(dat$PEP_ID,"_",dat$YEAR)
        dat           <- dat %>% select(ts_yr, everything()) 
    })
    # Bind the yearly subsets in one single dataset for each variable 
    finalData_var     <- do.call(rbind.fill, allData_var)
    
    # Eliminate rows which do not match timeseries of phenology
    timeseries_id     <- unique(pheno_soil_co2.df$ts_yr)
    finalData_var     <- finalData_var[which(finalData_var$ts_yr %in% timeseries_id),]
    
    # Export datasets
    write.table(finalData_var, paste(x,".csv",sep=""),sep=";",row.names=FALSE)
    print(paste("---the",x,"dataset has been exported---"))
}

## Step 2: format units

# Precipitation estimates are given as rate in kg m-2 s-1. We need mm d-1
# 1 kg of rain water spread over 1 square meter of surface is 1 mm in thickness
# there are 60X60X24=86400 seconds in one day.
# Therefore, 1 kg m-2 s-1 = 86400 mm d-1
prec.df <- fread("Precipitation.csv")
prec.df[,as.numeric(1:366)] <- prec.df[,as.numeric(1:366)]*864000
write.table(prec.df,"Precipitation.csv",sep=";",row.names=FALSE)
print("---the Precipitation dataset has been converted---")

# GLDAS-2 gives temperatures in Celsius, whereas GLDAS-2.1 gives temperature in Kelvin. We want Celsius
library(weathermetrics)
vn <- c('Soil temperature','Mean Temperature','Maximum Temperature','Minimum Temperature')
years_nas2 <- 2011:2015 
DataList <- replicate(length(vn),data.frame())
for (i in 1:length(vn)) {
    dat <- fread(paste(vn[i],".csv",sep=""))
    DataList[[i]] <- dat
    DataList[[i]][which(DataList[[i]]$YEAR %in% years_nas2),][,as.numeric(1:366)] <- kelvin.to.celsius(DataList[[i]][which(DataList[[i]]$YEAR %in% years_nas2),][,as.numeric(1:366)], round=2)
    write.table(DataList[[i]],paste0(vn[i],".csv"),sep=";",row.names=FALSE)
    print(paste("---the ",vn[i]," dataset has been re-exported---"))
}
rm(list=setdiff(ls(), c("pheno_soil.df", "model","main_folder")))


##----------------------------------------
## Calculate drivers of leaf senescence
##----------------------------------------

##----------------------------------------
## Calculate photoperiod

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
write.table(photo.df,"Photoperiod.csv",sep=";",row.names=FALSE)
print("---the Photoperiod dataset has been exported---")


##----------------------------------------
## Calculate climatic predictors

# temp_GS: average temperature during the growing season
# temp_aut2: average minimum temperature of the 2 months before leaf senescence
# temp_aut3: average minimum temperature of the 3 months before leaf senescence
# RD: rainy days per year (i.e. number of days with > 2 mm of precipitation)
# RD_summer: rainy days during summer (i.e. water supply during the driest season)
# HRD: heavy rainy days per year (number of days with > 20 mm of precipitation)
# HRD_aut: heavy rainy days during autumn  
# H35: extreme heat events (i.e. number of days with maximum temperature > 35 degC)
# FD: frost days (number of days with minimum temperature < 0 degC)
# FD_spring: early frost events (i.e. 2 months after leaf flushing).

# Import climatic datasets 
vn <- c('Mean Temperature','Minimum Temperature','Maximum Temperature','Precipitation')
DataList <- replicate(length(vn),data.frame())
for(i in 1:length(vn)) {
  data <- fread(paste0(vn[i],".csv"))
  data$ts_yr <- paste0(data$PEP_ID,"_",data$YEAR)
  DataList[[i]] <- data
}
DataList[[5]] <- pheno_soil_co2.df
DataList[[6]] <- photo.df

# Initialize dataset 
Factors.df <- data.frame()

# Get all time-points
timeseries_year <- unique(pheno_soil_co2.df$ts_yr)

# Loop through all time-points
for(ty in timeseries_year) {
  
  # Subset by time-point
  pheno.sub <- DataList[[5]][which(DataList[[3]]$ts_yr==ty),]
  T_mean.sub <- DataList[[1]][which(DataList[[1]]$ts_yr==ty),]
  T_min.sub <- DataList[[2]][which(DataList[[2]]$ts_yr==ty),]
  T_max.sub <- DataList[[3]][which(DataList[[1]]$ts_yr==ty),]
  prec.sub <- DataList[[4]][which(DataList[[2]]$ts_yr==ty),]
  
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
  
  # Order time-series
  z  <- zoo(daily_vals, days)
  
  # Go from daily values to mean monthly values 
  month <- function(x)format(x, '%Y-%m')
  monthly_vals <- as.data.frame(aggregate(z, by=month, FUN=mean))

  # N.B.:
  # for each time-point there might be multiple rows (i.e. different species)
  species.n <- 1:nrow(pheno.sub)

  for(sp in species.n) {
    
    # Subset per species
    pheno.sub <- DataList[[5]][which(DataList[[3]]$ts_yr==ty),]
    pheno.sub <- pheno.sub[sp,]     
    factors.sub <- pheno.sub %>% 
      select(timeseries,ts_yr)
    photo.sub <- DataList[[6]][which(DataList[[6]]$lat_yr==pheno.sub$lat_yr),]
    
    # Growing season GS
    DoY_out <- pheno.sub$DOY_out
    DoY_off <- pheno.sub$DOY_off
    GS_interval <- DoY_out:DoY_off

    # Get the months of predicted leaf-out and leaf-off
    month_out <- lubridate::month(as.Date(DoY_out,origin=days[1]))
    month_off <- lubridate::month(as.Date(DoY_off,origin=days[1]))

    # Get the top-3 driest months during the growing season
    monthDrys <- sort(monthly_vals$Prec, index.return=TRUE, decreasing=FALSE)$ix
    monthDrys <- monthDrys[which(monthDrys<=month_off & monthDrys>=month_out)]
    monthDry3 <- monthDrys[1:3]
    
    # Calculate the average mean temperature of the growing season
    temp_GS <- T_mean.sub %>% 
      select(as.character(1:366)) %>% 
      select(as.character(GS_interval))
    factors.sub$temp_GS <- mean(as.numeric(temp_GS))
    
    # Calculate the average minimum temperature of the 2 months before leaf senescence
    temp_aut2 <- T_min.sub %>% 
      select(as.character(1:366)) %>% 
      select(as.character((DoY_off-60):DoY_off))
    factors.sub$temp_aut2 <- mean(as.numeric(temp_aut2))

    # Calculate the average minimum temperature of the 3 months before leaf senescence
    temp_aut3 <- T_min.sub %>% 
      select(as.character(1:366)) %>% 
      select(as.character((DoY_off-90):DoY_off))
    factors.sub$temp_aut3 <- mean(as.numeric(temp_aut3))

    # Calculate the number of rainy days of the whole year
    RD <- daily_vals %>% 
      filter(Prec>=2)
    factors.sub$RD <- nrow(RD)

    # Calculate the number of rainy days of during driest months
    RD_summer <- daily_vals %>% 
      filter(MONTH %in% monthDry3) %>% 
      filter(Prec>=2)
    factors.sub$RD_summer <- nrow(RD_summer)
    
    # Calculate the number of days with heavy rain for the whole year
    HRD <- daily_vals %>% 
      filter(Prec>=20)
    factors.sub$HRD <- nrow(HRD)
    

    # Calculate the number of days with maximum temperature > 35 degC
    H35 <- T_max.sub %>% 
      select(as.character(1:366)) %>% 
      select(as.character(GS_interval)) 
    factors.sub$H35 <- length(which(H35>35))
    
    # Calculate the number of days with minimum temperature < 0 degC
    FD <- T_min.sub %>% 
      select(as.character(1:366)) %>% 
      select(as.character(GS_interval)) 
    factors.sub$FD <- length(which(FD<0))

    # Calculate the early frost events
    FD_spring <- T_min.sub %>% 
      select(as.character(1:366)) %>% 
      select(as.character(DoY_out:(DoY_out+60))) 
    factors.sub$FD_spring <- length(which(FD_spring<0))
    
    print(paste0("RUN: ",pheno.sub$timeseries," => ",which(timeseries_year==ty)," OF ",length(timeseries_year)))
    Factors.df <- rbind(Factors.df,factors.sub)
  }
}


##----------------------------------------
## Calculate cGSI
# cGSI = cumulative growing season index

vn <- c('Maximum Temperature','Minimum Temperature','Mean Temperature')
DataList <- replicate(length(vn),data.frame())
for(i in 1:length(vn)) {
  data <- fread(paste0(vn[i],".csv"))
  data <- as.data.frame(data)
  DataList[[i]] <- data
}
DataList[[4]] <- pheno_soil_co2.df
DataList[[5]] <- photo.df

# Add the time-point identifiers
for(i in c(1:3)) {
  DataList[[i]]$ts_yr  <- paste0(DataList[[i]]$PEP_ID,"_",DataList[[i]]$YEAR)
}

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
GSI <- data.frame()

# Get all time-points
timeseries_year <- unique(DataList[[4]]$ts_yr)

for(ty in timeseries_year) {
  
  # Subest by time-point
  pheno.sub <- DataList[[4]][which(DataList[[4]]$ts_yr==ty),]
  T_max.sub <- DataList[[1]][which(DataList[[1]]$ts_yr==ty),]
  T_min.sub <- DataList[[2]][which(DataList[[2]]$ts_yr==ty),]
  T_mean.sub <- DataList[[3]][which(DataList[[3]]$ts_yr==ty),]
  
  # N.B.:
  # for each time-point there might be multiple rows (i.e. different species)
  species.n <- 1:nrow(pheno.sub)
  
  for(sp in species.n) {
    
    # Subset by species
    pheno.sub <- DataList[[4]][which(DataList[[4]]$ts_yr==ty),]
    pheno.sub <- pheno.sub[sp,]     
    GSI.sub <- pheno.sub %>% 
      select(timeseries,ts_yr,Species,PEP_ID,YEAR)
    
    # Subset photoperiod
    photo.sub         <- DataList[[5]][which(DataList[[5]]$lat_yr==pheno.sub$lat_yr),]
    photo.sub         <- as.data.frame(photo.sub)
    
    # Estimate VPD parameters based on plant-functional type (PFT)
    if(pheno.sub$PFT=="BNL") {
      VPD_min       <- 0.61
      VPD_max       <- 3.1
    } else {
      VPD_min       <- 1.1
      VPD_max       <- 3.6
    }
    
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
    
    # Subset climatic variables from the starting GS
    Tmax <- T_max.sub %>% 
      select(as.character(1:366)) %>% 
      select(as.character(GS_interval)) 
    Tmin <- T_min.sub %>% 
      select(as.character(1:366)) %>% 
      select(as.character(GS_interval)) 
    Tmean <- T_mean.sub %>% 
      select(as.character(1:366)) %>% 
      select(as.character(GS_interval)) 
    photoperiod <- photo.sub %>% 
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
    GSI <- rbind(GSI,GSI.sub)
    print(paste0("RUN: ",which(ty==timeseries_year)," OF ",length(timeseries_year)))
  }
}


##----------------------------------------
## Calculate photosynthesis
# cA_tot = cumulative net photosynthetic rate during the growing season, including and excluding a water deficit index
# cA_tot- = cumulative net photosynthetic rate during the growing season, excluding and excluding a water deficit index

# Import datasets
# DataList[[1]] = net short-wave radiation [W/m^2]
# DataList[[2]] = net long-wave radiation [W/m^2]
# DataList[[3]] = precipitation [mm]
# DataList[[4]] = mean temperature [degC]
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
cq                  <- 4.6e-6 # conversion factor for solar radiation at 550 nm from J/m2 to E/m2 (E=mol quanta)
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
photosynthesis.cum  <- data.frame()
photosynthesis_daily.cum  <- data.frame()

# Get timeseries
timeseries_year <- unique(DataList[[6]]$ts_yr)

for(ty in timeseries_year) {
  
  # Subset by time-point
  pheno.sub <- DataList[[6]][which(DataList[[6]]$ts_yr==ty),]
  shortw.sub <- DataList[[1]][which(DataList[[1]]$ts_yr==ty),]
  longw.sub <- DataList[[2]][which(DataList[[2]]$ts_yr==ty),]
  prec.sub <- DataList[[3]][which(DataList[[3]]$ts_yr==ty),]
  T_mean.sub <- DataList[[4]][which(DataList[[4]]$ts_yr==ty),]
  
  # N.B.:
  # for each time-point there might be multiple rows (i.e. different species)
  species.n <- 1:nrow(pheno.sub)
  
  for(sp in species.n) {
    
    # Subset by species
    pheno.sub <- DataList[[6]][which(DataList[[6]]$ts_yr==ty),]
    pheno.sub <- pheno.sub[sp,]     
    photosynthesis.sub <- pheno.sub %>% 
      select(timeseries,Species,PEP_ID,YEAR)
    
    # Subset photoperiod
    photo.sub <- DataList[[5]][which(DataList[[5]]$lat_yr==pheno.sub$lat_yr),]
    photo.sub <- as.data.frame(photo.sub)
    
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
    
    # Subset climatic variables from the starting GS
    ShortW <- shortw.sub %>% 
      select(as.character(1:366)) %>% 
      select(as.character(GS_interval)) 
    LongW <- longw.sub %>% 
      select(as.character(1:366)) %>% 
      select(as.character(GS_interval)) 
    Prec <- prec.sub %>% 
      select(as.character(1:366)) %>% 
      select(as.character(GS_interval)) 
    Tmean <- T_mean.sub %>% 
      select(as.character(1:366)) %>% 
      select(as.character(GS_interval)) 
    photoperiod <- photo.sub %>% 
      select(as.character(1:366)) %>% 
      select(as.character(GS_interval)) 
    
    ## Daily Net Photosynthesis rate (dA_n) and water stress factor (dw) are calculated daily and then accumulated by summation
    dA_tot_year <- vector()
    dA_gross_year <- vector()
    dR_year <- vector()
    dA_totw_year <- vector()
    
    for(day in GS_interval) {
      day <- as.character(day)
      
      ## Net photosynthesis rate (ref. Sitch et al. 2003)
      
      # apar: daily integral of absorbed photosynthetically active radiation (PAR), J m-2 d-1
      # Eqn 4, Haxeltine & Prentice 1996
      # alphaa: scaling factor for absorbed PAR at ecosystem, versus leaf, scale
      # nearly half of short-wave radiation is PAR --> mean annual value of 0.473 observed for the irradiance ratio in the PAR (ref. Papaioannou et al. 1993) plus 8% reflected and transmitted
      # convert in J/m^-2 day: the power in watts (W) is equal to the energy in joules (J), divided by the time period in seconds (s): 
      # --> 1 Watt = 1 Joule/second, therefore j = W*86400
      apar <- alphaa*ShortW[,day]*86400  
      
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
      
      # ca: intercellular partial pressure of CO_2, Pa
      # from Gt C per year to pCO2 yearly values for the growing season
      # first convert to ambient pressure, ppm (ca):
      # To convert from ppm to gigatonne of carbon:
      # the conversion tables of the Carbon Dioxide Information Analysis Center advise that: 
      # 1 part per million of atmospheric CO2 is equivalent to 2.13 Gigatonnes Carbon. 
      # Using the correction for emissions/concentrstions: 
      # 1ppm = 7.81 Gigatonnes of Carbon Dioxide.
      co2 <- pheno.sub %>% 
        select(CO2)
      ca <- 7.81*as.numeric(co2)
      
      # Convert ambient CO2 level, ca, from mole fraction to partial pressure, Pa
      pa <- ca*p
      
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
      vm <- (1.0/b)*(c1/c2)*((2.0*theta-1.0)*s-(2.0*theta*s-c2)*sigma)*apar*cmass*cq

      # je: PAR-limited photosynthesis rate, molC m-2 h-1
      # Eqn 3, Haxeltine & Prentice 1996
      # Convert je from daytime to hourly basis
      je <- c1*apar*cmass*cq/photoperiod[,day]
      
      # jc: rubisco-activity-limited photosynthesis rate, molC m-2 h-1
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
      if(pheno.sub$PFT=="TBL") { 
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
      if (pheno.sub$PFT=="TBL") { 
        root1 <- 0.7  
        root2 <- 0.3
      } else {
        root1 <- 0.9
        root2 <- 0.1
      }
      
      # melt: snowmelt
      # above -2 degC the snowpack begins to melt at a maximum rate of melt
      if(Tmean[,day] > -2){
        melt <- (Tmean[,day] + 2)*k_melt 
      }
      
      # w_max: soil texture-dependent difference between field capacity and wilting point, %
      w_max <- pheno.sub %>% 
        select(w_max)
      w_max <- as.numeric(w_max)
      
      # E_bare: evapotranspiration from bare soil
      E_bare <- E_pot*(root1*(pheno.sub$w1/w_max)) 
      
      # perc1: daily percolation from the upper to the lower soil layer
      w1.rel <- pheno.sub$w1/(pheno.sub$w1 + pheno.sub$w2)
      w2.rel <- pheno.sub$w2/(pheno.sub$w1 + pheno.sub$w2)
      perc1 <- pheno.sub$k_perc*w1.rel^2
      
      # perc2: daily percolation from the lower soil layer to the groundwater
      perc2 <- pheno.sub$k_perc*w2.rel^2
      
      # betas: rates of water extraction for the upper and the lower soil layers
      beta1 <- (root1*pheno.sub$w1)/(root1*pheno.sub$w1 + (1-root1)*pheno.sub$w2)
      beta2 <- 1-beta1
      
      # deltaw: daily changes in the water content of both layers; take into account:
      # INs: snowmelt, precipitation
      # OUTs: evapotranspiration, percolation, runoff
      if(is.na(pheno.sub$ETA)) {
        pheno.sub$ETA <- E_pot
      } else {
        pheno.sub$ETA <- pheno.sub$ETA
      }
      deltaw1 <- prec.sub[,day] + melt - E_bare - pheno.sub$ETA*beta1 - perc1
      deltaw2 <- perc1 + melt - pheno.sub$ETA*beta2 - perc2    
      
      # w1 = soil moisture for the upper layer (0-10cm)
      # w2 = soil moisture for the bottom layer (10-40cm)
      w1 <- pheno.sub$w1 + deltaw1
      w2 <- pheno.sub$w2 + deltaw2
      
      # Store values for the next day
      pheno.sub$w1 <- w1
      pheno.sub$w2 <- w2
      
      # wr: relative soil moisture wr
      # ratio between current soil water content and plant-available water capacity
      # wr is computed for both soil layers by weighting their relative soil water contents (w1, w2) with the fraction of roots present in the respective layer
      wr <- root1*(pheno.sub$w1/w_max) + root2*(pheno.sub$w2/w_max)
      
      # E_supply: plant- and soil-limited supply function 
      E_supply <- as.numeric(E_max*wr)
      
      # ETA: evapotranspiration 
      ETA <- min(E_supply,E_demand)
      pheno.sub$ETA <- ETA                 
      
      # dw: daily water stress factor
      dw <- min(1,(E_supply/E_demand))
      
      # dA_totw: daily net photosynthesis modified by water stress factor
      dA_totw <- adt*dw
      
      # Add daily result to the yearly vector
      dA_totw_year <- c(dA_totw_year,dA_totw)
    }
    
    # Calculate the cumulative values for the growing season
    cA_tot <- sum(dA_tot_year)
    cA_gross <- sum(dA_gross_year)
    cR_d <- sum(dR_year)
    cA_totw <- sum(dA_totw_year)
    
    # Store results of daily values
    photosynthesis_daily.sub <- photosynthesis.sub[rep(seq_len(nrow(photosynthesis.sub)),365),]
    photosynthesis_daily.sub$DoY <- 1:365 
    photosynthesis_daily.sub$dA_tot_year <- 0
    photosynthesis_daily.sub$dA_totw_year <- 0
    photosynthesis_daily.sub$dA_gross_year <- 0
    photosynthesis_daily.sub$dR_year <- 0
    photosynthesis_daily.sub[GS_interval,]$dA_tot_year <- dA_tot_year
    photosynthesis_daily.sub[GS_interval,]$dA_totw_year <- dA_totw_year
    photosynthesis_daily.sub[GS_interval,]$dA_gross_year <- dA_gross_year
    photosynthesis_daily.sub[GS_interval,]$dR_year <- dR_year
    
    # Store results of cumulative values
    photosynthesis.sub$cA_n <- cA_tot
    photosynthesis.sub$cA_nw <- cA_totw
    photosynthesis.sub$cA_gross <- cA_gross
    photosynthesis.sub$cR_d <- cR_d
    
    # Bind final datasets
    photosynthesis.cum  <- rbind(photosynthesis.cum,photosynthesis.sub)
    photosynthesis_daily.cum  <- rbind(photosynthesis_daily.cum,photosynthesis_daily.sub)
    
    print(paste0("RUN: ",ty," => ",which(timeseries_year==ty)," OF ",length(timeseries_year)))
    return(photosynthesis_daily.sub)
  }
}

# Export dataset
write.table(photosynthesis_daily.cum,"Photosynthesis_daily.csv",sep=";",row.names=FALSE)


##----------------------------------------
## Dataset of Autumn phenology drivers

# Add coordinates
stations <- fread("DataMeta_1_PhenologyObs_PEP725_Stations.csv")
stations <- stations %>% 
  select(PEP_ID,LON,LAT)

# Add time-spatial labels
pheno.df <- fread("DataMeta_2_PhenologyObs_PEP725_CleanData.csv")
drivers.df <- pheno.df %>% 
  filter(phenology == "leaf.out") %>% 
  select(timeseries,PEP_ID,Species,YEAR)
drivers.df <-  merge(drivers.df,stations,by="PEP_ID")

# Add phenological observations
DoY_off.df <- pheno.df %>% 
  filter(phenology == "leaf.off") %>% 
  select(DoY,mean_DoYoff)
drivers.df$DoY_off <- DoY_out.df$DoY_off 
DoY_out.df <- pheno.df %>% 
  filter(phenology == "leaf.out") %>% 
  select(DoY,mean_DoYout)
drivers.df$DoY_out <- DoY_out.df$DoY_out
drivers.df$autumn_anomaly <- DoY_off.df$DoY_off - DoY_off.df$mean_DoYoff
drivers.df$spring_anomaly <- DoY_out.df$DoY_out - DoY_out.df$mean_DoYout

# Add drivers of autumn phenology
drivers.df$temp_GS <- Factors.df$temp_GS
drivers.df$temp_aut2 <- Factors.df$temp_aut2
drivers.df$temp_aut3 <- Factors.df$temp_aut3
drivers.df$RD <- Factors.df$temp_RD
drivers.df$RD_summer <- Factors.df$RD_summer
drivers.df$HRD <- Factors.df$HRD
drivers.df$HD35 <- Factors.df$HD35
drivers.df$FD <- Factors.df$FD
drivers.df$FD_spring <- Factors.df$FD_spring
drivers.df$cGSI <- GSI$cGSI
drivers.df$`cA_tot-w` <- photosynthesis.cum$cA_tot
drivers.df$cA_tot <- photosynthesis.cum$cA_totw

# Data wrangling
drivers.df <- drivers.df %>% 
  select(timeseries, everything()) %>% 
  arrange(-desc(timeseries))

# Export dataset
write.table(drivers.df,"DataMeta_3_Drivers.csv",sep=";",row.names=FALSE)


##----------------------------------------
## References

# Haxeltine, A., & Prentice, I. C. BIOME3: An equilibrium terrestrial biosphere model based on ecophysiological constraints, resource availability, and competition among plant functional types. Glob Biogeochem Cy. 10, 693-709 (1996).
# Jolly, W. M., Nemani, R. & Running, S. W. A generalized, bioclimatic index to predict foliar phenology in response to climate. Glob. Chang. Biol. 11, 619-632 (2005).
# Smith B., Prentice I. C., & Sykes M. T. Representation of vegetation dynamics in the modelling of terrestrial ecosystems: comparing two contrasting approaches within European climate space. Glob Ecol Biogeogr. 10, 621-37 (2001).