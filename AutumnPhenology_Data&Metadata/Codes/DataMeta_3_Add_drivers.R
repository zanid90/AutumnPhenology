#---
#title: Increased growing-season productivity drives earlier autumn leaf senescence in temperate trees (Zani et al. 2020 Science)
#author: Constantin Zohner, Deborah Zani
#date: "last updated May 17, 2022"

#R code generating the drivers



### Driver names
#- Anet: Net day-time photosynthesis
#- GSI: Growing-season index
#- Tls: Leafy season mean temperature (Liu et al. 2018)
#- LPI: Low precipitation index
#- VPD: Vapour pressure deficit index
#- Tautumn: Autumn temperature



#####################
# Required packages #
#####################



require(data.table)
require(ncdf4)
require(raster)
require(tidyverse)
require(sp)
require(rpmodel)
require(purrr)
require(pbmcapply)
require(zoo)
require(lubridate)
require(weathermetrics)



#######################################################################################################################
#######################################################################################################################



################
# Define Paths #
################


# 1. Input
##########

#Climate
GLDAS_path       = "/Users/consti/Desktop/PhD/Publication_material/17_Autumn_phenology_tier2/PEP_analysis/Analysis/Analysis_input/Drivers/GLDAS"

#Photoperiod
clim_path        = "/Users/consti/Desktop/PhD/Publication_material/17_Autumn_phenology_tier2/PEP_analysis/Analysis/Analysis_input/Drivers"

#Soil
SoilTexture_path = "/Users/consti/Desktop/PhD/Publication_material/17_Autumn_phenology_tier2/PEP_analysis/Analysis/Analysis_input/Drivers/SoilTexture"

#Phenology and drivers
Data_path        = "/Users/consti/Desktop/PhD/Publication_material/000_Zani_et_al_reanalysis/2_Combined_analysis/Analysis_input"


# 2. Output
###########

PEP_drivers_files   = "/Users/consti/Desktop/PhD/Publication_material/000_Zani_et_al_reanalysis/2_Combined_analysis/Analysis_input/Driver_files/Individual_files"
PEP_drivers_final   = "/Users/consti/Desktop/PhD/Publication_material/000_Zani_et_al_reanalysis/2_Combined_analysis/Analysis_input"
PEP_drivers_missing = "/Users/consti/Desktop/PhD/Publication_material/000_Zani_et_al_reanalysis/2_Combined_analysis/Analysis_input/Driver_files/Missing_observations"



#######################################################################################################################
#######################################################################################################################



#################
## Import data ##
#################



## Phenology data
#################

PEP.df = fread(paste(Data_path, "DataMeta_2_PhenologyObs_PEP725_CleanData.csv", sep="/"))%>%
  #order table
  arrange(species, pep_id, year) %>% 
  dplyr::select(-c(V1)) 


## LAI data
###########

LAI.df = fread(paste(Data_path, "LAI_data.csv", sep="/")) %>% 
  dplyr::select(-c(V1))
PEP.df = merge(PEP.df, LAI.df, by=c("pep_id","species","timeseries","year"))#merge with PEP725 data


## CO2 data
###########

CO2.df = fread(paste(Data_path, "CO2_data.csv", sep="/"))


## Photoperiod
##############

photo.df = fread(paste(clim_path, "Photoperiod.csv", sep="/"))


## Soil Texture
###############

SoilTexture.df = fread(paste(SoilTexture_path, "SoilTexture.csv", sep="/"))


## Import daily climatic datasets from GLDAS
############################################

#define climate variables
vn <- c('Daily_Mean_Data_Tair_f_inst','Daily_Min_Data_Tair_f_inst','Daily_Max_Data_Tair_f_inst',
        'Daily_Data_Rainf_f_tavg','Daily_Data_Qair_f_inst',
        'Daily_Data_SoilMoi0_10cm_inst','Daily_Data_SoilMoi10_40cm_inst',
        'Daily_Data_Swnet_tavg','Daily_Data_Lwnet_tavg','Daily_Data_SWdown_f_tavg')
#create empty list
DataList <- replicate(length(vn),data.frame())
#loop through climate variables
for(i in 1:length(vn)) {
  #read data
  data = fread(paste0(GLDAS_path, "/", vn[i],".csv")) %>%
    #add site x year identifier
    mutate(site_year = paste0(pep_id, '_', year)) %>%
    #order table
    dplyr::select(site_year, pep_id, year, lat, lon, everything()) 
  #add table to list
  DataList[[i]] <- data
}
#add names to list
names(DataList)=vn
# Note: Precipitation is given as rate in mm d-1.


## Get average monthly precipitation per site
#############################################

PrcpMonths.df = DataList[[4]] %>% 
  #long format
  pivot_longer(., -c(site_year, pep_id, year, lat, lon), names_to = "DOY", values_to = "prcp") %>% 
  #delete NAs
  filter(!is.na(prcp)) %>% 
  #add date and month info
  mutate(date  = as.Date(as.numeric(DOY), origin = paste0(year-1, "-12-31")),
         Month = lubridate::month(date)) %>% 
  #sum precipitation for each month, site and year
  group_by(pep_id,year,Month) %>% 
  summarize(prcp_sum = sum(prcp)) %>% 
  ungroup() %>% 
  #get the mean for each site and month
  group_by(pep_id,Month) %>% 
  summarize(prcp_mean = mean(prcp_sum)) %>% 
  ungroup()



#######################################################################################################################
#######################################################################################################################



##############################
## Add Soil texture and PFT ##
##############################



#Add w_max (soil texture-dependent difference between field capacity and wilting point [%])
PEP.df = merge(PEP.df, SoilTexture.df[,c('pep_id','w_max')], by='pep_id')

#Add plant functional type info
PEP.df$PFT <- "TBL" # T-BL-SG: Temperate broad-leaved summergreen tree
PEP.df[PEP.df$species=='Larix',]$PFT <- "BNL" # B-NL-SG: Boreal needle-leaved summergreen tree

#remove stuff
rm(SoilTexture.df)



#######################################################################################################################
#######################################################################################################################



###############
## Constants ##
###############



## Constants in the Photosynthesis module
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
cmass               <- 12.0107 # molecular mass of C [g mol-1]
cq                  <- 2.04e-6 # accounts for the energy content of incoming short-wave radiation and the fraction of photosynthetically active radiation in total shortwave radiation; conversion factor for solar radiation from J m-2 to mol m-2 (Meek et al. 1984, Stocker et al. 2020)
n0                  <- 7.15 # leaf N concentration (mg/g) not involved in photosynthesis
m                   <- 25.0 # corresponds to parameter p in Eqn 28, Haxeltine & Prentice 1996
t0c3                <- 250.0 # base temperature (K) in Arrhenius temperature response function for C3 plants
e0                  <- 308.56 # parameter in Arrhenius temp response function
tk25                <- 298.15 # 25 deg C in Kelvin
tmc3                <- 45.0 # maximum temperature for C3 photosynthesis

## Constants in the Water balance module
gamma               <- 65 # psychrometer constant gamma [Pa/K]
L                   <- 2.5*10^6 # latent heat of vaporization of water L [J/kg]
d1                  <- 0.5 # thickness of upper soil layer [m]
d2                  <- 1 # thickness of lower soil layer [m]    
a_m                 <- 1.391 # maximum Priestley-Taylor coefficient a_m 
g_m                 <- 3.26 # scaling conductance g_m [mm/s]
k_melt              <- 3 # rate of snowmelt [mm/???C d]
E_max               <- 5 # maximum transpiration rate that can be sustained under well-watered conditions E_max [mm/d]



#######################################################################################################################
#######################################################################################################################



######################
## Helper functions ##
######################



##############################################
# Temperate inhibition function from LPJ-GUESS 
##############################################

temp_opt.fun <- function(temp) {
  x1        <- 5
  x2        <- 18
  x3        <- 25
  x4        <- 45
  k1        <- 2*log((1/0.99)-1)/(x1-x2)
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

#Response plot
plot(temp_opt.fun(c(0:45)), type="b", xlab=expression("Temperature"~(degree*C)), ylab="Photosynthessis response", col="red")


######################
# Photoperiod function
######################

# Photosynthetic capacity (Vcmax, Jmax) declines after summer solstice (Bauerle et al. 2012) 
# photo = photoperiod 
# photo_min = day length at which photosynthetic capacity reaches zero. 
# photo_max = maximum day length during the growing season --> maximum photosynthetic capacity
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

#Example with maximum day-length of 16 hours
photoperiod = c(8,11,16)
response=sapply(photoperiod, photoperiod.fun, photo_min=11, photo_max=16)
plot(response~photoperiod, type="b", xlab="Daylength (h)", ylab="Photosynthessis response", col="orange")


########################################
# Vapour Pressure Deficit (VPD) function
########################################

# VPD = vapour pressure deficit [kPa]
# T_min & T_max = minimum and maximum daily temperature [C]
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

#Response plot
VPD = c(0,900,4100,5000)
response = sapply(VPD, VPD.fun, VPD_min=900, VPD_max=4100)
plot(response~VPD, type="b", xlab="VPD (Pa)", ylab="Photosynthesis response", col="blue")


#####################
# Convert degC to kPa
#####################

degC_to_kPa.fun <- function(temp) {
  out       <- 0.6108*exp((17.27*temp)/(temp+237.3)) 
  return(out)
}


#######################################
# Convert specific to relative humidity
#######################################

qair2rh <- function(qair, temp, press = 1013.25){
  es <-  6.112 * exp((17.67 * temp)/(temp + 243.5))
  e <- qair * press / (0.378 * qair + 0.622)
  rh <- e / es
  rh[rh > 1] <- 1
  rh[rh < 0] <- 0
  return(rh)
}



#######################################################################################################################
#######################################################################################################################



##############################
## Merge datasets into list ##
##############################



# Identifier 1 (all site x year combinations)
PEP.df$site_year = paste0(PEP.df$pep_id,"_",PEP.df$year)

# Identifier 2 (all timeseries x year combinations)
PEP.df$ts_yr     = paste0(PEP.df$timeseries,"_",PEP.df$year)
timeseries_year  = unique(PEP.df$ts_yr)

# add PEP data (+plant functional type label) and photoperiod to list
DataList[[11]] = photo.df
DataList[[12]] = CO2.df
DataList[[13]] = PEP.df
DataList[[14]] = PrcpMonths.df

rm(photo.df, CO2.df, data, PEP.df, PrcpMonths.df)
names(DataList)=c(vn,"photoperiod",'CO2',"PEP","PrcpMonths")
names(DataList)



############################################################
## Loop through each observation using parallel computing ##
############################################################



parallelCalc <- function(timeseries_years){ 
  
  #############################################
  # Subset input data by individual observation
  #############################################
  
  #phenology data
  pheno.sub  <- DataList[[13]][which(DataList[[13]]$ts_yr==timeseries_years),]
  
  #daily mean temperature
  TMEAN <- DataList[[1]][which(DataList[[1]]$site_year==pheno.sub$site_year),]%>% 
    dplyr::select(as.character(1:366))
  
  #Skip timeseries for which there is no data
  if (nrow(TMEAN)==0) {
    write.table(pheno.sub, file=paste0(PEP_drivers_missing, '/', timeseries_years, '.csv'), sep=',', row.names = F, col.names = T)
  } else {
    
    #daily minimum temperature
    TMIN      <- DataList[[2]][which(DataList[[2]]$site_year==pheno.sub$site_year),]%>% 
      dplyr::select(as.character(1:366))
    
    #daily maximum temperature
    TMAX      <- DataList[[3]][which(DataList[[3]]$site_year==pheno.sub$site_year),]%>% 
      dplyr::select(as.character(1:366))
    
    #precipitation
    PRCP      <- DataList[[4]][which(DataList[[4]]$site_year==pheno.sub$site_year),]%>% 
      dplyr::select(as.character(1:366))
    
    #mean monthly precipitation across all years
    PRCP.site <- DataList[[14]][which(DataList[[14]]$pep_id==pheno.sub$pep_id),]
    
    #air humidity 
    QAIR      <- DataList[[5]][which(DataList[[5]]$site_year==pheno.sub$site_year),]%>% 
      dplyr::select(as.character(1:366))
    
    #soil moisture (<10cm)
    MOIST10   <- DataList[[6]][which(DataList[[6]]$site_year==pheno.sub$site_year),]%>% 
      dplyr::select(as.character(1:366))
    
    #soil moisture (10-40 cm)
    MOIST40   <- DataList[[7]][which(DataList[[7]]$site_year==pheno.sub$site_year),]%>% 
      dplyr::select(as.character(1:366))
    
    #net short-wave radiation
    SWRAD     <- DataList[[8]][which(DataList[[8]]$site_year==pheno.sub$site_year),]%>% 
      dplyr::select(as.character(1:366))
    
    #net long-wave radiation
    LWRAD     <- DataList[[9]][which(DataList[[9]]$site_year==pheno.sub$site_year),]%>% 
      dplyr::select(as.character(1:366))
    
    #short-wave radiation down
    SWRADdown <- DataList[[10]][which(DataList[[10]]$site_year==pheno.sub$site_year),]%>% 
      dplyr::select(as.character(1:366))
    
    #day length
    PHOTO     <- DataList[[11]][which(DataList[[11]]$lat==pheno.sub$lat),][1]%>% 
      dplyr::select(as.character(1:366))
    
    #CO2 (monthly)
    CO2       <- as.data.frame(t(DataList[[12]][which(DataList[[12]]$site_year==pheno.sub$site_year),]%>%
                                   dplyr::select(as.character(1:12)))) %>%
      rename(CO2 = V1)%>%
      mutate(Month = as.numeric(1:12))
    
    
    # -------------------------------------------------------------- 
    # Short Name        ; Long Name                     ; Unit [d-1]
    # -------------------------------------------------------------- 
    # SWRAD             ; Net short wave radiation flux ; W m-2
    # LWRAD             ; Net long-wave radiation flux  ; W m-2
    # SWRADdown         ; Downward short wave radiation ; W m-2
    # PRCP              ; Precipitation rate            ; mm d-1
    # MOIST10           ; Soil moisture @ 0-10cm        ; kg m-2
    # MOIST40           ; Soil moisture @ 10-40cm       ; kg m-2
    # QAIR              ; Specific humidity             ; kg kg-1
    # TMIN              ; Minimum Air Temperature       ; degC
    # TMEAN             ; Mean Air Temperature          ; degC
    # TMAX              ; Maximum Air Temperature       ; degC
    # CO2               ; CO2 concentration             ; ppm
    # -------------------------------------------------------------- 
    
    
    #################################################################################################################
    
    
    ###############################
    # Create table of daily climate
    ###############################
    
    # Generate sub-dataframe to store results
    factors.sub <- pheno.sub %>% 
      dplyr::select(pep_id,species,timeseries,year,lat,lon,alt,leaf_out,leaf_off) %>%
      mutate(CO2 = CO2$CO2[6])
    
    # Define the current year in calendar units
    year      <- as.character(pheno.sub$year)
    start_doy <- paste(year,"-01-01", sep="") 
    end_doy   <- paste(year,"-12-31", sep="")
    days      <- seq(as.Date(start_doy), as.Date(end_doy), by="days")
    
    #create table
    daily_vals <- data.frame(Year     = year,
                             Month    = 0,
                             Day      = 0,
                             Tmin     = as.numeric(TMIN), 
                             Tmean    = as.numeric(TMEAN), 
                             Tmax     = as.numeric(TMAX), 
                             SWrad    = as.numeric(SWRAD),
                             LWrad    = as.numeric(LWRAD),
                             SWradDown= as.numeric(SWRADdown),
                             Moist10  = as.numeric(MOIST10),
                             Moist40  = as.numeric(MOIST40),
                             Prcp     = as.numeric(PRCP), 
                             Qair     = as.numeric(QAIR),
                             Photo    = as.numeric(PHOTO))
    
    #Add climate variables and data wrangling
    daily_vals = daily_vals %>%
      filter(!is.na(Tmean)) %>% #delete NAs
      mutate(
        #add month and day identifiers
        Month = lubridate::month(as.Date(days,origin=days[1])),
        Day   = lubridate::day(as.Date(days,origin=days[1])),
        #Precipitation > 2mm
        PrcpDays   = as.numeric(Prcp>2),
        #relative humidity
        RH         = qair2rh(Qair, Tmean)*100,
        #dewpoint temperature
        Tdew  = weathermetrics::humidity.to.dewpoint(t = Tmean,
                                                     rh = RH,
                                                     temperature.metric = "celsius")) %>%
      #Add CO2
      left_join(CO2, by = "Month")
    
    #set NAs to 0
    daily_vals[is.na(daily_vals)] <- 0.0001
    
    
    #################################################################################################################
    
    
    #####################
    # Get important dates
    #####################
    
    
    # longest day of year (summer solstice)
    solstice = which(daily_vals$Photo==max(daily_vals$Photo))[1] 
    
    # Daylength 12h (September equinox)
    DL12 = which(daily_vals$Photo[172:nrow(daily_vals)]<=12)[1] + 171                                   
    
    # Daylength 11h
    DL11 = which(daily_vals$Photo[172:nrow(daily_vals)]<=11)[1] + 171     
    
    # leaf-out
    DOY_out <- pheno.sub$leaf_out
    
    # Senescence date
    DOY_off <- pheno.sub$leaf_off_mean
    
    # month of leaf-out
    Month_out = as.numeric(format(as.Date(pheno.sub$leaf_out_mean, origin = "1970-01-01"), "%m"))
    
    # month of leaf-off
    Month_off = as.numeric(format(as.Date(pheno.sub$leaf_off_mean, origin = "1970-01-01"), "%m"))
    
    
    #################################################################################################################
    
    
    ################################
    # Add Summer Precipitation (LPI)
    ################################
    
    
    # get three driest months
    dry.months = PRCP.site %>% 
      filter(Month>=Month_out,
             Month<=Month_off)%>%
      top_n(-3,prcp_mean) %>% 
      pull(Month)
      
    #number of days >2mm precipitation during three driest months of growing season
    factors.sub$LPI = sum(daily_vals %>% 
      filter(Month %in% dry.months) %>% 
      pull(PrcpDays))
    
    
    #################################################################################################################
    
    
    
    ################################
    ## Photosynthesis calculation ##
    ################################
    
    
    
    # GSI, Daily Net Photosynthesis rate (dA_n) and water stress factor (dw) are calculated daily 
    # and then accumulated by summation
    
    
    # Initialize vector to store daily values
    
    #growing season index
    iGSI_year         <- vector()
    
    #vapor pressure deficit
    VPD_year          <- vector()
    
    #Anet
    dA_tot_year         <- vector()
    dA_totw_year        <- vector()
    
    
    # Loop through days of the growing season
    for(i in 1:nrow(daily_vals)) {
      
      # set VPD min and max
      VPD_min <- 0.9
      VPD_max <- 4.1
      
      # Estimate photoperiod thresholds based on the maximum and minimum values of the growing season
      photo_min <- 11
      photo_max <- max(daily_vals$Photo) 
      
      # e_s: saturation vapour pressure [kPa] 
      e_s <- (degC_to_kPa.fun(temp=daily_vals$Tmax[i]) + degC_to_kPa.fun(temp=daily_vals$Tmin[i])) / 2
      
      # e_a: derived from dewpoint temperature [kPa]
      e_a <- degC_to_kPa.fun(temp=daily_vals$Tdew[i])
      
      # VPD: Vapour pressure deficit [kPa]
      VPD  <- e_s-e_a
      VPD_year <- c(VPD_year,VPD)
      
      # apply vapor pressure deficit function
      iVPD <- VPD.fun(VPD, VPD_min, VPD_max)
      
      # iOpt_temp: response to optimal temperature (Gompertz function)
      iOpt_temp  <- temp_opt.fun(daily_vals$Tmean[i])
      
      # iDL: Day-length response of photosynthetic capacity (Bauerle et al. 2012)
      if(i %in% c(1:solstice)) {iDL = 1} else {
        iDL = photoperiod.fun(daily_vals$Photo[i], photo_min, photo_max)}
      
      
      ################################
      ## Growing Season Index (GSI) ##
      ################################
      
      # modified from Jolly et al. 2005
      iGSI <- as.numeric(iVPD * iOpt_temp * iDL)
      
      # Add to the cumulative cGSI
      iGSI_year <- c(iGSI_year, iGSI)
      
      #----------------------------------------------------------------------------------------------
      
      ########################################
      ## Net day-time photosynthesis (Anet) ##
      ########################################
      
      # Net photosynthesis rate (PHOTOSYNTHESIS-CONDUCTANCE MODEL, ref. Sitch et al. 2003)
      
      # The fraction of PAR absorbed by leaves (fapar) is derived by species-specific leaf area index (LAI)
      # using the Lamber-Beer law to calculate the radiation reaching the soil (Rad_soil)
      # therefore, fapar is the remaining proportion
      # Eqn 27, Prentice et al. 1993
      fapar <- 1 - (exp(-.5 * pheno.sub$LAI))
      
      # convert in J/m^-2 day: the power in watts (W) is equal to the energy in joules (J), 
      # divided by the time period in seconds (s): 
      # --> 1 Watt = 1 Joule/second, therefore j = W*86400
      # inter-annual LAI variation simulated via LPJ-GUESS
      apar        <- fapar  * daily_vals$SWradDown[i] * cq * 60 * 60 * 24
      
      # Calculate temperature inhibition function limiting photosynthesis at low and high temperatures 
      # (ref. Sitch et al. 2002)
      tstress <- temp_opt.fun(daily_vals$Tmean[i])
      
      # Calculate catalytic capacity of rubisco, Vm, assuming optimal (non-water-stressed) value for lambda, 
      # i.e. lambdamc3
      # adjust kinetic parameters for their dependency on temperature 
      # i.e. relative change in the parameter for a 10 degC change in temperature
      # Eqn 22, Haxeltine & Prentice 1996a
      ko  <- ko25*q10ko^((daily_vals$Tmean[i]-25.0)/10.0)   # Michaelis constant of rubisco for O2
      kc  <- kc25*q10kc^((daily_vals$Tmean[i]-25.0)/10.0)   # Michaelis constant for CO2
      tau <- tau25*q10tau^((daily_vals$Tmean[i]-25.0)/10.0) # CO2/O2 specificity ratio
      
      # gammastar: CO_2 compensation point [CO2 partial pressure, Pa]   
      # Eqn 8, Haxeltine & Prentice 1996
      gammastar <- po2/(2.0*tau)
      
      # Convert ambient CO2 level from mole fraction to partial pressure, Pa
      pa <- daily_vals$CO2[i]*p
      
      # p_i: non-water-stressed intercellular CO2 partial pressure, Pa
      # Eqn 7, Haxeltine & Prentice 1996
      p_i <- pa*lambdamc3 
      
      # Calculate coefficients
      # Eqn 4, Haxeltine & Prentice 1996
      c1 <- tstress*alphac3*((p_i-gammastar)/(p_i+2.0*gammastar))
      
      # Eqn 6, Haxeltine & Prentice 1996
      c2 <- (p_i-gammastar)/(p_i+kc*(1.0+po2/ko)) 
      b  <- bc3 # choose C3 value of b for Eqn 10, Haxeltine & Prentice 1996
      t0 <- t0c3 # base temperature for temperature response of rubisco
      
      # Eqn 13, Haxeltine & Prentice 1996
      s <- (24.0 / daily_vals$Photo[i] ) * b
      
      # Eqn 12, Haxeltine & Prentice 1996
      sigma <- sqrt(max(0.0,1.0-(c2-s)/(c2-theta*s)))
      
      # vm: optimal rubisco capacity, gC m-2 d-1 
      # Eqn 11, Haxeltine & Prentice 1996
      # cmass: the atomic weight of carbon, used in unit conversion from molC to g 
      # cq: conversion factor from apar [J m-2] to photosynthetic photon flux density [mol m-2]
      vm        <- (1.0/b)*(c1/c2)*((2.0*theta-1.0)*s-(2.0*theta*s-c2)*sigma)*apar*cmass
      
      # je: PAR-limited photosynthesis rate, gC m-2 h-1
      # Eqn 3, Haxeltine & Prentice 1996
      # Convert je from daytime to hourly basis
      if(daily_vals$Photo[i]==0) {
        je        <- 0
      } else {
        je        <- c1*apar * cmass / daily_vals$Photo[i]
      }
      
      # jc: rubisco-activity-limited photosynthesis rate, gC m-2 h-1
      # Eqn 5, Haxeltine & Prentice 1996
      jc        <- c2*vm/24.0
      
      # agd: daily gross photosynthesis, gC m-2 d-1
      # Eqn 2, modified with k_shape (theta)
      if(je<1e-10 | jc<=1e-10) {
        agd        <- 0
      } else {
        agd <- (je + jc - sqrt((je + jc)^2.0 - 4.0 * theta * je * jc)) / (2.0 * theta) * daily_vals$Photo[i] * iDL
      }
      
      # rd: daily leaf respiration, gC m-2 d-1
      # Eqn 10, Haxeltine & Prentice 1996
      rd        <- b*vm
      
      # and: daily net photosynthesis (at leaf level), gC m-2 d-1
      and        <- agd-rd
      
      # adt: total daytime net photosynthesis, gC m-2 d-1
      # Eqn 19, Haxeltine & Prentice 1996
      adt        <- and        + (1.0 - daily_vals$Photo[i] / 24.0) * rd
      
      # Store the daily result in the yearly vector
      dA_tot_year        <- c(dA_tot_year       , adt)  
      
      
      ## Water Stress Factor (ref. Gerten et al. 2004)
      ################################################
      
      # soil is treated as a simple bucket consisting of two layers with fixed thickness
      
      # Calculate potential evapotranspiration (ETA) rate, E_pot, mm d-1
      
      # delta: rate of increase of the saturation vapour pressure with temperature
      delta <- (2.503*10^6 * exp((17.269 * daily_vals$Tmean[i]) / (237.3 + daily_vals$Tmean[i]))) / (237.3 + daily_vals$Tmean[i])^2
      
      # R_n: instantaneous net radiation, W m-2 = net short-wave radiation flux + net long-wave flux 
      R_n <- daily_vals$SWrad[i] + daily_vals$LWrad[i] 
      
      # E_eq: equilibrium EvapoTranspiration
      # from seconds to day
      E_eq <- 24 * 3600 * (delta / (delta + gamma)) * (R_n / L) 
      
      # E_pot: potential EvapoTranspiration = equilibrium ETA * Priestley-Taylor coefficient 
      E_pot <- E_eq*a_m
      
      # ratio: stomata-controlled ratio between intercellular and ambient CO2 partial pressure in the absence of water limitation
      ratio <- p_i/pa # ca. 0.8
      
      # g_min: minimum canopy conductance, mm s-1
      # depends on PFT
      if(pheno.sub$PFT=="TBL") { 
        g_min <- 0.5*3600*24 # from seconds to day
      } else {
        g_min <- 0.3*3600*24
      }
      
      # g_pot: non-water-stressed potential canopy conductance, mm s-1
      g_pot <- g_min + ((1.6*adt)/((pa/p)*(1-ratio)))
      
      # E_demand: atmospheric demand 
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
      
      # relative soil moisture wr:
      # ratio between current soil water content and plant-available water capacity
      # wr ratio is computed for both soil layers by 
      # weighting their relative soil water contents (w1, w2) 
      # with the fraction of roots present in the respective layer
      w1  <- daily_vals$Moist10[i]
      w2  <- daily_vals$Moist40[i]
      
      # soil texture-dependent difference between field capacity and wilting point w_max [%]
      w_max <- pheno.sub$w_max
      wr <- root1*(w1/w_max) + root2*(w2/w_max)
      
      # E_supply: plant- and soil-limited supply function 
      E_supply <- as.numeric(E_max*wr)
      
      # dw: daily water stress factor
      dw <- min(1,(E_supply/E_demand))
      dw[dw<0]=0
      
      # dA_totw: daily net photosynthesis modified by water stress factor
      dA_totw        <- adt*dw
      
      # Add daily result to the yearly vector
      dA_totw_year     <- c(dA_totw_year, dA_totw)
      
    } # END loop through days of the growing season
    
    #set negative values to zero
    VPD_year[VPD_year<=0] = 0.001
    
    #----------------------------------------------------------------------------------------------
    
    #Store the results
    ##################
    
    daily_vals = daily_vals %>%
      mutate(VPD           = VPD_year *1000,      #VPD (in Pa)
             GSI           = iGSI_year,           #Growing-season index (GSI)
             AnetNW        = dA_tot_year,         #net daytime photosynthesis (no water function)
             Anet          = dA_totw_year         #net daytime photosynthesis (water-stressed)
      )  %>%
      rename(Tls=Tmean) 
    
    
    ###############################################################################################################
    
    
    ###################
    ## Store drivers ##
    ###################
    
    
    ############################
    ## Growing season drivers ##
    ############################
    
    #define variables
    variable.names = c('Anet', 'AnetNW', 'GSI',
                       'Tls', 'VPD')
    
    #---------------------------------------------------------------------------------------------------------
    
    for(i in 1:length(variable.names)) {
      
      #choose variable (daily values)
      variable = daily_vals[,variable.names[i]]
      
      #---------------------------------------------------------------------------------------------------------
      
      # Name variables
      ################
      
      # Off...date of timeseries average senescence
      # DL12/11...date at which day length falls below 12h/11h
      varname.Off     <- paste(variable.names[i], "Off",sep=".")
      varname.DL12    <- paste(variable.names[i], "DL12", sep=".")
      varname.DL11    <- paste(variable.names[i], "DL11", sep=".")
      
      #---------------------------------------------------------------------------------------------------------
      
      # Create columns
      ################
      
      
      #####################################
      # Cumulative photosynthesis and VPD #
      #####################################
      
      if(variable.names[i] %in% c('Anet','AnetNW','GSI','VPD')){
        
        factors.sub = factors.sub %>%
          mutate(!!varname.Off     := sum(variable[DOY_out:DOY_off]),
                 !!varname.DL12    := sum(variable[DOY_out:DL12]),
                 !!varname.DL11    := sum(variable[DOY_out:DL11]) )
      } 
      
      
      ###############################################################################
      # Leafy season mean temperature [Tls] (Liu et al. 2018, Ecological Modelling) #
      ###############################################################################
      
      if(variable.names[i] %in% c("Tls")){
        
        factors.sub = factors.sub %>%
          mutate(!!varname.Off     := mean(variable[DOY_out:DOY_off]),
                 !!varname.DL12    := mean(variable[DOY_out:DL12]),
                 !!varname.DL11    := mean(variable[DOY_out:DL11]) )
      }
    }
    
    
    #--------------------------------------------------------------------------
    
    
    ########################
    ## Autumn temperature ##
    ########################
    
    ## Calculate the average preseason temperature 60 days prior to mean senescence date
    factors.sub = factors.sub %>%
      mutate(Tautumn = mean(daily_vals$Tmin[(DOY_off-60):DOY_off]))
    
    
    #################################################################################################################
    
    
    # Safe the table
    write.table(factors.sub, file=paste0(PEP_drivers_files, '/', timeseries_years, '.csv'), sep=',', row.names =F, col.names = T)
    
  }
}



#######################################################################################################################
#######################################################################################################################



##################
## Run the Loop ##
##################



#initialize the loop
outputlist <- pbmclapply(timeseries_year[1:50], parallelCalc, mc.cores=12, mc.preschedule=T)

#check how many files there are
length(list.files(path=PEP_drivers_files, pattern='.csv'))
length(list.files(path=PEP_drivers_final, pattern='.csv'))
length(list.files(path=PEP_drivers_missing, pattern='.csv'))

#Rbind files
climate.factors.table = rbindlist(lapply(list.files(path = PEP_drivers_files), 
                                         function(n) fread(file.path(PEP_drivers_files, n)))) %>%
  #delete timeseries with fewer than 15 years
  group_by(timeseries) %>%
  filter(n() >= 15) %>% 
  ungroup()



#######################################################################################################################
#######################################################################################################################



###################
## Safe the data ##
###################



#Safe table
write.csv(climate.factors.table, paste(PEP_drivers_final, "DataMeta_3_Drivers_subset.csv", sep="/"))

#Remove individual files
do.call(file.remove, list(list.files(PEP_drivers_files, 
                                     full.names = TRUE)))





#####################
## Reproducibility ##	
#####################



## date time
Sys.time()

## session info
sessionInfo()



#######################################################################################################################
############################################################# THE END #################################################
#######################################################################################################################