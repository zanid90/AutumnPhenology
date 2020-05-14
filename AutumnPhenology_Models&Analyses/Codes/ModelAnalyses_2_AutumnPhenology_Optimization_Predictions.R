###################################
## Model Optimization & Predictions

# From the text:
# See Materials and Methods
# Autumn phenology modelling

# Define directory path
setwd(".../AutumnPhenology/AutumnPhenology_Models&Analyses/Data/")

# Load libraries
library(data.table)
library(tidyverse)
library(phenor)
library(car)
library(Metrics)
library(broom)
library(dismo)


##----------------------------------------
## MODELS OF LEAF SENESCENCE (and drivers):

## First-generation:
# CDD (chilling temperature) - Dufrene et al. (2005)
# DM1 and DM2 (chilling temperature, autumn daylength) - Delpierre et al. (2009)
# TPM (chilling temperature, autumn daylength) - Lang et al. (2019)
## Second-generation:
# SIAM (chilling temperature, autumn daylength, spring anomaly) - Keenan and Richardson (2015)
# TDM and TPDM (chilling temperature, autumn daylength, growing season temperature / + water stress) - Liu et al. (2019)
## PIA:
# PIA_gsi (chilling temperature, autumn daylength, leaf flushing date, growing season mean temperature, daylength, vapour pressure deficit)
# PIA-/+ (chilling temperature, autumn daylength, leaf flushing date, growing season mean temperature, daylength, precipitation, net radiation, CO2 concentration, -/+ water stress)

# Define functions of Autumn phenology Models
# Modified from https://github.com/khufkens/phenor/blob/master/R/phenology_models.R
CDD.model = function(par, data){
  # exit the routine as some parameters are missing
  if (length(par) != 2){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  T_base = par[1]
  F_crit = par[2]
  
  # create forcing/chilling rate vector
  Rf = data$Tmini - T_base
  Rf[Rf > 0] = 0
  
  # photoperiod-dependent start-date for chilling accumulation (t0)
  # t0 is defined as the first day when daily minimum temperature is lower than a temperature threshold (T_base) 
  # after the date of the peak multiyear average daily minimum temperature, namely the 200th day of year
  t0 <- vector()
  for(c in 1:ncol(data$Tmini)) {
    interval = 1:366
    t0A = interval[which(data$Tmini[,c] < T_base)]
    ind1 = min(which(t0A > 200))
    t0A = t0A[ind1]
    t0 = c(t0,t0A)
  }
  # nullify values before the t0
  for(c in 1:ncol(data$Tmini)){
    Rf[1:t0[c],c] = 0 
  }
  
  # predict date of leaf.off according to optimized F_crit
  doy = apply(Rf,2, function(xt){
    doy = which(cumsum(xt) <= F_crit)[1]
  })
  
  return(doy)
}
DM1.model = function(par, data){
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  T_base = par[1]
  P_base = par[2]
  F_crit = par[3]
  
  # create forcing/chilling rate vector at the day level
  Rf = (data$Tmini - T_base)*(data$Li/P_base) # lengthening photoperiod promoting leaf senescence
  Rf[Rf > 0] = 0
  
  # photoperiod-dependent start-date for chilling accumulation (t0)
  # t0 is defined as the first day when photoperiod is shorter than the photoperiod threshold (P_base)
  # after the date of the longest photoperiod (summer solstice), namely, the 173rd day of year
  t0 <- vector()
  for(c in 1:ncol(data$Tmini)) {
    interval = 1:366
    t0A = interval[which(data$Li[,c] < P_base)]
    ind1 = min(which(t0A > 173))
    t0A = t0A[ind1]
    t0 = c(t0,t0A)
  }
  
  # nullify values before the t0
  for(c in 1:ncol(data$Tmini)){
    Rf[1:t0[c],c] = 0 #nullify values before the date of leaf.out
  }
  
  # calculate the summation along the year (interval = 1:366) and derive the date of leaf.off
  # DOY of budburst criterium
  doy = apply(Rf,2, function(xt){
    doy = which(cumsum(xt) <= F_crit)[1]
  })
  
  return(doy)
}
DM2.model = function(par, data){
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  T_base = par[1]
  P_base = par[2]
  F_crit = par[3]
  
  # create forcing/chilling rate vector at the day level
  Rf = (data$Tmini - T_base)*(1-(data$Li/P_base)) # shortening photoperiod promoting leaf senescence
  Rf[Rf > 0] = 0
  
  # photoperiod-dependent start-date for chilling accumulation (t0)
  # t0 is defined as the first day when photoperiod is shorter than the photoperiod threshold (P_base)
  # after the date of the longest photoperiod (summer solstice), namely, the 173rd day of year
  t0 <- vector()
  for(c in 1:ncol(data$Tmini)) {
    interval = 1:366
    t0A = interval[which(data$Li[,c] < P_base)]
    ind1 = min(which(t0A > 173))
    t0A = t0A[ind1]
    t0 = c(t0,t0A)
  }
  
  # nullify values before the t0
  for(c in 1:ncol(data$Tmini)){
    Rf[1:t0[c],c] = 0 #nullify values before the date of leaf.out
  }
  
  # calculate the summation along the year (interval = 1:366) and derive the date of leaf.off
  # DOY of budburst criterium
  doy = apply(Rf,2, function(xt){
    doy = which(cumsum(xt) <= F_crit)[1]
  })
  
  return(doy)
}
TPM.model = function(par, data){
  # exit the routine as some parameters are missing
  if (length(par) != 4){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the par argument for readability
  P_base = par[1]
  a = par[2]
  b = par[3]
  F_crit = par[4]
  
  # create forcing/chilling rate vector at the day level
  Rf = 1/(1+exp(a*(data$Tmini*data$Li-b)))
  
  # photoperiod-dependent start-date for chilling accumulation (t0)
  # t0 is defined as the first day when photoperiod is shorter than the photoperiod threshold (P_base)
  # after the date of the longest photoperiod (summer solstice), namely, the 173rd day of year
  t0 <- vector()
  for(c in 1:ncol(data$Tmini)) {
    interval = 1:366
    t0A = interval[which(data$Li[,c] < P_base)]
    ind1 = min(which(t0A > 173))
    t0A = t0A[ind1]
    t0 = c(t0,t0A)
  }
  
  # nullify values before the t0
  for(c in 1:ncol(data$Tmini)){
    Rf[1:t0[c],c] = 0 #nullify values before the date of leaf.out
  }
  
  # calculate the summation along the year and derive the date of leaf.off
  # DOY of budburst criterium
  doy = apply(Rf,2, function(xt){
    doy = which(cumsum(xt) >= F_crit)[1]
  })
  
  return(doy)
}
SecondGen_PIA.models = function(par, predictor, data) {
  # exit the routine as some parameters are missing
  if (length(par) != 5 & length(par) != 6){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  P_base = as.numeric(par[1])
  a = as.numeric(par[2])
  b = as.numeric(par[3])
  c = as.numeric(par[4])
  if(length(par)==5) {
    d = as.numeric(par[5])
    pred = predictor
  }
  if(length(par)==6) {
    d = as.numeric(par[5])
    e = as.numeric(par[6])
    pred1 = predictor[1]
    pred2 = predictor[2]
  }
  
  # create forcing/chilling rate vector at the day level
  Rf = 1/(1+exp(a*(data$Tmini*data$Li-b)))
  
  # photoperiod-dependent start-date for chilling accumulation (t0)
  # t0 is defined as the first day when photoperiod is shorter than the photoperiod threshold (P_base)
  # after the date of the longest photoperiod (summer solstice), namely, the 173rd day of year
  t0 <- vector()
  for(col in 1:ncol(data$Tmini)) {
    interval = 1:366
    t0A = interval[which(data$Li[,col] < 173)]
    ind1 = min(which(t0A > 173))
    t0A = t0A[ind1]
    t0 = c(t0,t0A)
  }
  
  # nullify values before the t0
  for(col in 1:ncol(data$Tmini)){
    Rf[1:t0[col],col] = 0 
  }
  
  if(length(par)==5) {
    # add predictor at the end of the matrix-columns
    Rf = rbind(Rf,predictor)
    
    # predict date of leaf.off
    doy = apply(Rf,2, function(xt){
      doy = which(cumsum(xt[1:366]) >= c+d*xt[367])[1]
    }) 
  }
  if(length(par)==6) {
    # add predictors at the end of the matrix-columns
    Rf = rbind(Rf,pred1)
    Rf = rbind(Rf,pred2)
    
    # predict date of leaf.off
    doy = apply(Rf,2, function(xt){
      doy = which(cumsum(xt[1:366]) >= c+d*xt[367]+e*xt[368])[1]
    }) 
  }
  
  return(doy)
}


##----------------------------------------
## Import & Format data

# Phenological observations
# DoY_off: leaf senescence transition dates
# DoY_out: leaf flushing transition dates
all.df <- fread("DataMeta_3_Drivers.csv") 
pheno.df <- all.df %>% 
  select("timeseries","PEP_ID","LON","LAT","Species","YEAR","DoY_off","DoY_out")

# Predictors
# DoY_out: leaf flushing transition dates
# temp_GS: growing season temperature
# RD_summer: number of rainy days (precipitation >2mm) during the driest months
# cGSI: cumulative Growing Season Index
# cA_tot: cumulative net photosynthesis modified by water-stress factor
# cA_tot-w: cumulative net photosynthesis without water-stress factor
preds.df <- all.df %>% 
  select("timeseries","PEP_ID","LON","LAT","Species","YEAR","DoY_out","temp_GS","RD_summer","cGSI","cA_tot","cA_tot-w")

# Minimum temperature and photoperiod
tmin.df <- fread("Minimum Temperature.csv")
photo.df <- fread("Photoperiod.csv")
tmin.df$ts_yr <- paste0(tmin.df$PEP_ID,"_",tmin.df$YEAR)
photo.df$ts_yr <- paste0(photo.df$PEP_ID,"_",photo.df$YEAR)

# Define model names
models   <- c("CDD","DM1","DM2","TPM","SIAM","TDM","TPDM","PIA_gsi","PIA-","PIA+")

# Define species
species <- unique(pheno.df$Species)

# Define timeseries
timeseries <- unique(pheno.df$timeseries)

# Prepare input datasets for PHENOR package
tmin.df <- tmin.df %>% 
  select(-c(YEAR,PEP_ID,LAT,LON))
photo.df <- photo.df %>% 
  select(-c(YEAR,PEP_ID,LAT))
phenor_input <- merge(tmin.df,photo.df,by="ts_yr")
phenor_input <- merge(pheno.df,input,by="ts_yr")
rm(all.df,pheno.df,tmin.df,photo.df)

# Create a DataList to store all subsets for each species
DataList <- replicate(length(species), data.frame())
names(DataList) <- paste0("DataList","_",species)
SiteList <- replicate(length(species), data.frame())

# Initiate an external loop to subset for each species
for(sp in 1:length(species)) {
  
  # Subset phenological dataset for each species
  pheno_sp.sub <- phenor_input[which(phenor_input$Species==species[sp]),]
  pheno_sp.sub <- as.data.frame(pheno_sp.sub)
  
  # Find the sites per species
  sites <- unique(pheno_sp.sub$PEP_ID)
  SiteList[[sp]] <- sites
  
  # Add empty subsets
  DataList[[sp]] <- replicate(length(sites),data.frame())
  names(DataList[[sp]]) <- paste0("data_",sites,"_",pheno_sp.sub[1,]$Species)
  
  # Counter for sites
  count  <- 0 
  
  # Loop for each site
  for(site in sites){
    index <- which(pheno_sp.sub$PEP_ID==site)
    data = list("doy" = as.vector(pheno_sp.sub[index,]$DoY_out),
                "site" = as.vector(paste0(pheno_sp.sub[index,]$PEP_ID,"_",pheno_sp.sub[1,]$Species)),
                "location" = t(pheno_sp.sub[index,c("LON","LAT")]), 
                "year" = as.vector(pheno_sp.sub[index,]$YEAR),
                "Ti" = NULL,
                "Tmini" = t(pheno_sp.sub[index,paste0(1:366,".x")]), 
                "Tmaxi" = NULL,
                "Li" = t(pheno_sp.sub[index,paste0(1:366,".y")]),
                "SPEI" = NULL,
                "VPDi" = NULL,
                "transition_dates" = as.vector(pheno_sp.sub[index,]$DoY_off),
                "georeferencing" = NULL
    )
    
    # Store each site-specific dataframe in the DataList of the corresponding species
    count <- count+1
    DataList[[sp]][[count]] <- data
  }
  print(paste0(species[sp]," DONE!"))
}
rm(index,data,pheno_sp.sub,sites,site,count,sp)


##----------------------------------------
## Predictions

# Initialize dataframe
DoYoff_Preds.df <- data.frame() 
opt_pars.df <- data.frame() 

for(sp in 1:length(species)) {
  
  # Initialize datasets per species
  DoYoff_Preds.sp <- data.frame()
  opt_pars.sp <- data.frame()
  
  # Calculate number of sites per species
  sites <- SiteList[[sp]]
  
  for(site in 1:length(sites)) {
    
    # Subset according to timeseries
    ts <- paste0(species[sp],"_",sites[site])
    data.sub <- DataList[[sp]][[site]]
    preds.sub <- preds.df %>% 
      filter(timeseries==ts)
    
    # Initialize sub-dataframes to store results
    DoYoff_Preds.sub <- data.frame(timeseries=data.sub$site, Species=species[sp], YEAR=dat.sub$year)
    DoYoff_Preds.sub$Obs_DoYoff <- data.sub$transition_dates
    opt_pars.sub <- data.frame(timeseries=data.sub$site, Species=species[sp], PEP_ID=sites)
    
    ## Parameter optimization and Prediction of leaf senescence dates
    # PHENOR package (Hufkenset al., 2018)
    
    ## CDD model
    optimal_pars <- optimize_parameters(par = NULL,
                                        data = data.sub,
                                        cost = rmse,
                                        model = "CDD.model",
                                        method = "GenSA",
                                        lower = c(15,-3000),
                                        upper = c(30,0),
                                        control = list(max.call = 40000))
    opt_pars.sub$Tbase_CDD <- optimal_pars[1]
    opt_pars.sub$Fcrit_CDD <- optimal_pars[2]
    DoYoff_Preds.sub$Pred_DoYoff_CDD <- estimate_phenology(par = optimal_pars$par,
                                                           data = data.sub,
                                                           model = "CDD.model")
    
    ## DM1 model
    optimal_pars <- optimize_parameters(par = NULL,
                                        data = data.sub,
                                        cost = rmse,
                                        model = "DM1.model",
                                        method = "GenSA",
                                        lower = c(15,11,-2000),
                                        upper = c(30,16,0),
                                        control = list(max.call = 40000))
    opt_pars.sub$Tbase_DM1 <- optimal_pars[1]
    opt_pars.sub$Pbase_DM1 <- optimal_pars[2]
    opt_pars.sub$Fcrit_DM1 <- optimal_pars[3]
    DoYoff_Preds.sub$Pred_DoYoff_DM1 <- estimate_phenology(par = optimal_pars$par,
                                                           data = data.sub,
                                                           model = "DM1.model")
    
    ## DM2 model
    optimal_pars <- optimize_parameters(par = NULL,
                                        data = data.sub,
                                        cost = rmse,
                                        model = "DM2.model",
                                        method = "GenSA",
                                        lower = c(15,11,-2000),
                                        upper = c(30,16,0),
                                        control = list(max.call = 40000))
    opt_pars.sub$Tbase_DM2 <- optimal_pars[1]
    opt_pars.sub$Pbase_DM2 <- optimal_pars[2]
    opt_pars.sub$Fcrit_DM2 <- optimal_pars[3]
    DoYoff_Preds.sub$Pred_DoYoff_DM2 <- estimate_phenology(par = optimal_pars$par,
                                                           data = data.sub,
                                                           model = "DM2.model")
    
    ## TPM model
    optimal_pars <- optimize_parameters(par = NULL,
                                        data = data.sub,
                                        cost = rmse,
                                        model = "TPM.model",
                                        method = "GenSA",
                                        lower = c(11,0.02,100,0),
                                        upper = c(16,0.1,250,200),
                                        control = list(max.call = 40000))
    opt_pars.sub$Pbase_TPM <- optimal_pars[1]
    opt_pars.sub$a_TPM <- optimal_pars[2]
    opt_pars.sub$b_TPM <- optimal_pars[3]
    opt_pars.sub$Fcrit_TPM <- optimal_pars[4]
    DoYoff_Preds.sub$Pred_DoYoff_TPM <- estimate_phenology(par = optimal_pars$par,
                                                           data = data.sub,
                                                           model = "TPM.model")
    
    ## SIAM model
    optimal_pars <- optimize_parameters(par = NULL,
                                        predictor = preds.sub$DoY_out,
                                        data = data.sub,
                                        cost = rmse,
                                        model = "SecondGen_PIA.models",
                                        method = "GenSA",
                                        lower = c(11,0.02,100,0,0),
                                        upper = c(16,0.1,250,300,1),
                                        control = list(max.call = 40000))
    opt_pars.sub$Pbase_SIAM <- optimal_pars[1]
    opt_pars.sub$a_SIAM <- optimal_pars[2]
    opt_pars.sub$b_SIAM <- optimal_pars[3]
    opt_pars.sub$c_SIAM <- optimal_pars[4]
    opt_pars.sub$d_SIAM <- optimal_pars[5]
    DoYoff_Preds.sub$Pred_DoYoff_SIAM <- estimate_phenology(par = optimal_pars$par,
                                                            predictor = preds.sub$DoY_out,
                                                            data = data.sub,
                                                            model = "SecondGen_PIA.models")
    
    ## TDM model
    optimal_pars <- optimize_parameters(par = NULL,
                                        predictor = preds.sub$temp_GS,
                                        data = data.sub,
                                        cost = rmse,
                                        model = "SecondGen_PIA.models",
                                        method = "GenSA",
                                        lower = c(11,0.02,100,0,0),
                                        upper = c(16,0.1,250,300,1),
                                        control = list(max.call = 40000))
    opt_pars.sub$Pbase_TDM <- optimal_pars[1]
    opt_pars.sub$a_TDM <- optimal_pars[2]
    opt_pars.sub$b_TDM <- optimal_pars[3]
    opt_pars.sub$c_TDM <- optimal_pars[4]
    opt_pars.sub$d_TDM <- optimal_pars[5]
    DoYoff_Preds.sub$Pred_DoYoff_TDM <- estimate_phenology(par = optimal_pars$par,
                                                           predictor = preds.sub$temp_GS,
                                                           data = data.sub,
                                                           model = "SecondGen_PIA.models")
    
    ## TPDM model
    optimal_pars <- optimize_parameters(par = NULL,
                                        predictor = c(preds.sub$temp_GS,preds.sub$RD_summer),
                                        data = data.sub,
                                        cost = rmse,
                                        model = "SecondGen_PIA.models",
                                        method = "GenSA",
                                        lower = c(11,0.02,100,0,0,0),
                                        upper = c(16,0.1,250,300,1,1),
                                        control = list(max.call = 40000))
    opt_pars.sub$Pbase_TPDM <- optimal_pars[1]
    opt_pars.sub$a_TPDM <- optimal_pars[2]
    opt_pars.sub$b_TPDM <- optimal_pars[3]
    opt_pars.sub$c_TDPM <- optimal_pars[4]
    opt_pars.sub$d_TPDM <- optimal_pars[5]
    opt_pars.sub$e_TPDM <- optimal_pars[6]
    DoYoff_Preds.sub$Pred_DoYoff_TPDM <- estimate_phenology(par = optimal_pars$par,
                                                            predictor = c(preds.sub$temp_GS,preds.sub$RD_summer),
                                                            data = data.sub,
                                                            model = "SecondGen_PIA.models")
    
    ## PIA_gsi model
    optimal_pars <- optimize_parameters(par = NULL,
                                        predictor = preds.sub$cGSI,
                                        data = data.sub,
                                        cost = rmse,
                                        model = "SecondGen_PIA.models",
                                        method = "GenSA",
                                        lower = c(11,0.02,100,0,0),
                                        upper = c(16,0.1,250,300,1),
                                        control = list(max.call = 40000))
    opt_pars.sub$Pbase_PIAgsi <- optimal_pars[1]
    opt_pars.sub$a_PIAgsi <- optimal_pars[2]
    opt_pars.sub$b_PIAgsi <- optimal_pars[3]
    opt_pars.sub$c_PIAgsi <- optimal_pars[4]
    opt_pars.sub$d_PIAgsi <- optimal_pars[5]
    DoYoff_Preds.sub$Pred_DoYoff_PIAgsi <- estimate_phenology(par = optimal_pars$par,
                                                              predictor = preds.sub$cGSI,
                                                              data = data.sub,
                                                              model = "SecondGen_PIA.models")
    
    ## PIA- model
    optimal_pars <- optimize_parameters(par = NULL,
                                        predictor = preds.sub$`cA_tot-w`,
                                        data = data.sub,
                                        cost = rmse,
                                        model = "SecondGen_PIA.models",
                                        method = "GenSA",
                                        lower = c(11,0.02,100,0,0),
                                        upper = c(16,0.1,250,300,1),
                                        control = list(max.call = 40000))
    opt_pars.sub$`Pbase_PIA-` <- optimal_pars[1]
    opt_pars.sub$`a_PIA-` <- optimal_pars[2]
    opt_pars.sub$`b_PIA-` <- optimal_pars[3]
    opt_pars.sub$`c_PIA-` <- optimal_pars[4]
    opt_pars.sub$`d_PIA-` <- optimal_pars[5]
    DoYoff_Preds.sub$`Pred_DoYoff_PIA-` <- estimate_phenology(par = optimal_pars$par,
                                                              predictor = preds.sub$`cA_tot`,
                                                              data = data.sub,
                                                              model = "SecondGen_PIA.models")
    
    ## PIA+ model
    optimal_pars <- optimize_parameters(par = NULL,
                                        predictor = preds.sub$cA_tot,
                                        data = data.sub,
                                        cost = rmse,
                                        model = "SecondGen_PIA.models",
                                        method = "GenSA",
                                        lower = c(11,0.02,100,0,0),
                                        upper = c(16,0.1,250,300,1),
                                        control = list(max.call = 40000))
    opt_pars.sub$`Pbase_PIA+` <- optimal_pars[1]
    opt_pars.sub$`a_PIA+` <- optimal_pars[2]
    opt_pars.sub$`b_PIA+` <- optimal_pars[3]
    opt_pars.sub$`c_PIA+` <- optimal_pars[4]
    opt_pars.sub$`d_PIA+` <- optimal_pars[5]
    DoYoff_Preds.sub$`Pred_DoYoff_PIA+` <- estimate_phenology(par = optimal_pars$par,
                                                              predictor = preds.sub$cA_tot,
                                                              data = data.sub,
                                                              model = "SecondGen_PIA.models")
    
    # Autumn anomalies
    DoYoff_Preds.sub$meansite_DoYoff <- mean(data.sub$transition_dates)
    DoYoff_Preds.sub$Obs_AnomDoYoff <- data.sub$transition_dates - DoYoff_Preds.sub$meansite_DoYoff
    DoYoff_Preds.sub$Preds_AnomDoYoff_CDD <- DoYoff_Preds.sub$Pred_DoYoff_CDD - DoYoff_Preds.sub$meansite_DoYoff
    DoYoff_Preds.sub$Preds_AnomDoYoff_DM1 <- DoYoff_Preds.sub$Pred_DoYoff_DM1 - DoYoff_Preds.sub$meansite_DoYoff
    DoYoff_Preds.sub$Preds_AnomDoYoff_DM2 <- DoYoff_Preds.sub$Pred_DoYoff_DM2 - DoYoff_Preds.sub$meansite_DoYoff
    DoYoff_Preds.sub$Preds_AnomDoYoff_TPM <- DoYoff_Preds.sub$Pred_DoYoff_TPM - DoYoff_Preds.sub$meansite_DoYoff
    DoYoff_Preds.sub$Preds_AnomDoYoff_SIAM <- DoYoff_Preds.sub$Pred_DoYoff_SIAM - DoYoff_Preds.sub$meansite_DoYoff
    DoYoff_Preds.sub$Preds_AnomDoYoff_TDM <- DoYoff_Preds.sub$Pred_DoYoff_TDM - DoYoff_Preds.sub$meansite_DoYoff
    DoYoff_Preds.sub$Preds_AnomDoYoff_TPDM <- DoYoff_Preds.sub$Pred_DoYoff_TPDM - DoYoff_Preds.sub$meansite_DoYoff
    DoYoff_Preds.sub$Preds_AnomDoYoff_PIAgsi <- DoYoff_Preds.sub$Pred_DoYoff_PIAgsi - DoYoff_Preds.sub$meansite_DoYoff
    DoYoff_Preds.sub$`Preds_AnomDoYoff_PIA-` <- DoYoff_Preds.sub$`Pred_DoYoff_PIA-` - DoYoff_Preds.sub$meansite_DoYoff
    DoYoff_Preds.sub$`Preds_AnomDoYoff_PIA+` <- DoYoff_Preds.sub$`Pred_DoYoff_PIA+` - DoYoff_Preds.sub$meansite_DoYoff
    
    # Bind site-datasets
    DoYoff_Preds.sp <- rbind(DoYoff_Preds.sp,DoYoff_Preds.sub)
    opt_pars.sp <- rbind(opt_pars.sp,opt_pars.sub)
    print(paste0("RUNNING: ",DoYoff_Preds.sub$timeseries," ",site," OF ",length(sites)))
  }
  # Bind final datasets
  DoYoff_Preds.df <- rbind(DoYoff_Preds.df,DoYoff_Preds.sp)
  opt_pars.df <- rbind(opt_pars.df,opt_pars.sp)
}
write.table(DoYoff_Preds.df,"ModelAnalysis_1_Predicted_DoYoff.csv",sep=";",row.names = F)
write.table(opt_pars.df,"ModelAnalysis_2_OptimalParameters.csv",sep=";",row.names = F)


##----------------------------------------
## References

# Dufrêne, E. et al. Modelling carbon and water cycles in a beech forest: Part I: Model description and uncertainty analysis on modelled NEE. Ecol. Modell. 185, 407-436 (2005).
# Delpierre, N. et al. Modelling interannual and spatial variability of leaf senescence for three deciduous tree species in France. Agric. For. Meteorol. 149, 938-948 (2009).
# Keenan, T. F. & Richardson, A. D. The timing of autumn senescence is affected by the timing of spring phenology: Implications 434 for predictive models. Glob. Chang. Biol. 21, 2634-2641 (2015).
# Lang, W., Chen, X., Qian, S., Liu, G. & Piao, S. A new process-based model for predicting autumn phenology: How is leaf senescence controlled by photoperiod and temperature coupling? Agric. For. Meteorol. 268, 124-135 (2019).
# Liu, G., Chen, X., Fu, Y. & Delpierre, N. Modelling leaf coloration dates over temperate China by considering effects of leafy season climate. 460 Ecol. Modell. 394, 34-43 (2019).
# Hufkens, K., Basler, D., Milliman, T., Melaas, E. K. & Richardson, A. D. An integrated phenology modelling framework in r. Methods Ecol. Evol. 9, 1276-1285 (2018).