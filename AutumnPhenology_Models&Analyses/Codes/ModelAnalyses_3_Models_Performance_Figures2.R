####################
## Model Performance

# From the text:
# We evaluated the quality of the model predictions with three performance indexes 
# of observed (y-axis) vs. predicted (x-axis) values of autumn anomaly: 
# 1) the overall fit (R2 values);
# 2) the error fit (root mean square error, RMSE); 
# 3) the accuracy of the fit (slope values). 
# Additionally, the site-specific performance of the models was evaluated using 
# 5-fold cross-validation over the observation period at the time-series scale ("site-specific validation").

# Define directory path
setwd(".../AutumnPhenology/AutumnPhenology_Models&Analyses/Data/")

# Load libraries
library(data.table)
library(tidyverse)
library(broom)
library(Metrics)
library(lmodel2)
library(ggplot2)
library(phenor)

# Define auxiliary functions
se <- function(x) sqrt(var(x)/length(x))


##----------------------------------------
## Calculation of Performance Indexes

# Import data
pred_DoYoff <- fread("ModelAnalysis_1_Predicted_DoYoff.csv")
pred_DoYoff <- pred_DoYoff %>% 
  filter(Species != "Betula pubescens")

# Define model names
model.names   <- c("CDD","DM1","DM2","TPM",
                   "SIAM","TDM","TPDM",
                   "PIA_gsi","PIA-","PIA+")

# Format data for calculation 
pred_DoYoff <- pred_DoYoff %>% 
  select(timeseries,Obs_DoYoff,
         Pred_DoYoff_CDD,Pred_DoYoff_DM1,Pred_DoYoff_DM2,Pred_DoYoff_TPM,
         Pred_DoYoff_SIAM,Pred_DoYoff_TDM,Pred_DoYoff_TPDM,
         Pred_DoYoff_PIAgsi,`Pred_DoYoff_PIA-`,`Pred_DoYoff_PIA+`)
colnames(pred_DoYoff) <- c("timeseries","observations",model.names)
pred_DoYoff <- pred_DoYoff %>% 
  pivot_longer(c(-timeseries,-observations),names_to="Model",values_to="predictions")

# Define timeseries
timeseries <- unique(pred_DoYoff$timeseries)

# Initialize dataframe to store results
performance_stats <- data.frame()

# Loop across timeseries
for(ts in timeseries) {
  
  # Subset dataset according to timeseries
  sub_ts <- pred_DoYoff %>% 
    filter(timeseries==ts)
  
  # Calculate R2 per Model
  sub_R2 <- sub_ts
  sub_R2$Model <- paste0("R2_",sub_R2$Model)
  sub_R2$Model <- factor(sub_R2$Model, levels=paste0("R2_",model.names))
  sub_R2 <- sub_R2 %>%
    select(-timeseries) %>% 
    split(.$Model) %>%
    map(~lm(observations~predictions, data=.x)) %>% 
    map(glance) %>% 
    map_df("r.squared") 
  
  # Calculate RMSE
  sub_RMSE <- sub_ts
  sub_RMSE$Model <- paste0("RMSE_",sub_RMSE$Model)
  sub_RMSE$Model <- factor(sub_RMSE$Model, levels=paste0("RMSE_",model.names))
  sub_RMSE <- sub_RMSE %>%
    select(-timeseries) %>% 
    split(.$Model) %>% 
    map_df(~rmse(.$predictions,.$observations))
  
  # Calculate slope
  sub_slope <- sub_ts
  sub_slope$Model <- paste0("slope_",sub_slope$Model)
  sub_slope$Model <- factor(sub_slope$Model, levels=paste0("slope_",model.names))
  sub_slope <- sub_slope %>%
    select(-timeseries) %>% 
    split(.$Model) %>%
    map(~lm(observations~predictions, data=.x)) %>% 
    map(coefficients) %>% 
    map_df(2)

  # Data wrangling 
  sub_ts <- sub_R2 %>% 
    cbind(sub_RMSE) %>% 
    cbind(sub_slope) %>% 
    mutate(Species=strsplit(ts,"_")[[1]][2]) %>% 
    select(Species,everything()) %>% 
    mutate(timeseries=ts) %>% 
    select(timeseries,everything())
  
  # Store results
  performance_stats <- rbind(performance_stats,sub_ts)
  
  print(paste0("TIMESERIES EXECUTED: ",which(ts==timeseries)," OF ",length(timeseries)))
}

# Export dataset
write.table(performance_stats,"ModelAnalysis_3_Performance_R2_RMSE_slope.csv",sep=";",row.names=F)


##----------------------------------------
## Cross-Validation

## Clean working environment
rm(list=ls())

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

## Import data

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

# Optimal parameters per timeseries per model
opt_pars.df <- fread("ModelAnalysis_2_OptimalParameters.csv") 

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


## 5-fold Cross-Validation

# Initialize dataframe
Xval.df <- data.frame() 

for(sp in 1:length(species)) {
  
  # Initialize dataset per species
  Xval.sp <- data.frame()
  
  # Calculate number of sites per species
  sites <- SiteList[[sp]]
  
  for(site in 1:length(sites)) {
    
    # Subset according to timeseries
    ts <- paste0(species[sp],"_",sites[site])
    data.sub <- DataList[[sp]][[site]]
    preds.sub <- preds.df %>% 
      filter(timeseries==ts)
    opt_pars.sub <- opt_pars.df %>% 
      filter(timeseries==ts)
    
    # Initialize sub-dataframe to store results
    Xval.sub <- data.frame(timeseries=data.sub$site, Species=species[sp])
    
    # Initialize vectors to store predictions
    predictions_CDD <- vector(length = length(data.sub$site))
    predictions_DM1 <- vector(length = length(data.sub$site))
    predictions_DM2 <- vector(length = length(data.sub$site))
    predictions_TPM <- vector(length = length(data.sub$site))
    predictions_SIAM <- vector(length = length(data.sub$site))
    predictions_TDM <- vector(length = length(data.sub$site))
    predictions_TPDM <- vector(length = length(data.sub$site))
    predictions_PIAgsi <- vector(length = length(data.sub$site))
    `predictions_PIA-` <- vector(length = length(data.sub$site))
    `predictions_PIA+` <- vector(length = length(data.sub$site))
    
    # k-fold partitioning of the data set for model testing purposes. 
    # Each record is randomly assigned to a group. Group numbers are between 1 and k
    ks <- kfold(data.sub$year, k=5)
    for(i in 1:5){
      # get the train and test subsets
      train = list("doy" = data.sub$doy[ks!=i],
                   "site" = data.sub$site[ks!=i],
                   "location" = data.sub$location[,ks!=i], 
                   "year" = data.sub$year[ks!=i],
                   "Ti" = NULL,
                   "Tmini" = data.sub$Tmini[,ks!=i],
                   "Tmaxi" = NULL,
                   "Li" = data.sub$Li[,ks!=i],
                   "SPEI" = NULL,
                   "VPDi" = NULL,
                   "transition_dates" = data.sub$transition_dates[ks!=i],
                   "georeferencing" = NULL
      )
      test  = list("doy" = data.sub$doy[ks==i],
                   "site" = data.sub$site[ks==i],
                   "location" = data.sub$location[,ks==i], 
                   "year" = data.sub$year[ks==i],
                   "Ti" = NULL,
                   "Tmini" = data.sub$Tmini[,ks==i],
                   "Tmaxi" = NULL,
                   "Li" = data.sub$Li[,ks==i],
                   "SPEI" = NULL,
                   "VPDi" = NULL,
                   "transition_dates" = data.sub$transition_dates[ks==i],
                   "georeferencing" = NULL
      )
      
      ## Parameter optimization and Prediction of leaf senescence dates
      # PHENOR package (Hufkenset al., 2018)
      
      ## CDD model
      # Optimize parameters using the train dataset 
      # use opt_pars.sub to restrict the parameter range and speed up the validation
      optimal_pars <- optimize_parameters(par = NULL,
                                          data = train,
                                          cost = rmse,
                                          model = "CDD.model",
                                          method = "GenSA",
                                          lower = c(opt_pars.sub$Tbase_CDD-3,opt_pars.sub$Fcrit_CDD-10), 
                                          upper = c(opt_pars.sub$Tbase_CDD+3,0),
                                          control = list(max.call = 40000))
      # Predict dates of leaf senescence of the test dataset using the train-optimized parameters
      predictions_CDD[which(ks==i)] <- estimate_phenology(par = optimal_pars$par,
                                                          data = test,
                                                          model = "CDD.model")
      
      ## DM1 model
      # Optimize parameters using the train dataset 
      optimal_pars <- optimize_parameters(par = NULL,
                                          data = train,
                                          cost = rmse,
                                          model = "DM1.model",
                                          method = "GenSA",
                                          lower = c(opt_pars.sub$Tbase_DM1-3,opt_pars.sub$Pbase_DM1-3,opt_pars.sub$Fcrit_DM1-10), 
                                          upper = c(opt_pars.sub$Tbase_DM1+3,opt_pars.sub$Pbase_DM1-3,0),
                                          control = list(max.call = 40000))
      # Predict dates of leaf senescence of the test dataset using the train-optimized parameters
      predictions_DM1[which(ks==i)] <- estimate_phenology(par = optimal_pars$par,
                                                          data = test,
                                                          model = "DM1.model")
      
      ## DM2 model
      # Optimize parameters using the train dataset 
      optimal_pars <- optimize_parameters(par = NULL,
                                          data = train,
                                          cost = rmse,
                                          model = "DM2.model",
                                          method = "GenSA",
                                          lower = c(opt_pars.sub$Tbase_DM2-3,opt_pars.sub$Pbase_DM2-3,opt_pars.sub$Fcrit_DM2-10), 
                                          upper = c(opt_pars.sub$Tbase_DM2+3,opt_pars.sub$Pbase_DM2-3,0),
                                          control = list(max.call = 40000))
      # Predict dates of leaf senescence of the test dataset using the train-optimized parameters
      predictions_DM2[which(ks==i)] <- estimate_phenology(par = optimal_pars$par,
                                                          data = test,
                                                          model = "DM2.model")
      
      ## TPM model
      # Optimize parameters using the train dataset 
      optimal_pars <- optimize_parameters(par = NULL,
                                          data = train,
                                          cost = rmse,
                                          model = "TPM.model",
                                          method = "GenSA",
                                          lower = c(opt_pars.sub$Pbase_TPM-3,opt_pars.sub$a_TPM-0.02,opt_pars.sub$b_TPM-10,opt_pars.sub$Fcrit_TPM-10), 
                                          upper = c(opt_pars.sub$Pbase_TPM+3,opt_pars.sub$a_TPM+0.02,opt_pars.sub$b_TPM+10,0),
                                          control = list(max.call = 40000))
      # Predict dates of leaf senescence of the test dataset using the train-optimized parameters
      predictions_TPM[which(ks==i)] <- estimate_phenology(par = optimal_pars$par,
                                                          data = test,
                                                          model = "TPM.model")
      
      # Divide predictors between test and train datasets
      train_id <- train$site
      test_id <- test$site
      DoY_out.train <- preds.sub[which(preds.sub$timeseries==train_id),]$DoY_out
      DoY_out.test <- preds.sub[which(preds.sub$timeseries==test_id),]$DoY_out
      temp_GS.train <- preds.sub[which(preds.sub$timeseries==train_id),]$temp_GS
      temp_GS.test <- preds.sub[which(preds.sub$timeseries==test_id),]$temp_GS
      RD_summer.train <- preds.sub[which(preds.sub$timeseries==train_id),]$RD_summer
      RD_summer.test <- preds.sub[which(preds.sub$timeseries==test_id),]$RD_summer
      cGSI.train <- preds.sub[which(preds.sub$timeseries==train_id),]$cGSI
      cGSI.test <- preds.sub[which(preds.sub$timeseries==test_id),]$cGSI
      `cA_tot-w.train` <- preds.sub[which(preds.sub$timeseries==train_id),]$`cA_tot-w`
      `cA_tot-w.test` <- preds.sub[which(preds.sub$timeseries==test_id),]$`cA_tot-w`
      cA_tot.train <- preds.sub[which(preds.sub$timeseries==train_id),]$cA_tot
      cA_tot.test <- preds.sub[which(preds.sub$timeseries==test_id),]$cA_tot
      
      # Error check
      if(DoY_out.train!=train$doy) {
        print("ERROR k-fold partitioning!")
        break
      }
      
      ## SIAM model
      # Optimize parameters using the train dataset 
      optimal_pars <- optimize_parameters(par = NULL,
                                          predictor = DoY_out.train,
                                          data = train,
                                          cost = rmse,
                                          model = "SecondGen_PIA.models",
                                          method = "GenSA",
                                          lower = c(opt_pars.sub$Pbase_SIAM-3,opt_pars.sub$a_SIAM-0.02,opt_pars.sub$b_SIAM-10,opt_pars.sub$c_SIAM-20,opt_pars.sub$d_SIAM-0.3), 
                                          upper = c(opt_pars.sub$Pbase_SIAM+3,opt_pars.sub$a_SIAM+0.02,opt_pars.sub$b_SIAM+10,opt_pars.sub$c_SIAM+20,opt_pars.sub$d_SIAM+0.3),
                                          control = list(max.call = 40000))
      # Predict dates of leaf senescence of the test dataset using the train-optimized parameters
      predictions_SIAM[which(ks==i)] <- estimate_phenology(par = optimal_pars$par,
                                                           predictor = DoY_out.test,
                                                           data = test,
                                                           model = "SecondGen_PIA.models")
      
      ## TDM model
      # Optimize parameters using the train dataset 
      optimal_pars <- optimize_parameters(par = NULL,
                                          predictor = temp_GS.train,
                                          data = train,
                                          cost = rmse,
                                          model = "SecondGen_PIA.models",
                                          method = "GenSA",
                                          lower = c(opt_pars.sub$Pbase_TDM-3,opt_pars.sub$a_TDM-0.02,opt_pars.sub$b_TDM-10,opt_pars.sub$c_TDM-20,opt_pars.sub$d_TDM-0.3), 
                                          upper = c(opt_pars.sub$Pbase_TDM+3,opt_pars.sub$a_TDM+0.02,opt_pars.sub$b_TDM+10,opt_pars.sub$c_TDM+20,opt_pars.sub$d_TDM+0.3),
                                          control = list(max.call = 40000))
      # Predict dates of leaf senescence of the test dataset using the train-optimized parameters
      predictions_TDM[which(ks==i)] <- estimate_phenology(par = optimal_pars$par,
                                                          predictor = temp_GS.test,
                                                          data = test,
                                                          model = "SecondGen_PIA.models")
      
      ## TPDM model
      # Optimize parameters using the train dataset 
      optimal_pars <- optimize_parameters(par = NULL,
                                          predictor = c(temp_GS.train,RD_summer.train),
                                          data = train,
                                          cost = rmse,
                                          model = "SecondGen_PIA.models",
                                          method = "GenSA",
                                          lower = c(opt_pars.sub$Pbase_TPDM-3,opt_pars.sub$a_TPDM-0.02,opt_pars.sub$b_TPDM-10,opt_pars.sub$c_TPDM-20,opt_pars.sub$d_TPDM-0.3), 
                                          upper = c(opt_pars.sub$Pbase_TPDM+3,opt_pars.sub$a_TPDM+0.02,opt_pars.sub$b_TPDM+10,opt_pars.sub$c_TPDM+20,opt_pars.sub$d_TPDM+0.3),
                                          control = list(max.call = 40000))
      # Predict dates of leaf senescence of the test dataset using the train-optimized parameters
      predictions_TPDM[which(ks==i)] <- estimate_phenology(par = optimal_pars$par,
                                                           predictor = c(temp_GS.test,RD_summer.test),
                                                           data = test,
                                                           model = "SecondGen_PIA.models")
      
      ## PIA_gsi model
      # Optimize parameters using the train dataset 
      optimal_pars <- optimize_parameters(par = NULL,
                                          predictor = cGSI.train,
                                          data = train,
                                          cost = rmse,
                                          model = "SecondGen_PIA.models",
                                          method = "GenSA",
                                          lower = c(opt_pars.sub$Pbase_PIAgsi-3,opt_pars.sub$a_PIAgsi-0.02,opt_pars.sub$b_PIAgsi-10,opt_pars.sub$c_PIAgsi-20,opt_pars.sub$d_PIAgsi-0.3), 
                                          upper = c(opt_pars.sub$Pbase_PIAgsi+3,opt_pars.sub$a_PIAgsi+0.02,opt_pars.sub$b_PIAgsi+10,opt_pars.sub$c_PIAgsi+20,opt_pars.sub$d_PIAgsi+0.3),
                                          control = list(max.call = 40000))
      # Predict dates of leaf senescence of the test dataset using the train-optimized parameters
      predictions_PIAgsi[which(ks==i)] <- estimate_phenology(par = optimal_pars$par,
                                                             predictor = cGSI.test,
                                                             data = test,
                                                             model = "SecondGen_PIA.models")
      
      ## PIA- model
      # Optimize parameters using the train dataset 
      optimal_pars <- optimize_parameters(par = NULL,
                                          predictor = `cA_tot-w.train`,
                                          data = train,
                                          cost = rmse,
                                          model = "SecondGen_PIA.models",
                                          method = "GenSA",
                                          lower = c(opt_pars.sub$`Pbase_PIA-`-3,opt_pars.sub$`a_PIA-`-0.02,opt_pars.sub$`b_PIA-`-10,opt_pars.sub$`c_PIA-`-20,opt_pars.sub$`d_PIA-`-0.3), 
                                          upper = c(opt_pars.sub$`Pbase_PIA-`+3,opt_pars.sub$`a_PIA-`+0.02,opt_pars.sub$`b_PIA-`+10,opt_pars.sub$`c_PIA-`+20,opt_pars.sub$`d_PIA-`+0.3),
                                          control = list(max.call = 40000))
      # Predict dates of leaf senescence of the test dataset using the train-optimized parameters
      `predictions_PIA-`[which(ks==i)] <- estimate_phenology(par = optimal_pars$par,
                                                             predictor = `cA_tot-w.test`,
                                                             data = test,
                                                             model = "SecondGen_PIA.models")
      
      ## PIA+ model
      # Optimize parameters using the train dataset 
      optimal_pars <- optimize_parameters(par = NULL,
                                          predictor = cA_tot.train,
                                          data = train,
                                          cost = rmse,
                                          model = "SecondGen_PIA.models",
                                          method = "GenSA",
                                          lower = c(opt_pars.sub$`Pbase_PIA+`-3,opt_pars.sub$`a_PIA+`-0.02,opt_pars.sub$`b_PIA+`-10,opt_pars.sub$`c_PIA+`-20,opt_pars.sub$`d_PIA+`-0.3), 
                                          upper = c(opt_pars.sub$`Pbase_PIA+`+3,opt_pars.sub$`a_PIA+`+0.02,opt_pars.sub$`b_PIA+`+10,opt_pars.sub$`c_PIA+`+20,opt_pars.sub$`d_PIA+`+0.3),
                                          control = list(max.call = 40000))
      # Predict dates of leaf senescence of the test dataset using the train-optimized parameters
      `predictions_PIA+`[which(ks==i)] <- estimate_phenology(par = optimal_pars$par,
                                                             predictor = cA_tot.test,
                                                             data = test,
                                                             model = "SecondGen_PIA.models")
      
    }
    
    # Get the RMSE
    Xval.sub$RMSE_CDD <- rmse(predictions_CDD,data.sub$transition_dates)
    Xval.sub$RMSE_DM1 <- rmse(predictions_DM1,data.sub$transition_dates)
    Xval.sub$RMSE_DM2 <- rmse(predictions_DM2,data.sub$transition_dates)
    Xval.sub$RMSE_TPM <- rmse(predictions_TPM,data.sub$transition_dates)
    Xval.sub$RMSE_SIAM <- rmse(predictions_SIAM,data.sub$transition_dates)
    Xval.sub$RMSE_TDM <- rmse(predictions_TDM,data.sub$transition_dates)
    Xval.sub$RMSE_TPDM <- rmse(predictions_TPDM,data.sub$transition_dates)
    Xval.sub$RMSE_PIAgsi <- rmse(predictions_PIAgsi,data.sub$transition_dates)
    Xval.sub$`RMSE_PIA-` <- rmse(`predictions_PIA-`,data.sub$transition_dates)
    Xval.sub$`RMSE_PIA+` <- rmse(`predictions_PIA+`,data.sub$transition_dates)
    
    # Bind site-datasets
    Xval.sp <- rbind(Xval.sp,Xval.sub)
    print(paste0("RUNNING: ",Xval.sub$timeseries," ",site," OF ",length(sites)))
  }
  # Bind final datasets
  Xval.df <- rbind(Xval.df,Xval.sp)
}
write.table(Xval.df,"ModelAnalysis_4_CrossValidation_RMSE.csv",sep=";",row.names = F)


##----------------------------------------
## Model Performance

## FIGURE 3A
# Observed (Obs_AnomDoYoff) vs.
# Predicted (Pred_AnomDoYoff)
# leaf senescence anomalies, i.e., as deviation from the mean observed leaf-out date at each site
# of the best-performing first-generation (TPM) [Lang et al. (2019)], 
# second-generation (TPDM) [Liu et al. (2019)]
# and CarbLim models (PIAM)

# Import data
pred_DoYoff <- fread("ModelAnalysis_1_Predicted_DoYoff.csv")
pred_DoYoff <- pred_DoYoff %>% 
  filter(Species != "Betula pubescens")

# Select autumn anomalies from best-performing models
best_pred_DoYoff <- pred_DoYoff %>% 
  select(Obs_AnomDoYoff,Pred_AnomDoYoff_TPM,Pred_AnomDoYoff_TPDM,`Pred_AnomDoYoff_PIA+`)

# Format dataset for plotting
types <- c("First-generation (TPM)","Second-generation (TPDM)","CarbLim model (PIA+)")
colnames(best_pred_DoYoff) <- c("observations",types)
best_pred_DoYoff <- best_pred_DoYoff %>% 
  pivot_longer(-observations,names_to="model_type",values_to="predictions")
best_pred_DoYoff$model_type <- factor(best_pred_DoYoff$model_type, levels=types)

# Define palette
paletteBlueRed <- c("blue3","white","red3")

# Plot
fig_3a <- ggplot(best_pred_DoYoff, aes(x=predictions, y=observations)) +
  stat_bin2d(bins=235) +
  labs(x = "Predicted autumn anomaly",
       y = "Observed autumn anomaly") +
  coord_cartesian(xlim=c(-50,50), ylim=c(-50,50))+
  scale_fill_gradientn(colours = paletteBlueRed,
                       limits = c(10,2300),
                       breaks = c(100,1200,2200)) +
  geom_abline(slope=1, intercept=0,
              na.rm = FALSE, show.legend = NA, linetype="dashed") +
  theme(aspect.ratio=1,
        legend.position = "right",
        plot.title=element_text(hjust=.5),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=11),
        axis.title.x=element_text(size=14),
        axis.text.x=element_text(size=11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) + 
  facet_wrap(.~model_type, ncol=3)
fig_3a

# Add R^2 values
# Average across time-series
# See Figure 2B
dat_text <- data.frame(
  label = c("R2 = 0.21", "R2 = 0.62", "R2 = 0.78"), #check R2_models (Figure 3B)
  model_type = types
)
fig_3a  <- fig_3a + geom_text(
  data = dat_text,
  mapping = aes(x=-35, y=45, label=label)
)

# Calculate intercept plus slope by Standard Major Axis (SMA)
SMA_values <- best_pred_DoYoff %>%  
  group_by(model_type) %>% 
  summarise(int=lmodel2(observations ~ predictions)$regression.results$Intercept[3],
            slope=lmodel2(observations ~ predictions)$regression.results$Slope[3])
fig_3a <- fig_3a +
  geom_abline(data = SMA_values,
              mapping = aes(intercept=int,slope=slope),
              linetype = "solid")
fig_3a


## FIGURE 3B
# R^2 values

# Define auxiliary function to calculate standard error
se <- function(x) sqrt(var(x)/length(x))

# Import data
stat_models <- fread("ModelAnalysis_3_Performance_R2_RMSE_slope.csv")
stat_models <- stat_models %>% 
  filter(Species != "Betula pubescens")

# Define model names
model.names   <- c("CDD","DM1","DM2","TPM",
                   "SIAM","TDM","TPDM",
                   "PIA_gsi","PIA+")

# Define model types
types = c(rep("First-generation",4),
          rep("Second-generation",3),
          rep("PIA models",2))

# Select R^2 values
R2_models <- stat_models %>% 
  select(R2_CDD,R2_DM1,R2_DM2,R2_TPM,
         R2_SIAM,R2_TDM,R2_TPDM,
         R2_PIAgsi,`R2_PIA+`)

# Format dataset for plotting
colnames(R2_models) <- c(model.names)
R2_models <- R2_models %>% 
  pivot_longer(model.names,names_to="Model",values_to="value")
R2_models$Model <- factor(R2_models$Model, levels=model.names)
R2_models <- R2_models %>% 
  group_by(Model) %>% 
  summarise(R2_mean = mean(value),
            R2_se = se(value)) %>% 
  mutate(Type=types)
R2_models$Type <- factor(R2_models$Type, levels=c("First-generation","Second-generation","PIA models"))

# Define palette
trio <- c("#56B4E9","#E69F00","#74C476")

# Plot
fig_3b <- ggplot(R2_models, aes(x=Model, y=R2_mean)) +
  geom_bar(position=position_dodge(), stat="identity", width=.9,
           fill=c("#56B4E9","#56B4E9","#56B4E9","#56B4E9",
                  "#E69F00","#E69F00","#E69F00",
                  "#74C476","#74C476")) +
  geom_errorbar(aes(ymin=R2_mean-R2_se, ymax=R2_mean+R2_se),
                width=.2,
                position=position_dodge(.9)) +
  labs(y = expression(paste("Coefficient of determination ( ", R^2, ")"))) +
  scale_x_discrete(limits=model.names) +
  scale_color_manual(values=trio) +
  coord_cartesian(ylim=c(0.045,0.955))+
  theme(aspect.ratio = 0.3,
        axis.title.y=element_text(size=14, vjust=1),
        axis.text.y=element_text(size=11),
        axis.ticks.y=element_line(size=.3),
        axis.title.x=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.border = element_rect(fill=NA,colour = "black")
  )
fig_3b


## FIGURE 3C
# RMSE values
# for predictive models (solid-color-boxes)
# cross-validation at the time-series level (solid-black-lines)
# and for NULL model (dashed-black-line)

# Select RMSE values
RMSE_models <- stat_models %>% 
  select(timeseries,RMSE_CDD,RMSE_DM1,RMSE_DM2,RMSE_TPM,
         RMSE_SIAM,RMSE_TDM,RMSE_TPDM,
         RMSE_PIAgsi,`RMSE_PIA+`)

# Import RMSE values from cross-validation
RMSE_xval <- fread("ModelAnalysis_4_CrossValidation_RMSE.csv")
RMSE_xval <- RMSE_xval %>% 
  filter(Species!="Betula pubescens")
RMSE_xval <- RMSE_xval %>% 
  select(-timeseries,-`RMSE_PIA-`)
colnames(RMSE_xval) <- c("Species",model.names)

# Calculate average xval-RMSE across timeseries
RMSE_xval <- RMSE_xval %>% 
  pivot_longer(-Species,names_to="Model",values_to="rmse_xval")
RMSE_xval$Model <- factor(RMSE_xval$Model, levels=model.names)
RMSE_xval <- RMSE_xval %>% 
  group_by(Model) %>% 
  summarise(xval_median=median(rmse_xval),
            xval_mean=mean(rmse_xval),
            xval_sd=sd(rmse_xval))

# Format dataset for plotting
colnames(RMSE_models) <- c("timeseries",model.names)
RMSE_models <- RMSE_models %>% 
  pivot_longer(-timeseries,names_to="Model",values_to="value")
RMSE_models$Model <- factor(RMSE_models$Model, levels=model.names)
RMSE_models <- RMSE_models %>% 
  group_by(Model) %>% 
  summarise(Min=min(value),
            Q1=quantile(value,.25),
            Avg=mean(value),
            Median=median(value),
            Q3=quantile(value,.75),
            Max=max(value)) %>% 
  mutate(xval_median=RMSE_xval$xval_median,
         xval_mean=RMSE_xval$xval_mean,
         xval_sd=RMSE_xval$xval_sd) %>% 
  mutate(Type=types)
RMSE_models$Type <- factor(RMSE_models$Type, levels=c("First-generation","Second-generation","PIA models"))

# Calculate RMSE according to null model
RMSE_null <- Metrics::rmse(rep(mean(pred_DoYoff$Obs_DoYoff),nrow(pred_DoYoff)), pred_DoYoff$Obs_DoYoff)

# Define palette
trio <- c("#56B4E9","#E69F00","#74C476")

# Plot
fig_3c <- ggplot(RMSE_models, aes(x=Model, y=xval_median, color=Type)) +
  # RMSE values for predictive models
  geom_errorbar(aes(ymin=Min,ymax=Max),linetype=1, width=0, size=0.4) + 
  geom_crossbar(aes(y=Median,ymin=Q1,ymax=Q3), fill="grey90", linetype=1, size=0.4) + 
  # RMSE values for cross-validation at the time-series level
  geom_crossbar(ymin=RMSE_models$xval_median,ymax=RMSE_models$xval_median, color="black", size=0.2) +
  labs(y="RMSE") +
  scale_x_discrete(limits=model.names) +
  scale_color_manual(values=trio) +
  coord_cartesian(ylim=c(0.045,14.955))+
  # RMSE for NULL model
  geom_hline(aes(yintercept=RMSE_null), linetype="dashed") +
  theme(aspect.ratio = 0.3,
        axis.title.y=element_text(size=14, vjust=1),
        axis.text.y=element_text(size=11),
        axis.ticks.y=element_line(size=.4),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black')
  )
fig_3c


## FIGURE 3D
# slope values

# Select slope values
slope_models <- stat_models %>%
  select(slope_CDD,slope_DM1,slope_DM2,slope_TPM,
         slope_SIAM,slope_TDM,slope_TPDM,
         slope_PIAgsi,`slope_PIA+`)

# Format dataset for plotting
colnames(slope_models) <- c(model.names)
slope_models <- slope_models %>%
  pivot_longer(model.names,names_to="Model",values_to="value")
slope_models$Model <- factor(slope_models$Model, levels=model.names)
slope_models <- slope_models %>%
  mutate(Type=rep(types,nrow(stat_models)))
slope_models$Type <- factor(slope_models$Type, levels=c("First-generation","Second-generation","PIA models"))

# Plot
fig_3d <- ggplot(slope_models, aes(x=Model, y=value, color=Type)) +
  geom_boxplot(fill="grey90",outlier.shape = NA) +
  labs(y="Slope values") +
  scale_x_discrete(limits=model.names) +
  scale_color_manual(values=trio) +
  coord_cartesian(ylim=c(0.045,1.955))+
  geom_hline(aes(yintercept=1), linetype="dashed") +
  theme(aspect.ratio = 0.3,
        axis.title.y=element_text(size=14, vjust=1),
        axis.text.y=element_text(size=11),
        axis.ticks.y=element_line(size=.3),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, size=11),
        axis.ticks.x=element_line(size=.3),
        legend.title=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x=unit(0.10,"cm"),
        legend.key=element_blank(),
        legend.position="bottom",
        legend.direction="horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black')
  )
fig_3d

# Combine and export image
library(egg)
library(grid)

jpeg("Figure3.jpeg", width = 5.8*4, height = 5.8*5, units = "cm", res=600)
ggarrange(fig_3a,fig_3b,fig_3c,fig_3d,
          nrow=4)
grid.text("a", x=unit(0.94, "npc"), y=unit(0.99, "npc"),
          gp=gpar(fontface="bold", fontsize=18))
grid.text("b", x=unit(0.94, "npc"), y=unit(0.71, "npc"),
          gp=gpar(fontface="bold", fontsize=18))
grid.text("c", x=unit(0.94, "npc"), y=unit(0.5, "npc"),
          gp=gpar(fontface="bold", fontsize=18))
grid.text("d", x=unit(0.94, "npc"), y=unit(0.285, "npc"),
          gp=gpar(fontface="bold", fontsize=18))
dev.off()


## SUPPLEMENTARY FIGURE 6 (Fig. S6)

## Model performance across time. 
## Root-mean-square errors (RMSEs) of 
## observed versus predicted autumn anomalies
predictions.df <- fread("ModelAnalysis_1_Predicted_DoYoff.csv")

# Best second-generation model
rmse_TPDM.df <- predictions.df %>%
  group_by(YEAR) %>% 
  summarise(value=Metrics::rmse(Obs_AnomDoYoff,Pred_AnomDoYoff_TPDM)) %>% 
  mutate(roll_value=caTools::runmean(value,window_width),
         roll_sd=caTools::runsd(value,window_width),
         Model="Second-generation (TPDM)")

# Best PIA model
rmse_PIA.df <- predictions.df %>%
  group_by(YEAR) %>% 
  summarise(value=Metrics::rmse(Obs_AnomDoYoff,`Pred_AnomDoYoff_PIA+`)) %>% 
  mutate(roll_value=caTools::runmean(value,window_width),
         roll_sd=caTools::runsd(value,window_width),
         Model="CarbLim Model (PIA+)") 
rmse_to_plot.df <- rbind(rmse_TPDM.df,rmse_PIA.df)
rmse_to_plot.df$Model <- factor(rmse_to_plot.df$Model, levels=c("First-generation (TPM)","Second-generation (TPDM)","CarbLim Model (PIA+)"))

# Plot RMSEs
fig_S6 <- ggplot(rmse_to_plot.df, aes(x=YEAR, y=roll_value, color=Model, linetype=Model))+
  labs(x = "Year", y = "RMSE (days)") +
  geom_line(size=0.75) +
  scale_linetype_manual(values=c("solid", "solid"))+
  geom_ribbon(aes(ymin=roll_value-roll_sd, ymax=roll_value+roll_sd, group=Model, fill=Model), linetype=0, alpha=0.07) +
  scale_color_manual(values=c("#E69F00","#74C476")) +
  scale_fill_manual(values=c("#E69F00","#74C476")) +
  coord_cartesian(xlim=c(min(rmse_to_plot.df$YEAR),max(rmse_to_plot.df$YEAR)),
                  ylim=c(0,15)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(aspect.ratio = 1/2.5, 
        legend.position=c(0.2,0.2),
        legend.background = element_rect(fill=NA, 
                                         size=0.5, linetype = "solid"),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        strip.text=element_text(size=10),
        strip.background=element_rect(size=0.35),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(linetype = "solid",fill=NA,colour = 'black',size = .3),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=12,vjust=1),
        axis.title.y=element_text(size=12,vjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black')
  ) 
fig_S6


##----------------------------------------
## Summary statistics

# Load datasets
stat_models <- fread("ModelAnalysis_3_Performance_R2_RMSE_slope.csv")
stat_models <- stat_models %>% 
  filter(Species != "Betula pubescens")
RMSE_xval <- fread("ModelAnalysis_4_CrossValidation_RMSE.csv")
RMSE_xval <- RMSE_xval %>% 
  filter(Species!="Betula pubescens")

# Calculate means and standard deviations across timeseries
stat.df <- stat_models[complete.cases(stat_models[,-c(1,2)]),-c(1,2)]
colMeans(stat.df)
sapply(stat.df, sd)
sapply(stat.df, se)
colMeans(RMSE_xval[,-c(1,2)])
sapply(RMSE_xval[,-c(1,2)], sd)
sapply(RMSE_xval[,-c(1,2)], se)



##----------------------------------------
## References

# Dufrêne, E. et al. Modelling carbon and water cycles in a beech forest: Part I: Model description and uncertainty analysis on modelled NEE. Ecol. Modell. 185, 407-436 (2005).
# Delpierre, N. et al. Modelling interannual and spatial variability of leaf senescence for three deciduous tree species in France. Agric. For. Meteorol. 149, 938-948 (2009).
# Keenan, T. F. & Richardson, A. D. The timing of autumn senescence is affected by the timing of spring phenology: Implications 434 for predictive models. Glob. Chang. Biol. 21, 2634-2641 (2015).
# Lang, W., Chen, X., Qian, S., Liu, G. & Piao, S. A new process-based model for predicting autumn phenology: How is leaf senescence controlled by photoperiod and temperature coupling? Agric. For. Meteorol. 268, 124-135 (2019).
# Liu, G., Chen, X., Fu, Y. & Delpierre, N. Modelling leaf coloration dates over temperate China by considering effects of leafy season climate. 460 Ecol. Modell. 394, 34-43 (2019).
# Hufkens, K., Basler, D., Milliman, T., Melaas, E. K. & Richardson, A. D. An integrated phenology modelling framework in r. Methods Ecol. Evol. 9, 1276-1285 (2018).
# Running, S., Mu, Q., Zhao, M. MOD17A3 MODIS/Terra Net Primary Production Yearly L4 Global 1 km SIN Grid V055. NASA EOSDIS L. Process. DAAC (2011).
# Hansen, M. C., et al. High-Resolution Global Maps of 21st-Century Forest Cover Change. Science (80-. ). 850, 2011-2014 (2013).