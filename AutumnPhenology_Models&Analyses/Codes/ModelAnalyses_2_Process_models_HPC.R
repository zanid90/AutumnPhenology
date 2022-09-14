########################################
## LEAF SENESCENCE models (and drivers):
########################################

## First-generation:
####################

# CDD (chilling temperature) - Dufrene et al. (2005)
# DM1 and DM2 (chilling temperature, autumn daylength) - Delpierre et al. (2009)
# TPM (chilling temperature, autumn daylength) - Lang et al. (2019)

## Second-generation:
#####################

# SIAM (chilling temperature, autumn daylength, spring anomaly) - Keenan and Richardson (2015)
# TDM and TPDM (chilling temperature, autumn daylength, growing season temperature / + water stress) - Liu et al. (2019)

## PIA models:
##############

# PIA (chilling temperature, autumn daylength, photosynthesis [LPJ model])
# PIAgsi (chilling temperature, autumn daylength, photosynthesis [growing-season index])


# Define functions of Autumn phenology Models
# Modified from https://github.com/khufkens/phenor/blob/master/R/phenology_models.R



#######################################################################################################################
#######################################################################################################################



#####################
# Required packages #
#####################



# Load libraries
library(data.table)
library(tidyverse)
#devtools::install_github("bluegreen-labs/phenor@v1.0")
library(phenor)
library(parallel)



#######################################################################################################################
#######################################################################################################################



################
# Define Paths #
################



# 1. Input
##########

Input_path    = "Input"


# 2. Output
###########

output_path1   = "Output/Predictions"
output_path2   = "Output/Parameters"



#######################################################################################################################
#######################################################################################################################



###################
# Model functions #
###################



################
# 1. CDD model #
################


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


################################################################################################################


################
# 2. DM1 model #
################


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


################################################################################################################


################
# 3. DM2 model #
################


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


################################################################################################################


################
# 4. TPM model #
################


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


################################################################################################################


###################
# 5. Ycrit models #
###################


########################
# 5.1. TPM Ycrit model #
########################


Ycrit.TPM = function(data, P_base, a, b){
  
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
  
  # add observed leaf-off dates at the end of the matrix-column
  Rf = rbind(Rf,data$transition_dates)
  
  # calculate the summation along the year and derive the date of leaf.off
  # DOY of budburst criterium
  F_crit = apply(Rf,2, function(xt){
    F_crit = cumsum(xt[1:366])[xt[367]]
    return(F_crit)
  })
  
  return(F_crit)
}


#####################
# 5. DM Ycrit model #
#####################


Ycrit.DM = function(data, P_base, T_base){
  
  # create forcing/chilling rate vector at the day level
  Rf = (data$Tmini - T_base)*(data$Li/P_base)
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
  
  # add observed leaf-off dates at the end of the matrix-column
  Rf = rbind(Rf,data$transition_dates)
  
  # calculate the summation along the year and derive the date of leaf.off
  # DOY of budburst criterium
  F_crit = apply(Rf,2, function(xt){
    F_crit = cumsum(xt[1:366])[xt[367]]
    return(F_crit)
  })
  
  return(F_crit)
}



################################################################################################################


###############################
# 6. Second generation models #
###############################


####################################
# 6.1. Second generation TPM model #
####################################


SecondGen_PIA.models.TPM = function(par, predictor, data) {
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
    t0A = interval[which(data$Li[,col] < P_base)]
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


###################################
# 6.2. Second generation DM model #
###################################


SecondGen_PIA.models.DM = function(par, predictor, data) {
  # exit the routine as some parameters are missing
  if (length(par) != 4 & length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  T_base = as.numeric(par[1])
  P_base = as.numeric(par[2])
  c = as.numeric(par[3])
  if(length(par)==4) {
    d = as.numeric(par[4])
    pred = predictor
  }
  if(length(par)==5) {
    d = as.numeric(par[4])
    e = as.numeric(par[5])
    pred1 = predictor[1]
    pred2 = predictor[2]
  }
  
  # create forcing/chilling rate vector at the day level
  Rf = (data$Tmini - T_base)*(data$Li/P_base)
  Rf[Rf > 0] = 0
  
  # photoperiod-dependent start-date for chilling accumulation (t0)
  # t0 is defined as the first day when photoperiod is shorter than the photoperiod threshold (P_base)
  # after the date of the longest photoperiod (summer solstice), namely, the 173rd day of year
  t0 <- vector()
  for(col in 1:ncol(data$Tmini)) {
    interval = 1:366
    t0A = interval[which(data$Li[,col] < P_base)]
    ind1 = min(which(t0A > 173))
    t0A = t0A[ind1]
    t0 = c(t0,t0A)
  }
  
  # nullify values before the t0
  for(col in 1:ncol(data$Tmini)){
    Rf[1:t0[col],col] = 0 
  }
  
  if(length(par)==4) {
    # add predictor at the end of the matrix-columns
    Rf = rbind(Rf,predictor)
    
    # predict date of leaf.off
    doy = apply(Rf,2, function(xt){
      doy = which(cumsum(xt[1:366]) <= c+d*xt[367])[1]
    }) 
  }
  if(length(par)==5) {
    # add predictors at the end of the matrix-columns
    Rf = rbind(Rf,pred1)
    Rf = rbind(Rf,pred2)
    
    # predict date of leaf.off
    doy = apply(Rf,2, function(xt){
      doy = which(cumsum(xt[1:366]) <= c+d*xt[367]+e*xt[368])[1]
    }) 
  }
  
  return(doy)
}



#######################################################################################################################
#######################################################################################################################



########################
# Import & Format data #
########################



#PEP drivers data
#################

all.df <- fread(paste(Input_path,"DataMeta_3_Drivers.csv",sep="/")) %>%
  #add timeseries length
  group_by(timeseries) %>%
  add_tally()%>%
  ungroup()


# Senescence observations
#########################

pheno.df <- all.df %>% 
  dplyr::select("timeseries","pep_id","lon","lat","species","year","leaf_off","leaf_out","n") %>%
  mutate(site_yr = paste0(pep_id, "_", year))%>%
  #sort from longest to shortest timeseries
  arrange(-n, timeseries)


# Predictors
############

preds.df <- all.df %>% 
  dplyr::select("timeseries","pep_id","lon","lat","species","year",
                "leaf_out","Tls.Off","LPI",
                "GSI.Off","Anet.Off","AnetNW.Off")
# leaf_out: leaf-out dates (Keenan et al. 2017)
# Tls: growing-season temperature (Liu et al. 2018)
# LPI: low precipitation index
# GSI: Growing-season Index
# Anet: Net photosynthesis
# AnetNW: Net photosynthesis excluding water stress function (no water)


# Daily minimum temperature
###########################

tmin.df <- fread(paste(Input_path,"Daily_Min_Data_Tair_f_inst.csv",sep="/")) %>%
  #add identifier
  mutate(site_yr = paste0(pep_id,"_", year)) %>%
  #remove columns
  dplyr::select(-c(year,pep_id,lat,lon))


# Photoperiod file
##################

photo.df <- fread(paste(Input_path,"Photoperiod.csv",sep="/"))%>%
  dplyr::select(-c(lat))


#######################################################################################################################


#merge tables
phenor_input = merge(pheno.df, tmin.df,  by="site_yr")
phenor_input = merge(phenor_input, photo.df, by="pep_id")

#delete columns from pheno data
pheno.df = pheno.df %>%
  dplyr::select("timeseries","pep_id","lon","lat","species","year","leaf_off")

#removse stuff
rm(all.df,tmin.df,photo.df)



#######################################################################################################################
#######################################################################################################################



############################
# Create phenoR input list #
############################



# Define timeseries
timeseries <- unique(pheno.df$timeseries)

# Create a DataList to store all subsets for each timeseries
DataList <- replicate(length(timeseries), data.frame())
names(DataList) <- paste0(timeseries)

# Counter for timeseries
count  <- 0

# Initiate an external loop to subset for each timeseries
for(timeseries.i in timeseries) {
  
  # Subset phenological data
  phenor.sub <- phenor_input[which(phenor_input$timeseries==timeseries.i),]
  
  data = list("doy" = as.vector(phenor.sub$leaf_out),
              "site" = timeseries.i,
              "location" = t(phenor.sub[,c("lon","lat")]), 
              "year" = as.vector(phenor.sub$year),
              "Ti" = NULL,
              "Tmini" = t(phenor.sub[,paste0(1:366,".x")]), 
              "Tmaxi" = NULL,
              "Li" = t(phenor.sub[,paste0(1:366,".y")]),
              "SPEI" = NULL,
              "VPDi" = NULL,
              "transition_dates" = as.vector(phenor.sub$leaf_off),
              "georeferencing" = NULL)
  
  # Store each site-specific dataframe in the DataList of the corresponding species
  count <- count+1
  DataList[[timeseries.i]] <- data
  #print(paste0(round(count/length(timeseries)*100,1),"% DONE!"))
}

rm(data,count,phenor_input)



#######################################################################################################################
#######################################################################################################################



###################
# Run Predictions #
###################



parallelCalc <- function(time.series){ 
  
  # Subset according to timeseries
  data.sub  <- DataList[[time.series]]
  preds.sub <- preds.df %>% filter(timeseries==time.series)
  
  # Initialize sub-dataframes to store results
  leafOff_Preds_PHENOR.sub <- pheno.df %>% filter(timeseries==time.series)
  opt_pars_PHENOR.sub <- preds.df %>% dplyr::select(timeseries, species, pep_id) %>% filter(row_number()==1)
  
  ###########################################################################################
  ###########################################################################################
  
  ## Parameter optimization and Prediction of leaf senescence dates
  # PHENOR package (Hufkens et al., 2018)
  
  ###########################################################################################
  ###########################################################################################
  
  ###############
  ## CDD model ##
  ###############
  
  optimal_pars <- optimize_parameters(par = NULL,
                                      data = data.sub,
                                      cost = rmse,
                                      model = "CDD.model",
                                      method = "GenSA",
                                      lower = c(10,-3000),
                                      upper = c(30,0)
                                      ,control = list(max.call = 40000)
  )
  
  opt_pars_PHENOR.sub$Tbase_CDD <- optimal_pars$par[1]
  opt_pars_PHENOR.sub$Fcrit_CDD <- optimal_pars$par[2]
  
  leafOff_Preds_PHENOR.sub$Pred_leafOff_CDD <- estimate_phenology(par = optimal_pars$par,
                                                                  data = data.sub,
                                                                  model = "CDD.model")
  
  ###########################################################################################
  
  ###############
  ## DM1 model ##
  ###############
  
  optimal_pars <- phenor::optimize_parameters(par = NULL,
                                              data = data.sub,
                                              cost = rmse,
                                              model = "DM1.model",
                                              method = "GenSA",
                                              lower = c(10,11,-2000),
                                              upper = c(30,16,0)
                                              ,control = list(max.call = 40000)
  )
  
  opt_pars_PHENOR.sub$Tbase_DM1 <- optimal_pars$par[1]
  opt_pars_PHENOR.sub$Pbase_DM1 <- optimal_pars$par[2]
  opt_pars_PHENOR.sub$Fcrit_DM1 <- optimal_pars$par[3]
  
  best_pars.DM <- optimal_pars$par
  
  leafOff_Preds_PHENOR.sub$Pred_leafOff_DM1 <- estimate_phenology(par = optimal_pars$par,
                                                                  data = data.sub,
                                                                  model = "DM1.model")
  
  ###########################################################################################
  
  ###############
  ## DM2 model ##
  ###############
  
  optimal_pars <- optimize_parameters(par = NULL,
                                      data = data.sub,
                                      cost = rmse,
                                      model = "DM2.model",
                                      method = "GenSA",
                                      lower = c(10,11,-2000),
                                      upper = c(30,16,0)
                                      ,control = list(max.call = 40000)
  )
  
  opt_pars_PHENOR.sub$Tbase_DM2 <- optimal_pars$par[1]
  opt_pars_PHENOR.sub$Pbase_DM2 <- optimal_pars$par[2]
  opt_pars_PHENOR.sub$Fcrit_DM2 <- optimal_pars$par[3]
  
  leafOff_Preds_PHENOR.sub$Pred_leafOff_DM2 <- estimate_phenology(par = optimal_pars$par,
                                                                  data = data.sub,
                                                                  model = "DM2.model")
  
  ###########################################################################################
  
  ###############
  ## TPM model ##
  ###############
  
  optimal_pars <- optimize_parameters(par = NULL,
                                      data = data.sub,
                                      cost = rmse,
                                      model = "TPM.model",
                                      method = "GenSA",
                                      lower = c(11,0.02,100,0),
                                      upper = c(16,0.1,250,200)
                                      ,control = list(max.call = 40000)
  )
  
  opt_pars_PHENOR.sub$Pbase_TPM <- optimal_pars$par[1]
  opt_pars_PHENOR.sub$a_TPM     <- optimal_pars$par[2]
  opt_pars_PHENOR.sub$b_TPM     <- optimal_pars$par[3]
  opt_pars_PHENOR.sub$Fcrit_TPM <- optimal_pars$par[4]
  
  best_pars.TPM <- optimal_pars$par
  
  leafOff_Preds_PHENOR.sub$Pred_leafOff_TPM <- estimate_phenology(par = optimal_pars$par,
                                                                  data = data.sub,
                                                                  model = "TPM.model")
  
  ###########################################################################################
  
  
  ##############################
  ## Second generation models ##
  ##############################
  
  
  # Calculate Y_crit with fitted TPM parameters (best First-Generation model)
  Y_crit.TPM <- Ycrit.TPM(data=data.sub, P_base=opt_pars_PHENOR.sub$Pbase_TPM, a=opt_pars_PHENOR.sub$a_TPM, b=opt_pars_PHENOR.sub$b_TPM)
  
  
  ###########################################################################################
  
  
  ########################
  # One predictor variable
  ########################
  
  
  variables = c("Anet.Off", "AnetNW.Off", "GSI.Off", "Tls.Off", "leaf_out")
  models    = c("PIA",      "PIAminus",   "PIAgsi",  "TDM",     "SIAM")
  
  for (i in 1:length(variables)){
    
    # get predictor variable
    predictor = preds.sub %>% pull(variables[i])
    
    # Calculate scaled site anomalies
    Anomaly = as.numeric(scale(predictor))
    
    #---------------------------------------------------------------------------------------------------------
    
    ###########
    # TPM-based
    ###########
    
    fit        <- lm(Y_crit.TPM ~ Anomaly)
    intercept  <- fit[[1]][1]
    slope      <- fit[[1]][2]
    opt_pars   <- c(best_pars.TPM[1:3],intercept,slope)
    boundaries <- c(2,100,100,100,100)
    lower_pars <- c()
    upper_pars <- c()
    for(n in 1:length(boundaries)) {
      lower_pars[n] <- opt_pars[n] - boundaries[n]
      upper_pars[n] <- opt_pars[n] + boundaries[n]
    }
    
    optimal_pars <- optimize_parameters(par = opt_pars,
                                        predictor = Anomaly,
                                        data = data.sub,
                                        cost = rmse,
                                        model = "SecondGen_PIA.models.TPM",
                                        method = "GenSA",
                                        lower = lower_pars,
                                        upper = upper_pars
                                        ,control = list(max.call = 40000)
    )
    
    # Name variables
    varname.Pbase <- paste("Pbase", models[i], sep="_")
    varname.a     <- paste("a",     models[i], sep="_")
    varname.b     <- paste("b",     models[i], sep="_")
    varname.c     <- paste("c",     models[i], sep="_")
    varname.d     <- paste("d",     models[i], sep="_")
    varname.Pred  <- paste("Pred_leafOff", models[i], sep="_")
    
    # Create columns
    opt_pars_PHENOR.sub = opt_pars_PHENOR.sub %>%
      mutate(!!varname.Pbase := optimal_pars$par[1],
             !!varname.a     := optimal_pars$par[2],
             !!varname.b     := optimal_pars$par[3],
             !!varname.c     := optimal_pars$par[4],
             !!varname.d     := optimal_pars$par[5] )
    
    leafOff_Preds_PHENOR.sub = leafOff_Preds_PHENOR.sub %>%
      mutate(!!varname.Pred := estimate_phenology(par = optimal_pars$par,
                                                  predictor = Anomaly,
                                                  data = data.sub,
                                                  model = "SecondGen_PIA.models.TPM") )
  }
  
  
  ###########################################################################################
  
  
  ########################
  #Two predictor variables
  ########################
  
  
  # get predictor variables
  predictor1 = preds.sub %>% pull("Tls.Off")
  predictor2 = preds.sub %>% pull("LPI")
  
  # Calculate scaled site anomalies      
  Anomaly1 <- as.numeric(scale(predictor1))
  Anomaly2 <- as.numeric(scale(predictor2))
  
  #---------------------------------------------------------------------------------------------------------
  
  ###########
  # TPM-based
  ###########
  
  fit        <- lm(Y_crit.TPM ~ Anomaly1 + Anomaly2)
  intercept  <- fit[[1]][1]
  slope      <- fit[[1]][2]
  slope2     <- fit[[1]][3]
  opt_pars   <- c(best_pars.TPM[1:3],intercept,slope,slope2)
  boundaries <- c(2,100,100,100,100,100)
  lower_pars <- c()
  upper_pars <- c()
  for(n in 1:length(boundaries)) {
    lower_pars[n] <- opt_pars[n] - boundaries[n]
    upper_pars[n] <- opt_pars[n] + boundaries[n]
  }
  
  optimal_pars <- optimize_parameters(par = opt_pars,
                                      predictor = c(Anomaly1, Anomaly2),
                                      data = data.sub,
                                      cost = rmse,
                                      model = "SecondGen_PIA.models.TPM",
                                      method = "GenSA",
                                      lower = lower_pars,
                                      upper = upper_pars
                                      ,control = list(max.call = 40000)
  )
  
  # Create columns
  opt_pars_PHENOR.sub = opt_pars_PHENOR.sub %>%
    mutate(Pbase_TPDM = optimal_pars$par[1],
           a_TPDM     = optimal_pars$par[2],
           b_TPDM     = optimal_pars$par[3],
           c_TPDM     = optimal_pars$par[4],
           d_TPDM     = optimal_pars$par[5],
           e_TPDM     = optimal_pars$par[6])
  
  leafOff_Preds_PHENOR.sub = leafOff_Preds_PHENOR.sub %>%
    mutate(Pred_leafOff_TPDM = estimate_phenology(par = optimal_pars$par,
                                                  predictor = c(Anomaly1, Anomaly2),
                                                  data = data.sub,
                                                  model = "SecondGen_PIA.models.TPM") )
  
  
  ###########################################################################################
  ###########################################################################################
  
  
  # Safe the tables
  write.table(leafOff_Preds_PHENOR.sub, 
              file=paste0(output_path1, '/', time.series, '.csv'), sep=',', row.names =F, col.names = T)
  
  write.table(opt_pars_PHENOR.sub,      
              file=paste0(output_path2, '/', time.series, '.csv'), sep=',', row.names =F, col.names = T)
  
  print(paste0("time series: ",time.series," done"))
}



#######################################################################################################################
#######################################################################################################################



##################
## Run the Loop ##
##################



#initialize the loop
outputlist <- mclapply(timeseries, parallelCalc, mc.cores=24, mc.preschedule=F)



#######################################################################################################################
############################################################# THE END #################################################
#######################################################################################################################