############################
## Future Autumn Projections

# Define directory paths
setwd(".../AutumnPhenology/AutumnPhenology_Data&Metadata/Data/")

# Load libraries
library(data.table)
library(tidyverse)
library(caTools)

# Define functions of Autumn phenology Models
CDD.models = function(T_base, F_crit, data){  
  
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
    return(doy)
  })
  
  return(doy)
}
TPM.models = function(P_base, a, b, F_crit, data){
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
  
  # calculate the summation along the year (interval = 1:366) and derive the date of leaf.off
  # DOY of budburst criterium
  doy = apply(Rf,2, function(xt){
    doy = which(cumsum(xt) >= F_crit)[1]
  })
  
  return(doy)
}
TPM.Ycrit = function(data, P_base, a, b, stop){
  
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
  Rf = rbind(Rf,stop)
  
  # calculate the summation along the year and derive the date of leaf.off
  # DOY of budburst criterium
  F_crit = apply(Rf,2, function(xt){
    F_crit = cumsum(xt[1:366])[xt[367]]
    return(F_crit)
  })
  
  return(F_crit)
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
      return(doy)
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

# Helper function
se <- function(x) sqrt(var(x)/length(x))

# Which model for future spring phenology?
model <- "M1" # AT, M1, PM1

##----------------------------------------
## Import data

# Future spring phenology
future_spring.df <- fread(paste0("Future_SpringPhenology_soil_CO2_",model,".csv"))
future_spring.df <- future_spring.df %>% 
  filter(DOY_out < 250) %>% 
  select(-c(ETA)) 

# Future minimum temperature and photoperiod
future_tmin.df <- fread("Future Minimum Temperature.csv")
future_photo.df <- fread("FuturePhotoperiod.csv")
future_tmin.df$ts_yr <- paste0(future_tmin.df$PEP_ID,"_",future_tmin.df$YEAR)
future_photo.df$ts_yr <- paste0(future_photo.df$PEP_ID,"_",future_photo.df$YEAR)

# Predictors
# GST: growing season temperature
# RD_summer: number of rainy days (precipitation >2mm) during the driest months
# cGSI: cumulative Growing Season Index
# cAtot_w: cumulative net photosynthesis modified by water-stress factor
future_preds.df <- fread("DataMeta_6_FutureDrivers.csv") 

# Optimal parameters per timeseries per model
opt_pars.df <- fread("ModelAnalysis_2_OptimalParameters.csv") 

# Future Autumn phenology Models
models <- c("CDD","TPM","SIAM","TPDM","PIA_gsi","PIA+")

# Define species
species <- unique(future_spring.df$Species)

# Define timeseries
timeseries <- unique(future_spring.df$timeseries)

# Prepare input datasets for PHENOR package
future_tmin.df <- future_tmin.df %>% 
  select(-c(YEAR,PEP_ID,LAT,LON))
future_photo.df <- future_photo.df %>% 
  select(-c(YEAR,PEP_ID,LAT))
phenor_input <- merge(future_tmin.df,future_photo.df,by="ts_yr")
phenor_input <- merge(future_spring.df,phenor_input,by="ts_yr")
rm(future_spring.df,future_tmin.df,future_photo.df)

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
    data = list("doy" = as.vector(pheno_sp.sub[index,]$DOY_out),
                "site" = as.vector(paste0(pheno_sp.sub[index,]$PEP_ID,"_",pheno_sp.sub[1,]$Species)),
                "location" = t(pheno_sp.sub[index,c("LON","LAT")]), 
                "year" = as.vector(pheno_sp.sub[index,]$YEAR),
                "Ti" = NULL,
                "Tmini" = t(pheno_sp.sub[index,paste0(1:366,".x")]), 
                "Tmaxi" = NULL,
                "Li" = t(pheno_sp.sub[index,paste0(1:366,".y")]),
                "SPEI" = NULL,
                "VPDi" = NULL,
                "transition_dates" = NULL,
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
## Future projections

# Initialize dataframe
future_proj.df <- data.frame() 

for(sp in 1:length(species)) {
  
  # Initialize dataset per species
  future_proj.sp <- data.frame()
  
  # Calculate number of sites per species
  sites <- SiteList[[sp]]
  
  for(site in 1:length(sites)) {
    
    # Subset according to timeseries
    ts <- paste0(species[sp],"_",sites[site])
    data.sub <- DataList[[sp]][[site]]
    preds.sub <- future_preds.df %>% 
      filter(timeseries==ts)
    opt_pars.sub <- opt_pars.df %>% 
      filter(timeseries==ts)
    
    # Initialize sub-dataframe to store results
    future_proj.sub <- data.frame(timeseries=data.sub$site, YEAR=data.sub$year, PEP_ID=sites[site], Species=species[sp])
    
    ## Future leaf senescence dates (DoY.off)
    # CDD model
    future_proj.sub$DoY.off_CDD  <- CDD.models(T_base=opt_pars.sub$Tbase_CDD, F_crit=opt_pars.sub$Fcrit_CDD, data=data.sub)

    # TPM model
    future_proj.sub$DoY.off_TPM  <- TPM.models(P_base=opt_pars.sub$Pbase_TPM, a=opt_pars.sub$a_TPM, b=opt_pars.sub$b_TPM, F_crit=opt_pars.sub$Fcrit_TPM, data=data.sub)

    # Estimate future Ycrit with TPM parameters (best First-Generation model)
    Y_crit <- TPM.Ycrit(data=data.sub, P_base=opt_pars.sub$Pbase_TPM, a=opt_pars.sub$a_TPM, b=opt_pars.sub$b_TPM, stop=future_proj.sub$DoY.off_TPM)
    
    # SIAM model
    future_proj.sub$DoY.off_SIAM  <- opt_pars.sub$a_SIAM + opt_pars.sub$b_SIAM*Y_crit + opt_pars.sub$c_SIAM*data.sub$doy

    # TPDM model
    future_proj.sub$DoY.off_TPDM  <- opt_pars.sub$a_TPDM + opt_pars.sub$b_TPDM*Y_crit + opt_pars.sub$c_TPDM*preds.sub$GST + opt_pars.sub$a_TPDM*preds.sub$RD_summer

    # PIA_gsi model
    future_proj.sub$DoY.off_PIA_gsi <- opt_pars.sub$a_PIAgsi + opt_pars.sub$b_PIAgsi*Y_crit + opt_pars.sub$c_PIAgsi*preds.sub$cGSI

    # PIA+ model
    future_proj.sub$`DoY.off_PIA+`  <- opt_pars.sub$`a_PIA+` + opt_pars.sub$`b_PIA+`*Y_crit + opt_pars.sub$`c_PIA+`*preds.sub$cA_tot

    ## Future Growing Season Length (GSL)
    future_proj.sub$GSL_CDD <- future_proj.sub$DoY.off_CDD - data.sub$doy
    future_proj.sub$GSL_TPM <- future_proj.sub$DoY.off_TPM - data.sub$doy
    future_proj.sub$GSL_SIAM <- future_proj.sub$DoY.off_SIAM - data.sub$doy
    future_proj.sub$GSL_TPDM <- future_proj.sub$DoY.off_TPDM - data.sub$doy
    future_proj.sub$GSL_PIA_gsi <- future_proj.sub$DoY.off_PIA_gsi - data.sub$doy
    future_proj.sub$`GSL_PIA+` <- future_proj.sub$`DoY.off_PIA+` - data.sub$doy
    
    # Bind site-datasets
    future_proj.sp <- rbind(future_proj.sp,future_proj.sub)
    print(paste0("RUNNING: ",future_proj.sub$timeseries," ",site," OF ",length(sites)))
  }
  # Bind final datasets
  future_proj.df <- rbind(future_proj.df,future_proj.sp)
}
write.table(future_proj.df,"ModelAnalysis_5_FutureAutumnProjections.csv",sep=";",row.names = F)


##----------------------------------------
## Plot

# Define model names 
model.names   <- c("CDD","TPM","SIAM","TPDM","PIA_gsi","PIA+")

# Import data and format
future.df <- fread("ModelAnalysis_5_FutureAutumnProjections.csv") 
future.df <- future.df %>% 
  as.data.frame()
doyoff.df <- future.df[,grep("DoY.off",colnames(future.df))]
GSL.df <- future.df[,grep("GSL",colnames(future.df))]
future_doyoff.df <- doyoff.df %>% 
  mutate(YEAR=future.df$YEAR)
future_GSL.df <- GSL.df %>% 
  mutate(YEAR=future.df$YEAR)
colnames(future_doyoff.df) <- c(model.names,"YEAR")
colnames(future_GSL.df) <- c(model.names,"YEAR")

# Calculate average senescence date/growing season length across time-series
future_doyoff.df <- future_doyoff.df %>% 
  pivot_longer(-YEAR,names_to="Model",values_to="FutProj_DoYoff") %>% 
  group_by(Model,YEAR) %>% 
  summarise(mean_DoY=mean(FutProj_DoYoff,na.rm=T))
future_GSL.df <- future_GSL.df %>% 
  pivot_longer(-YEAR,names_to="Model",values_to="FutProj_GSL") %>% 
  group_by(Model,YEAR) %>% 
  summarise(mean_GSL=mean(FutProj_GSL,na.rm=T))

# Calculate senescence/growing season length projections as 15-year moving averages 
window_width <- 15
future_doyoff.df <- future_doyoff.df %>% 
  group_by(Model) %>% 
  mutate(roll_DoY=runmean(mean_DoY,window_width),
         sd_DoY=runsd(mean_DoY,window_width))
future_doyoff.df$Model <- factor(future_doyoff.df$Model,levels=model.names)
future_GSL.df <- future_GSL.df %>% 
  group_by(Model) %>% 
  mutate(roll_GSL=runmean(mean_GSL,window_width),
         sd_GSL=runsd(mean_GSL,window_width))
future_GSL.df$Model <- factor(future_GSL.df$Model,levels=model.names)

# Define palette
paletteSoftRainbow <- c("#ABDAF4","#2b5a74",
                        "#e6d900","#e66600",
                        "#C7E9C0","#006D2C")

# Plot future leaf senescence dates
fig_4a <- ggplot(future_doyoff.df, aes(x=YEAR, y=round(roll_DoY), colour=Model))+
  labs(x = "Year", y = "Leaf senescence (DoY)") +
  geom_line(size=.5) +
  coord_cartesian(xlim=c(2022,2096.5), ylim=c(260,307))+
  scale_colour_manual(values=paletteSoftRainbow) +
  theme(aspect.ratio = 1, 
        legend.position = "none",
        legend.background = element_rect(fill=NA, 
                                         size=0.5, linetype = "solid"),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(linetype = "solid",fill=NA,colour = 'black',size = .3),
        axis.ticks.y=element_line(size=.3),
        axis.ticks.x=element_blank(),
        axis.line.x=element_line(size=.3),
        axis.line.y=element_line(size=.3),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=6),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8,vjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black')
  ) + 
  geom_ribbon(aes(ymin=roll_DoY-sd_DoY, ymax=roll_DoY+sd_DoY, group=Model), 
              linetype=0, colour="grey", alpha=0.07) +
  geom_text(aes(x=2025,y=305,label="a"),size=4,color="black",fontface="bold")
fig_4a

# Plot future growing season length
fig_4c <- ggplot(future_GSL.df, aes(x=YEAR, y=round(roll_GSL), colour=Model))+
  labs(x = "Year", y = "Growing season length (days)") +
  geom_line(size=.5) +
  coord_cartesian(xlim=c(2022,2096.5), ylim=c(162,215))+
  scale_colour_manual(values=paletteSoftRainbow) +
  theme(aspect.ratio = 1, 
        legend.position = "none",
        legend.background = element_rect(fill=NA, 
                                         size=0.5, linetype = "solid"),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(linetype = "solid",fill=NA,colour = 'black',size = .3),
        axis.ticks.y=element_line(size=.3),
        axis.ticks.x=element_line(size=.3),
        axis.line.x=element_line(size=.3),
        axis.line.y=element_line(size=.3),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.title.x=element_text(size=8,vjust=1),
        axis.title.y=element_text(size=8,vjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black')
  ) + 
  geom_ribbon(aes(ymin=roll_GSL-sd_GSL, ymax=roll_GSL+sd_GSL, group=Model), 
              linetype=0, colour="grey", alpha=0.07) +
  geom_text(aes(x=2025,y=213,label="c"),size=4,color="black",fontface="bold")
fig_4c


## Estimated delays in leaf senescence by the end of the 21st century (2080-2100) 
# compared to the average senescence dates between 1990-2010

# Select future senescence dates/growing season length (2080-2100 range)
fut_DoYs <- doyoff.df %>% 
  mutate(timeseries=future.df$timeseries) %>% 
  mutate(Species=future.df$Species) %>% 
  mutate(YEAR=future.df$YEAR) %>% 
  filter(YEAR %in% 2080:2100) %>% 
  select(-YEAR)
colnames(fut_DoYs) <- c(model.names,"timeseries","Species")
fut_GSLs <- GSL.df %>% 
  mutate(timeseries=future.df$timeseries) %>% 
  mutate(Species=future.df$Species) %>% 
  mutate(YEAR=future.df$YEAR) %>% 
  filter(YEAR %in% 2080:2100) %>% 
  select(-YEAR)
colnames(fut_GSLs) <- c(model.names,"timeseries","Species")

# Calculate average senescence date/growing season length for the time range
fut_DoYs <- fut_DoYs %>%
  pivot_longer(-c(timeseries,Species),names_to="Model",values_to="fut_DoYoff")
fut_DoYs$Model <- factor(fut_DoYs$Model, levels=model.names)
fut_DoYs <- fut_DoYs %>%
  group_by(timeseries,Species,Model) %>%
  summarise(fut_mean=mean(fut_DoYoff,na.rm=T))
fut_GSLs <- fut_GSLs %>%
  pivot_longer(-c(timeseries,Species),names_to="Model",values_to="fut_GSL")
fut_GSLs$Model <- factor(fut_GSLs$Model, levels=model.names)
fut_GSLs <- fut_GSLs %>%
  group_by(timeseries,Species,Model) %>%
  summarise(fut_mean=mean(fut_GSL,na.rm=T))

# Select predictions of past senescence data (1990-2010 range) for the 6 models
past_DoYs <- fread("ModelAnalysis_1_Predicted_DoYoff.csv")
past_DoYs <- past_DoYs %>%
  filter(YEAR %in% 1990:2010) %>%
  select(timeseries,Pred_DoYoff_CDD,Pred_DoYoff_TPM,
         Pred_DoYoff_SIAM,Pred_DoYoff_TPDM,
         Pred_DoYoff_PIAgsi,`Pred_DoYoff_PIA+`) %>% 
  arrange(-desc(timeseries))
colnames(past_DoYs) <- c("timeseries",model.names)

# Calculate GSL
past_spring <- fread("DataMeta_2_PhenologyObs_PEP725_CleanData.csv")
past_spring <- past_spring %>% 
  filter(phenology=="leaf.out") %>% 
  filter(YEAR %in% 1990:2010) %>%
  filter(timeseries %in% unique(past_DoYs$timeseries)) %>% 
  select(timeseries,DoY) %>% 
  arrange(-desc(timeseries))
past_GSLs <- past_DoYs
past_GSLs[,-1] <- past_GSLs[,-1] - past_spring$DoY

# Calculate averages across-timeseries
past_DoYs <- past_DoYs %>%
  pivot_longer(-timeseries,names_to="Model",values_to="past_DoYoff")
past_DoYs$Model <- factor(past_DoYs$Model,levels=model.names)
past_DoYs <- past_DoYs %>%
  group_by(timeseries,Model) %>%
  summarise(past_mean=mean(past_DoYoff,na.rm=T))
past_GSLs <- past_GSLs %>%
  pivot_longer(-timeseries,names_to="Model",values_to="past_GSL")
past_GSLs$Model <- factor(past_GSLs$Model,levels=model.names)
past_GSLs <- past_GSLs %>%
  group_by(timeseries,Model) %>%
  summarise(past_mean=mean(past_GSL,na.rm=T))

# Join future and past datasets and calculate senescence delays/growing season length change
futpast_DoYs <- fut_DoYs %>%
  inner_join(past_DoYs,by=c("timeseries","Model")) %>%
  mutate(delay=fut_mean-past_mean)
futpast_DoYs <- futpast_DoYs[complete.cases(futpast_DoYs),]
futpast_GSLs <- fut_GSLs %>%
  inner_join(past_GSLs,by=c("timeseries","Model")) %>%
  mutate(change=fut_mean-past_mean)
futpast_GSLs <- futpast_GSLs[complete.cases(futpast_GSLs),]

# Calculate delay/change across timeseries
delay_doyoff <- futpast_DoYs %>%
  group_by(Model) %>%
  summarise(delay_ts=mean(delay),
            delay_se=se(delay))
change_GSL <- futpast_GSLs %>%
  group_by(Model) %>%
  summarise(change_ts=mean(change),
            change_se=se(change))

# Plot leaf senescence delay
fig_4b <- ggplot(data=delay_doyoff, aes(x=Model, y=delay_ts)) +
  geom_bar(position=position_dodge(), stat="identity", width=.9,
           fill=paletteSoftRainbow) +
  geom_errorbar(aes(ymin=delay_ts-delay_se, ymax=delay_ts+delay_se),
                width=.2,
                position=position_dodge(.9)) +
  ylab("Leaf senescence delay (days)") +
  coord_cartesian(ylim=c(-8,22)) +
  scale_y_continuous(position = 'right', expand=c(0,0)) +
  theme(aspect.ratio = 1,
        legend.position = "none",
        legend.background = element_rect(fill=NA,
                                         size=0.5, linetype = "solid"),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(linetype = "solid",fill=NA,colour = 'black',size = .3),
        axis.ticks.y=element_line(size=.3),
        axis.ticks.x=element_blank(),
        axis.line.x=element_line(size=.3),
        axis.line.y=element_line(size=.3),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=6),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8,vjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black')
  ) +
  geom_text(aes(x=6,y=19.5,label="b"),size=4,color="black",fontface="bold")
fig_4b

# Plot growing season length change
fig_4d <- ggplot(data=change_GSL, aes(x=Model, y=change_ts)) +
  geom_bar(position=position_dodge(), stat="identity", width=.9,
           fill=paletteSoftRainbow) +
  geom_errorbar(aes(ymin=change_ts-change_se, ymax=change_ts+change_se),
                width=.2,
                position=position_dodge(.9)) +
  ylab("Growing season length change (days)") +
  coord_cartesian(ylim=c(0,37)) +
  scale_y_continuous(position = 'right', expand=c(0,0)) +
  theme(aspect.ratio = 1,
        legend.position = "none",
        legend.background = element_rect(fill=NA,
                                         size=0.5, linetype = "solid"),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(linetype = "solid",fill=NA,colour = 'black',size = .3),
        axis.ticks.y=element_line(size=.3),
        axis.ticks.x=element_line(size=.3),
        axis.line.x=element_line(size=.3),
        axis.line.y=element_line(size=.3),
        axis.text.x=element_text(size=6, angle=45, hjust=1),
        axis.text.y=element_text(size=6),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8,vjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black')
  ) +
  geom_text(aes(x=6,y=34,label="d"),size=4,color="black",fontface="bold")
fig_4d


## SUPPLEMENTARY FIGURE 7 (Fig. S7)
# Future projections of autumn senescence dates (a) and growing season length (b) for six Central European species

# Calculate average senescence date/growing season length per Species
future_doyoff.sp <- doyoff.df %>% 
  mutate(Species=future.df$Species) %>% 
  mutate(YEAR=future.df$YEAR)
colnames(future_doyoff.sp) <- c(model.names,"Species","YEAR")
future_doyoff.sp <- future_doyoff.sp %>% 
  pivot_longer(-c(Species,YEAR),names_to="Model",values_to="FutProj_DoYoff") %>% 
  group_by(Species,Model,YEAR) %>% 
  summarise(mean_DoY=mean(FutProj_DoYoff,na.rm=T))
future_GSL.sp <- GSL.df %>% 
  mutate(Species=future.df$Species) %>% 
  mutate(YEAR=future.df$YEAR)
colnames(future_GSL.sp) <- c(model.names,"Species","YEAR")
future_GSL.sp <- future_GSL.sp %>% 
  pivot_longer(-c(Species,YEAR),names_to="Model",values_to="FutProj_GSL") %>% 
  group_by(Species,Model,YEAR) %>% 
  summarise(mean_GSL=mean(FutProj_GSL,na.rm=T))

# Calculate senescence projections/growing season lengths as 15-year moving averages 
window_width <- 15
future_doyoff.sp <- future_doyoff.sp %>% 
  group_by(Species,Model) %>% 
  mutate(roll_DoY=runmean(mean_DoY,window_width),
         sd_DoY=runsd(mean_DoY,window_width))
future_doyoff.sp$Model <- factor(future_doyoff.sp$Model,levels=model.names)
future_GSL.sp <- future_GSL.sp %>% 
  group_by(Species,Model) %>% 
  mutate(roll_GSL=runmean(mean_GSL,window_width),
         sd_GSL=runsd(mean_GSL,window_width))
future_GSL.sp$Model <- factor(future_GSL.sp$Model,levels=model.names)

# Define palette
paletteSoftRainbow <- c("#ABDAF4","#2b5a74",
                        "#e6d900","#e66600",
                        "#C7E9C0","#006D2C")

fig_s7a <- ggplot(future_doyoff.sp, aes(x=YEAR, y=round(roll_DoY), colour=Species))+
  labs(x = "Year", y = "Leaf senescence (DoY)") +
  geom_line(size=.5) +
  coord_cartesian(xlim=c(2029,2096.5), ylim=c(245,315))+
  theme(aspect.ratio = 1, 
        legend.position = "top",
        legend.background = element_rect(fill=NA, 
                                         size=0.5, linetype = "solid"),
        legend.title=element_blank(),
        legend.text=element_text(size=6),
        legend.key.size=unit(0.3,"cm"),
        legend.spacing.x = unit(0.1, 'cm'),
        strip.text=element_text(size=6),
        strip.background=element_rect(size=0.35),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(linetype = "solid",fill=NA,colour = 'black',size = .3),
        axis.ticks.y=element_line(size=.3),
        axis.ticks.x=element_line(size=.3),
        axis.line.x=element_line(size=.3),
        axis.line.y=element_line(size=.3),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.title.x=element_text(size=8,vjust=1),
        axis.title.y=element_text(size=8,vjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black')
  ) + facet_wrap(~Model, ncol=3,nrow=2)
fig_s7a

fig_s7b <- ggplot(future_GSL.sp, aes(x=YEAR, y=round(roll_GSL), colour=Species))+
  labs(x = "Year", y = "Growing season length (days)") +
  geom_line(size=.5) +
  scale_y_continuous(position = 'right', expand=c(0,0)) +
  coord_cartesian(xlim=c(2029,2096.5), ylim=c(130,230))+
  theme(aspect.ratio = 1, 
        legend.position = "top",
        legend.background = element_rect(fill=NA, 
                                         size=0.5, linetype = "solid"),
        legend.title=element_blank(),
        legend.text=element_text(size=6),
        legend.key.size=unit(0.3,"cm"),
        legend.spacing.x = unit(0.1, 'cm'),
        strip.text=element_text(size=6),
        strip.background=element_rect(size=0.35),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(linetype = "solid",fill=NA,colour = 'black',size = .3),
        axis.ticks.y=element_line(size=.3),
        axis.ticks.x=element_line(size=.3),
        axis.line.x=element_line(size=.3),
        axis.line.y=element_line(size=.3),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.title.x=element_text(size=8,vjust=1),
        axis.title.y=element_text(size=8,vjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black')
  ) + facet_wrap(~Model, ncol=3,nrow=2)
fig_s7b


## SUPPLEMENTARY FIGURE 8 (Fig. S8)
# Future delays in autumn senescence dates (a) and increases in growing season length (b) for six Central European species

# Calculate senescence delay/growing season length change per species
delays_sp <- futpast_DoYs %>% 
  group_by(Species,Model) %>% 
  summarise(delay_sp=mean(delay),
            delay_se=se(delay))
change_sp <- futpast_GSLs %>% 
  group_by(Species,Model) %>% 
  summarise(change_sp=mean(change),
            change_se=se(change))

# Plot senescence delay per species
fig_s8a <- ggplot(data=delays_sp, aes(x=Species, y=delay_sp, fill=Model)) +
  geom_bar(position = position_dodge(preserve = "single"),
           stat="identity") + 
  geom_errorbar(aes(ymin=delay_sp-delay_se, ymax=delay_sp+delay_se),
                size=.3,    
                width=.5,
                position = position_dodge(width = 0.9, preserve = "single")) +
  ylab("Leaf senescence delay (days)") +
  scale_fill_manual(values=paletteSoftRainbow) +
  geom_hline(aes(yintercept=0), size=.3) +
  theme(aspect.ratio = 1, 
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(linetype = "solid",fill=NA,colour = 'black',size = .3),
        axis.ticks.y=element_line(size=.3),
        axis.ticks.x=element_line(size=.3),
        axis.line.x=element_line(size=.3),
        axis.line.y=element_line(size=.3),
        axis.text.x=element_text(size=6, angle=45, hjust=1),
        axis.text.y=element_text(size=6),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8,vjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black')
  )
fig_s8a

# Plot senescence delay per species
fig_s8b <- ggplot(data=change_sp, aes(x=Species, y=change_sp, fill=Model)) +
  geom_bar(position = position_dodge(preserve = "single"),
           stat="identity") + 
  geom_errorbar(aes(ymin=change_sp-change_se, ymax=change_sp+change_se),
                size=.3,    
                width=.5,
                position = position_dodge(width = 0.9, preserve = "single")) +
  ylab("Growing season length change (days)") +
  scale_y_continuous(position = 'right', expand=c(0,0)) +
  scale_fill_manual(values=paletteSoftRainbow) +
  coord_cartesian(ylim=c(-15,45)) +
  geom_hline(aes(yintercept=0), size=.3) +
  theme(aspect.ratio = 1, 
        legend.title=element_blank(),
        legend.text=element_text(size=6),
        legend.key.size=unit(0.3,"cm"),
        legend.background = element_rect(fill=NA, 
                                         size=0.5, linetype = "solid"),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(linetype = "solid",fill=NA,colour = 'black',size = .3),
        axis.ticks.y=element_line(size=.3),
        axis.ticks.x=element_line(size=.3),
        axis.line.x=element_line(size=.3),
        axis.line.y=element_line(size=.3),
        axis.text.x=element_text(size=6, angle=45, hjust=1),
        axis.text.y=element_text(size=6),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8,vjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black')
  )
fig_s8b


## SUPPLEMENTARY FIGURE 9 (Fig. S9)
## Variation in annual net photosynthesis in our dataset

## Sub-panel A
# Plot frequency distribution showing the total variation across years and sites
drivers.df <- fread("DataMeta_3_Drivers.csv")
dens_to_plot.df  <- drivers.df %>% 
  select(YEAR,cA_tot,`cA_tot-w`) %>%
  rename(`Non-water stressed`=cA_tot,`Water-stressed`=`cA_tot-w`) %>% 
  pivot_longer(-YEAR,names_to="Photosynthesis")

# Plot
fig_S9a <- ggplot(data=dens_to_plot.df, aes(x=value)) + 
  geom_histogram(aes(color = Photosynthesis, fill = Photosynthesis),
                 alpha = 0.4, position = "identity") +
  xlab(bquote('Net photosynthetic rate (gC '~m^-2~year^-1*')')) +
  scale_color_manual(values = c("#F8766D", "#619CFF"))+
  scale_fill_manual(values = c("#F8766D", "#619CFF")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(aspect.ratio = 1, 
        legend.position=c(0.25,0.89),
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
fig_S9a


## Sub-panel B
# Plot Temporal variation in annual photosynthesis

# Calculate 15-year moving averages
window_width <- 15

# Water-stressed
water.df <- drivers.df %>% 
  select(YEAR,`cA_tot-w`) %>%
  rename(`Water-stressed`=`cA_tot-w`) %>% 
  group_by(YEAR) %>% 
  summarise(value=mean(`Water-stressed`)) %>% 
  mutate(roll_photo=caTools::runmean(value,window_width),
         sd_photo=caTools::runsd(value,window_width),
         Photosynthesis="Water-stressed")

# Non-water-stressed
nwater.df <- drivers.df %>% 
  select(YEAR,cA_tot) %>%
  rename(`Non-water stressed`=cA_tot) %>% 
  group_by(YEAR) %>% 
  summarise(value=mean(`Non-water stressed`)) %>% 
  mutate(roll_photo=caTools::runmean(value,window_width),
         sd_photo=caTools::runsd(value,window_width),
         Photosynthesis="Non-water-stressed")
phototrend_to_plot.df <- rbind(water.df,nwater.df)
colnames(phototrend_to_plot.df)

# Plot
fig_9b <- ggplot(phototrend_to_plot.df, aes(x=YEAR, y=roll_photo, color=Photosynthesis))+
  ylab(bquote('Net photosynthetic rate (gC '~m^-2~year^-1*')')) +
  xlab("Year") +
  geom_line(size=0.75) +
  geom_ribbon(aes(ymin=roll_photo-sd_photo, ymax=roll_photo+sd_photo, group=Photosynthesis, fill=Photosynthesis), linetype=0, alpha=0.07) +
  scale_color_manual(values=c("#F8766D", "#619CFF")) +
  scale_fill_manual(values=c("#F8766D", "#619CFF")) +
  coord_cartesian(xlim=c(min(phototrend_to_plot.df$YEAR),max(phototrend_to_plot.df$YEAR))) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(aspect.ratio = 1/2.5, 
        legend.position=c(0.85,0.175),
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
fig_S9b