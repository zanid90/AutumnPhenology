#---
#title: Increased growing-season productivity drives earlier autumn leaf senescence in temperate trees (Zani et al. 2020 Science)
#author: Constantin Zohner, Deborah Zani
#date: "last updated May 24, 2022"

#R code creating the process-based model figures



## Model names

### PIA models
#- PIA+: Photosynthesis-influenced autumn phenology model (LPJ photosynthesis)
#- PIA-: Photosynthesis-influenced autumn phenology model (LPJ photosynthesis without water stress)
#- PIAgsi: Photosynthesis-influenced autumn phenology model (growing-season index photosynthesis)

### Second-generation models
#- SIAM: Spring-influenced autumn phenology model (Keenan & Richardson 2014)
#- TDM: Temperature-influenced Delpierre model (Liu et al. 2018)
#- TPDM: Temperature and precipitation-influenced Delpierre model (Liu et al. 2018)

### First-generation models
#- TPM: Low Temperature And Photoperiod Multiplicative model (Lang et al. 2019)
#- DM1: Delpierre model 1 (Delpierre et al. 2009)
#- DM2: Delpierre model 2 (Delpierre et al. 2009)
#- CDD: Cold-Degree-Day model (Dufrêne et al. 2005)



#####################
# Required packages #
#####################



require(tidyverse)
require(data.table)
require(ggplot2)
require(wesanderson)
require(patchwork)
require(broom)



##############################################################################################################################################
##############################################################################################################################################



################
## Plot theme ##
################



plotTheme1 = theme(
  legend.position   = "right",
  legend.background = element_blank(),
  legend.text       = element_text(color="black"),
  legend.title      = element_blank(),
  legend.key        = element_blank(),
  panel.grid.major  = element_blank(),
  panel.grid.minor  = element_blank(),
  panel.background  = element_blank(),
  panel.border      = element_rect(colour = "black", fill=NA),
  axis.line         = element_line(color = "black"),
  axis.text         = element_text(colour = "black"),
  strip.background  = element_rect(fill=NA),
  strip.text        = element_text(colour = 'black',face = "italic"),
  plot.title        = element_text(face="bold"))



##############################################################################################################################################
##############################################################################################################################################



###############
## Functions ##
###############



##################
# Model efficiency
##################

# R2 values represent the coefficient of determination relative to the 1:1 line of observed versus predicted values. 
# This is equivalent to a standardized mean squared error (see e.g. Ma et al. 2021; https://doi.org/10.1038/s41559-021-01485-1)
# negative R2 values indicate that the model performs worse than the mean
R2 <- function(xtrue, xpred){
  return(1-sum((xtrue-xpred)^2)/sum((xtrue-mean(xtrue))^2))
}   


#############
# Adjusted R2
#############

# the percentage of variation explained by only the independent variables that actually affect the dependent variable
# penalizes you for adding independent variables/parameters (k in the equation) that do not fit the model. 
R2.adj = function (R2,n,k) {
  R2.adj = 1 - ( (1-R2)*(n-1) / (n-k-1) )
  return(R2.adj)
}
#R2...R2 in relation to 1:1 line, see R2 function above
#n...number of points in data sample
#k...number of independent regressors, i.e. the number of variables in your model, excluding the constant


########################
# Root mean square error
########################

RMSE <- function(xtrue, xpred){
  return(sqrt(mean((xtrue-xpred)^2)))
}



##############################################################################################################################################
##############################################################################################################################################



#####################
## Set directories ##
#####################



#input
PhenoR_predictions_path = "/Users/consti/Desktop/PhD/Publication_material/000_Zani_et_al_reanalysis/2_Combined_analysis/Analysis_output/Process_models/Merged_files"

# Output
output_path             = "/Users/consti/Desktop/PhD/Publication_material/000_Zani_et_al_reanalysis/2_Combined_analysis/Analysis_output/Process_models/Figures"



##############################################################################################################################################
##############################################################################################################################################



#################
## Import data ##
#################



Predictions.df = fread(paste(PhenoR_predictions_path,"Process_model_predictions.csv",sep="/")) %>%
  dplyr::select(-V1) %>% 
  #add timeseries count
  group_by(timeseries)%>%
  add_tally()%>%
  ungroup() %>%
  #long format
  pivot_longer(., -c(timeseries, pep_id, lon, lat, species, year, leaf_off, n), 
               names_to = "model", values_to = "leaf_off_pred")%>%
  #create model names
  mutate(model = substring(model, 14),
         model_type = case_when(model %in% c("CDD","DM1","DM2","TPM") ~ "1st generation",
                                model %in% c("TPDM","TDM","SIAM") ~ "2nd generation",
                                model %in% c("PIAgsi","PIAminus","PIA") ~ "PIA model"))



##############################################################################################################################################
##############################################################################################################################################



#######################
## Model performance ##
#######################



################
# Get model info
################

ModelResults.df = Predictions.df %>%
  
  group_by(timeseries, model_type, model, n)%>%
  
  do({
    
    #run models
    ###########
    
    model.Obs.Pred = lm(leaf_off~leaf_off_pred, data=.)
    
    
    #create combined dataframe
    ##########################
    
    data.frame(
      RMSE = RMSE(.$leaf_off, .$leaf_off_pred),
      R2 = R2(.$leaf_off, .$leaf_off_pred),
      tidy(model.Obs.Pred),
      glance(model.Obs.Pred) )
  })%>%
  filter(!term=="(Intercept)")%>%
  ungroup() %>%
  
  #Add adjusted R2 based on number of free parameters  
  mutate(n.params = case_when(model=="CDD" ~ 2,
                              model %in% c("DM1","DM2") ~ 3,
                              model %in% c("TPM") ~ 4,
                              model %in% c("PIA","PIAminus","PIAgsi","SIAM","TDM") ~ 6,
                              model %in% c("TPDM") ~ 7),
         R2.adj = R2.adj(R2, nobs, n.params))


#############################
# Summarize results per model 
#############################

summary.df = ModelResults.df %>%
  group_by(model_type, model) %>%
  summarize(mean.R2adj  = mean(R2.adj), 
            lowCI.R2adj = t.test(R2.adj)$conf.int[1],
            hiCI.R2adj  = t.test(R2.adj)$conf.int[2],
            
            mean.R2    = mean(R2), 
            lowCI.R2   = t.test(R2)$conf.int[1],
            hiCI.R2    = t.test(R2)$conf.int[2],
            
            mean.rmse  = mean(RMSE), 
            lowCI.rmse = t.test(RMSE)$conf.int[1],
            hiCI.rmse  = t.test(RMSE)$conf.int[2],
            
            mean.slope  = mean(estimate, na.rm=T), 
            lowCI.slope = t.test(estimate)$conf.int[1],
            hiCI.slope  = t.test(estimate)$conf.int[2]
  ) %>%
  
  arrange(mean.R2, model) %>%
  mutate(model = factor(model, levels=with(., reorder(model, mean.R2, na.rm=T)), ordered=T))
         
summary.df


#######
# Plots
#######

# R2
summary.df$model
ModelLabs <- c("CDD", "DM2", "DM1", "TPM", "TPDM", "TDM", "SIAM",
               expression("PIA"["GSI"]), expression("PIA"^"-"), expression("PIA"^"+"))
               
R2.barplot = ggplot(data=summary.df, aes(x = model, y = mean.R2, fill=model_type)) + 
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = lowCI.R2, ymax = hiCI.R2), width=.2,
                position=position_dodge(.9)) +
  geom_hline(yintercept=0)+
  scale_fill_manual(values = c("#5BB4E5", '#E6A024', '#76C375'))+
  xlab("") + ylab(expression("Coefficient of determination ( R"^2~")") ) +
  coord_cartesian(ylim=c(0.023,.5))+
  scale_x_discrete(labels= ModelLabs)+
  plotTheme1 +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

#save plots as .pdf
ggsave(R2.barplot, file="Process_models_R2.pdf", path=output_path,
       width=6.5, height=3)

R2.barplot





#######################
## Model performance ##
#######################



model.increaseR2 = data.frame(length=c(15:60), SecondGen=NA,FirstGen=NA)
model.R2         = data.frame(length=c(15:60), SecondGen=NA,FirstGen=NA,PIA=NA)

for(i in (1:46)) {

  # Summarize results per model 
  #############################
  
  summary.df = ModelResults.df %>%
    
    #delete timeseries with fewer than 15-65 years
    filter(n >= i+14) %>%
    
    group_by(model_type, model) %>%
    summarize(mean.R2adj  = mean(R2.adj), 
              lowCI.R2adj = t.test(R2.adj)$conf.int[1],
              hiCI.R2adj  = t.test(R2.adj)$conf.int[2],
              .groups = "keep"
    ) %>%
    #order table
    arrange(mean.R2adj, model) %>%
    mutate(model = factor(model, levels=with(., reorder(model, mean.R2adj, na.rm=T)), ordered=T))
  
  
  # Increase in model accuracy compared to TPDM and SIAM models
  #############################################################
  
  
  # extract maximum adjusted R2s for first and second gen models and minimum R2s for PIA models
  FirstGen.R2adj  = max((summary.df %>% filter(model_type %in% c('1st generation')))$mean.R2adj)
  SecondGen.R2adj = max((summary.df %>% filter(model_type %in% c('2nd generation')))$mean.R2adj)
  PIA.R2adj       = min((summary.df %>% filter(model_type %in% c('PIA model')))$mean.R2adj)
  
  # Safe relative model performance increase for adjusted R2
  model.increaseR2[i,]$FirstGen  = (PIA.R2adj/FirstGen.R2adj-1)*100
  model.increaseR2[i,]$SecondGen = (PIA.R2adj/SecondGen.R2adj-1)*100
  
  # Safe absolute adjusted R2s
  model.R2[i,]$FirstGen  = FirstGen.R2adj
  model.R2[i,]$SecondGen = SecondGen.R2adj
  model.R2[i,]$PIA       = PIA.R2adj
  
  # print performance increase
  print(paste0(i+14,' years: ',round((PIA.R2adj/max(c(SecondGen.R2adj,FirstGen.R2adj))-1)*100,1),'-', 
               round((PIA.R2adj/min(c(SecondGen.R2adj,FirstGen.R2adj))-1)*100,1),'% increase in performance (adjusted R2)'))
}



#################################################################################
# Compare adjusted R2 of worst PIA model with best First- and Second-gen models #
#################################################################################



# Plot increase in PIA model performance relative to first- and second-gen models
#################################################################################


#unfold dataset  
model.increaseR2 = model.increaseR2 %>%
  pivot_longer(., -c(length), names_to = "model", values_to = "increase")%>%
  mutate(model = factor(model, levels=c("FirstGen", "SecondGen")))
model.increaseR2

#Plot
performance.plot.R2 = ggplot(data=model.increaseR2, aes(x = length, y = increase, group=model, colour=model)) + 
  geom_hline(yintercept=0)+
  geom_hline(yintercept=25,  color='lightgrey')+
  geom_hline(yintercept=50,  color='lightgrey')+
  geom_hline(yintercept=75,  color='lightgrey')+
  geom_hline(yintercept=-25, color='lightgrey')+
  geom_hline(yintercept=-50, color='lightgrey')+
  geom_hline(yintercept=-75, color='lightgrey')+
  
  geom_point(stat = "identity")+
  geom_line() +
  
  scale_color_manual(
    values = c("#5BB4E5", '#E6A024'))+
  
  xlab("Minimum length of time series (years)") + ylab("Increase in PIA model performance (%)") +
  coord_cartesian(ylim=c(-92,92),xlim=c(16.7,58.3))+
  plotTheme1+
  theme(legend.text.align = 0)
performance.plot.R2



# Plot PIA, first gen, and second gen model performance
#######################################################


#unfold dataset  
model.R2 = model.R2 %>%
  pivot_longer(., -c(length), names_to = "model", values_to = "Adjusted_R2") %>%
  mutate(model = factor(model, levels=c("FirstGen", "SecondGen", "PIA")))
model.R2

#Plot
R2.plot = ggplot(data=model.R2, aes(x = length, y = Adjusted_R2, group=model, colour=model)) + 
  geom_hline(yintercept=.05, color='lightgrey')+
  geom_hline(yintercept=.1,  color='lightgrey')+
  geom_hline(yintercept=.15, color='lightgrey')+
  geom_hline(yintercept=.2,  color='lightgrey')+
  
  geom_point(stat = "identity")+
  geom_line() +
  
  scale_color_manual(
    values = c("#5BB4E5", '#E6A024', '#76C375'))+
  
  xlab("") + ylab(expression("Adjusted R"^2)) +
  coord_cartesian(ylim=c(0.0095,.13),xlim=c(16.7,58.3))+
  plotTheme1+
  theme(legend.text.align = 0)
R2.plot



###########################################################
# Number of time series in relation to time series length #
###########################################################


#add year count to table
Count.df = Predictions.df %>%
  group_by(timeseries) %>%
  filter(row_number()==1)

#total time series
nrow(Count.df)

#time series with at least 15 years of data
nrow(Count.df %>%
       filter(n>14))

#count the time series
model.count = data.frame(length=c(15:60),cum.length=NA)
for (i in 1:46) {
  count = nrow(Count.df %>% filter(n>=i+14))
  model.count[i,2] = count
}

#Plot
count.plot = ggplot(data=model.count, aes(x = length, y = cum.length)) + 
  geom_point(stat = "identity")+
  geom_line() +
  xlab("Minimum length of time series (years)") + ylab("Number of time series") +
  coord_cartesian(xlim=c(16.7,58.3),ylim=c(600,14000))+
  plotTheme1
count.plot

#define plot layout
layout <- "
A
B"

#Merge plots
Fig3_Plot_R2 = R2.plot + performance.plot.R2 + 
  plot_layout(design = layout) + plot_annotation(tag_levels = 'A')&
  theme(plot.tag = element_text(face = 'bold'))

#save plots as .pdf
ggsave(Fig3_Plot_R2, file="Process_models_years_R2.pdf", path=output_path,
       width=6.5, height=6.5)
Fig3_Plot_R2





#####################
## Reproducibility ##	
#####################



## date time
Sys.time()

## session info
sessionInfo()



##############################################################################################################################################
#############################################################THE END##########################################################################
##############################################################################################################################################