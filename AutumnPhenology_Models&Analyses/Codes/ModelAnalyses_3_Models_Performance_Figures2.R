####################
## Model Performance

# Contains:
# > Calculation of performance indexes (R2, RMSE, slope)
# > 5-fold site-specific validation
# > Code for Figures 2A-2B-2C-2D

# Define directory path
setwd(".../AutumnPhenology/AutumnPhenology_Models&Analyses/Data/")

# Load libraries
library(data.table)
library(tidyverse)
library(broom)
library(Metrics)


##----------------------------------------
# Calculation of Performance Indexes

# From the text:
# We evaluated the quality of the model predictions with three performance indexes 
# of observed (y-axis) vs. predicted (x-axis) values of autumn anomaly: 
# 1) the overall fit (R2 values);
# 2) the error fit (root mean square error, RMSE); 
# 3) the accuracy of the fit (slope values). 
# Additionally, the site-specific performance of the models was evaluated using 
# 5-fold cross-validation over the observation period at the time-series scale ("site-specific validation").

# Import data
pred_DoYoff <- fread("ModelAnalysis_1_Predicted_DoYoff.csv")

# Define model names
model.names   <- c("CDD","DM1","DM2","TPM",
                   "SIAM","TDM","TPDM",
                   "GSIAM","PIAM-W","PIAM")

# Format data for calculation 
pred_DoYoff <- pred_DoYoff %>% 
  select(timeseries,Obs_DoYoff,
         Pred_DoYoff_CDD,Pred_DoYoff_DM1,Pred_DoYoff_DM2,Pred_DoYoff_TPM,
         Pred_DoYoff_SIAM,Pred_DoYoff_TDM,Pred_DoYoff_TPDM,
         Pred_DoYoff_GSIAM,`Pred_DoYoff_PIAM-W`,Pred_DoYoff_PIAM)
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
## FIGURE 2A - 2B - 2C - 2D
## Model Performance

## MODELS OF LEAF SENESCENCE (and drivers):

## First-generation:
# CDD (chilling temperature) - Dufr�ne et al. (2005)
# DM1 and DM2 (chilling temperature, autumn daylength) - Delpierre et al. (2009)
# TPM (chilling temperature, autumn daylength) - Lang et al. (2019)
## Second-generation:
# SIAM (chilling temperature, autumn daylength, spring anomaly) - Keenan and Richardson (2015)
# TDM and TPDM (chilling temperature, autumn daylength, growing season temperature / + water stress) - Liu et al. (2019)
## CarbLim:
# GSIAM (chilling temperature, autumn daylength, leaf flushing date, growing season mean temperature, daylength, vapour pressure deficit)
# PIAM (chilling temperature, autumn daylength, leaf flushing date, growing season mean temperature, daylength, precipitation, net radiation, CO2 concentrationwater stress)

# Load libraries
library(data.table)
library(tidyverse)
library(lmodel2)
library(ggplot2)
library(Metrics)


################################
## FIGURE 2A
## Observed (Obs_AnomDoYoff) vs.
## Predicted (Pred_AnomDoYoff)
## leaf senescence anomalies, i.e., as deviation from the mean observed leaf-out date at each site
## of the best-performing first-generation (TPM) [Lang et al. (2019)], 
## second-generation (TPDM) [Liu et al. (2019)]
## and CarbLim models (PIAM)

# Import data
pred_DoYoff <- fread("ModelAnalysis_1_Predicted_DoYoff.csv")

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
fig_2a <- ggplot(best_pred_DoYoff, aes(x=predictions, y=observations)) +
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
fig_2a

# Add R^2 values
# Average across time-series
# See Figure 2B
dat_text <- data.frame(
  label = c("R2 = 0.21", "R2 = 0.62", "R2 = 0.78"), #check R2_models (Figure 2B)
  model_type = types
)
fig_2a  <- fig_2a + geom_text(
  data = dat_text,
  mapping = aes(x=-35, y=45, label=label)
)

# Calculate intercept plus slope by Standard Major Axis (SMA)
SMA_values <- best_pred_DoYoff %>%  
  group_by(model_type) %>% 
  summarise(int=lmodel2(observations ~ predictions)$regression.results$Intercept[3],
            slope=lmodel2(observations ~ predictions)$regression.results$Slope[3])
fig_2a <- fig_2a +
  geom_abline(data = SMA_values,
              mapping = aes(intercept=int,slope=slope),
              linetype = "solid")
fig_2a


################################
## FIGURE B
## R^2 values

# Import data
stat_models <- fread("Models_Performance_R2_RMSE_slope_v2.csv")

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
            R2_sd = sd(value)) %>% 
  mutate(Type=types)
R2_models$Type <- factor(R2_models$Type, levels=c("First-generation","Second-generation","PIA models"))

# Define palette
trio <- c("#56B4E9","#E69F00","#74C476")

# Plot
fig_2b <- ggplot(R2_models, aes(x=Model, y=R2_mean)) +
  geom_bar(position=position_dodge(), stat="identity", width=.9,
           fill=c("#56B4E9","#56B4E9","#56B4E9","#56B4E9",
                  "#E69F00","#E69F00","#E69F00",
                  "#74C476","#74C476")) +
  geom_errorbar(aes(ymin=R2_mean-R2_sd, ymax=R2_mean+R2_sd),
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
fig_2b


################################
## FIGURE C
## RMSE values
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
RMSE_null <- rmse(rep(mean(pred_DoYoff$Obs_DoYoff),nrow(pred_DoYoff)), pred_DoYoff$Obs_DoYoff)

# Define palette
trio <- c("#56B4E9","#E69F00","#74C476")

# Plot
fig_2c <- ggplot(RMSE_models, aes(x=Model, y=xval_median, color=Type)) +
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
fig_2c


################################
## FIGURE D
## slope values

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
fig_2d <- ggplot(slope_models, aes(x=Model, y=value, color=Type)) +
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
fig_2d