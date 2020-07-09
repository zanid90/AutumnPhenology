################################
## Correlation and Path Analyses

## From the text:
# See Materials and Methods
# Correlation and path analyses

# Define directory paths
setwd(".../AutumnPhenology/AutumnPhenology_Data&Metadata/Data/")

# Load libraries
library(data.table)
library(tidyverse)

# Define auxiliary functions
se <- function(x) sqrt(var(x)/length(x))

# Import data
drivers.df <- fread("DataMeta_3_Drivers.csv")


##----------------------------------------
## Explore correlations among variables
library(ggcorrplot)

# DoY_off = autumn leaf senescence dates vs.

## Spring phenology (Fu et al. 2014):
# DoY_out = spring leaf-out dates

## Environmental predictors of leaf senescence (Xie et al. 2018; Liu et al. 2019):
# temp_GS = average temperature during the growing season
# temp_aut2 = average minimum temperature during autumn (2 months before leaf senescence)
# temp_aut3 = average minimum temperature during autumn (3 months before leaf senescence; Fu et al. 2014)
# RD = annual rainy days (number of days with precipitation > 2 mm)
# RD_summer = rainy days during summer (water supply during the driest season)
# HRD = annual heavy rainy days (number of days with precipitation > 20 mm)
# HRD_aut = heavy rainy days during autumn
# HD35 = extreme heat event (number of days with maximum temperature > 35�C)
# FD = frost days (number of days with minimum temperature < 0�C)
# FD_spring = early frost events
# ci = intercellular CO2 concentration (atmospheric CO2 concentration multiplied by vapour pressure deficit function)

## Photosynthesis activity
# cGSI = cumulative/annuals Growing Season Index (function of temperature, day-length and vapour pressure deficit)
# cA_tot-w = annual net photosynthetic activity
# cA_tot = accounting for water supply during the growing season

# Compute a correlation matrix
corr <- drivers.df %>% 
  select(-c(timeseries,PEP_ID,LON,LAT,Species,YEAR,
             autumn_anomaly,spring_anomaly,temp_aut3,CO2)) %>% 
  cor() %>% # default method = "pearson"
  round(.,2)

# Visualize the correlation matrix
# SUPPLEMENTARY FIGURE 1 (Fig. S1)
fig_s1 <- ggcorrplot(corr,
                     type = "lower", # Get the lower triangle
                     colors = c("red3", "white", "blue3"), 
                     legend.title = "",
                     outline.col = "white",
                     lab = TRUE) # Add correlation coefficients
fig_s1
ggsave(filename = "FigureS1.jpeg",
       device = "jpeg",
       width = 5.8*4, height=5.8*4, units = "cm",
       dpi = 600)


##----------------------------------------
## Partial Correlations
library(ppcor)

# autumn_anomaly = autumn leaf senescence dates vs.
# NPP_anomaly = photosynthesis activity accounting for water deficit
# tempGS_anomaly = average temperature during the growing season (Xie et al. 2018; Liu et al. 2019)
# RD_summer = rainy days during summer (water supply during the driest season) (Xie et al. 2018; Liu et al. 2019)
# spring_anomaly = spring leaf-out dates (Fu et al. 2014)
# tempaut2_anomaly = average minimum temperature during autumn (2 months before leaf senescence; Fu et al. 2014)
# anomalies = 

# Select variables
varsel <- c("DoY_off", "cA_tot","DoY_out","temp_GS","temp_aut2","HRD","ci")

# Calculate site-anomalies
drivers_sel_anomalies.df <- drivers.df %>% 
  select(c(timeseries,Species,PEP_ID,YEAR,varsel)) %>% 
  group_by(PEP_ID) %>% 
  mutate(autumn_anomaly = mean(DoY_off),
         spring_anomaly = mean(DoY_out),
         tempGS_anomaly = mean(temp_GS),
         tempaut2_anomaly = mean(temp_aut2),
         prec_anomaly = mean(HRD),
         ci_anomaly = mean(ci),
         NPP_anomaly = mean(cA_tot)) %>% 
  ungroup(PEP_ID)
drivers_sel_anomalies.df$autumn_anomaly <- drivers_sel_anomalies.df$DoY_off - drivers_sel_anomalies.df$autumn_anomaly
drivers_sel_anomalies.df$spring_anomaly <- drivers_sel_anomalies.df$DoY_out - drivers_sel_anomalies.df$spring_anomaly
drivers_sel_anomalies.df$tempGS_anomaly <- drivers_sel_anomalies.df$temp_GS - drivers_sel_anomalies.df$tempGS_anomaly
drivers_sel_anomalies.df$tempaut2_anomaly <- drivers_sel_anomalies.df$temp_aut2 - drivers_sel_anomalies.df$tempaut2_anomaly
drivers_sel_anomalies.df$prec_anomaly <- drivers_sel_anomalies.df$HRD - drivers_sel_anomalies.df$prec_anomaly
drivers_sel_anomalies.df$ci_anomaly <- drivers_sel_anomalies.df$ci - drivers_sel_anomalies.df$ci_anomaly
drivers_sel_anomalies.df$NPP_anomaly <- drivers_sel_anomalies.df$cA_tot - drivers_sel_anomalies.df$NPP_anomaly


## Partial correlation - Time-series level
# Control for all variables
pcorr.tm <- drivers_sel_anomalies.df %>%
  group_by(timeseries) %>%
  summarise(pcorr_cAtot = pcor.test(autumn_anomaly, NPP_anomaly, c(ci_anomaly,tempGS_anomaly,prec_anomaly,tempaut2_anomaly,spring_anomaly))[,1],
            pcorr_CO2 = pcor.test(autumn_anomaly, ci_anomaly, c(NPP_anomaly,tempGS_anomaly,prec_anomaly,tempaut2_anomaly,spring_anomaly))[,1],
            pcorr_tempGS = pcor.test(autumn_anomaly, tempGS_anomaly, c(ci_anomaly,NPP_anomaly,prec_anomaly,tempaut2_anomaly,spring_anomaly))[,1],
            pcorr_prec = pcor.test(autumn_anomaly, prec_anomaly, c(ci_anomaly,tempGS_anomaly,NPP_anomaly,tempaut2_anomaly,spring_anomaly))[,1],
            pcorr_tempaut2 = pcor.test(autumn_anomaly, tempaut2_anomaly, c(ci_anomaly,tempGS_anomaly,prec_anomaly,NPP_anomaly,spring_anomaly))[,1],
            pcorr_spring = pcor.test(autumn_anomaly, spring_anomaly, c(ci_anomaly,tempGS_anomaly,prec_anomaly,tempaut2_anomaly,NPP_anomaly))[,1]
  )

# Average values across timeseries
meanpcorr_across_timeseries <- pcorr.tm %>%
  summarise(`NPP anomaly` = mean(pcorr_cAtot),
            `CO2 anomaly` = mean(pcorr_CO2),
            `Summer temperature anomaly` = mean(pcorr_tempGS),
            `Precipitation anomaly` = mean(pcorr_prec),
            `Spring anomaly` = mean(pcorr_spring),
            `Autumn temperature anomaly` = mean(pcorr_tempaut2)) %>%
  add_column(timeseries = "Across timeseries") %>%
  pivot_longer(cols=-timeseries,names_to="predictors",values_to="mean_values")

# Standard errors across timeseries
SE_across_timeseries <- pcorr.tm %>%
  summarise(`NPP anomaly` = se(pcorr_cAtot),
            `CO2 anomaly` = se(pcorr_CO2),
            `Summer temperature anomaly` = se(pcorr_tempGS),
            `Precipitation anomaly` = se(pcorr_prec),
            `Spring anomaly` = se(pcorr_spring),
            `Autumn temperature anomaly` = se(pcorr_tempaut2)) %>%
  pivot_longer(cols=1:6,names_to="predictors",values_to="SE_values")
meanpcorr_across_timeseries <- meanpcorr_across_timeseries %>%
  left_join(SE_across_timeseries,by="predictors")

# Set the order of predictors to plot
meanpcorr_across_timeseries$predictors <- factor(meanpcorr_across_timeseries$predictors,
                                                 levels=c("NPP anomaly",
                                                          "CO2 anomaly",
                                                          "Summer temperature anomaly",
                                                          "Precipitation anomaly",
                                                          "Spring anomaly",
                                                          "Autumn temperature anomaly"
                                                 ))

# Plot
# FIGURE 1B
# Define palettes
paletteBlueRed <- c("blue3","white","red3")
plot_colors <-  colorRampPalette(paletteBlueRed)(6)

fig_1b <- ggplot(data=meanpcorr_across_timeseries, aes(x=predictors, y=mean_values)) +
  geom_bar(position=position_dodge(), stat="identity", width=.9,
           fill=plot_colors) +
  geom_errorbar(aes(ymin=mean_values-SE_values, ymax=mean_values+SE_values),
                size=.3,
                width=.5,
                position=position_dodge(.9)) +
  labs(y = "Partial correlation coefficient") +
  coord_cartesian(ylim=c(-0.60,0.30)) +
  scale_y_continuous(breaks = c(-0.45,-0.30,-0.15,0,0.15,0.30), position="right") +
  geom_hline(aes(yintercept=0), size=.3) +
  theme(aspect.ratio = 1,
        axis.title.y=element_text(size=11, vjust=1),
        axis.text.y=element_text(size=7),
        axis.ticks.y=element_line(size=.3),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, size=7),
        axis.ticks.x=element_line(size=.3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black')
  )
fig_1b


## Partial correlation - Species level
pcorr.sp <- drivers_sel_anomalies.df %>%
  group_by(Species,timeseries) %>%
  summarise(pcorr_cAtot = pcor.test(autumn_anomaly, NPP_anomaly, c(ci_anomaly,tempGS_anomaly,prec_anomaly,tempaut2_anomaly,spring_anomaly))[,1],
            pcorr_CO2 = pcor.test(autumn_anomaly, ci_anomaly, c(NPP_anomaly,tempGS_anomaly,prec_anomaly,tempaut2_anomaly,spring_anomaly))[,1],
            pcorr_tempGS = pcor.test(autumn_anomaly, tempGS_anomaly, c(ci_anomaly,NPP_anomaly,prec_anomaly,tempaut2_anomaly,spring_anomaly))[,1],
            pcorr_prec = pcor.test(autumn_anomaly, prec_anomaly, c(ci_anomaly,tempGS_anomaly,NPP_anomaly,tempaut2_anomaly,spring_anomaly))[,1],
            pcorr_tempaut2 = pcor.test(autumn_anomaly, tempaut2_anomaly, c(ci_anomaly,tempGS_anomaly,prec_anomaly,NPP_anomaly,spring_anomaly))[,1],
            pcorr_spring = pcor.test(autumn_anomaly, spring_anomaly, c(ci_anomaly,tempGS_anomaly,prec_anomaly,tempaut2_anomaly,NPP_anomaly))[,1]
  )

# Average values per species
meanpcorr_per_species <- pcorr.sp %>%
  summarise(`NPP anomaly` = mean(pcorr_cAtot),
            `CO2 anomaly` = mean(pcorr_CO2),
            `Summer temperature anomaly` = mean(pcorr_tempGS),
            `Precipitation anomaly` = mean(pcorr_prec),
            `Spring anomaly` = mean(pcorr_spring),
            `Autumn temperature anomaly` = mean(pcorr_tempaut2)) %>%
  pivot_longer(cols=-Species,names_to="predictors",values_to="mean_values")

# Standard errors per species
SE_per_species <- pcorr.sp %>%
  summarise(`NPP anomaly` = se(pcorr_cAtot),
            `CO2 anomaly` = se(pcorr_CO2),
            `Summer temperature anomaly` = se(pcorr_tempGS),
            `Precipitation anomaly` = se(pcorr_prec),
            `Spring anomaly` = se(pcorr_spring),
            `Autumn temperature anomaly` = se(pcorr_tempaut2)) %>%
  pivot_longer(cols=2:7,names_to="predictors",values_to="SE_values")
meanpcorr_per_species$SE_values <- SE_per_species$SE_values

# Set the order of predictors to plot
meanpcorr_per_species$predictors <- factor(meanpcorr_per_species$predictors,
                                              levels=c("NPP anomaly",
                                                       "CO2 anomaly",
                                                       "Summer temperature anomaly",
                                                       "Precipitation anomaly",
                                                       "Spring anomaly",
                                                       "Autumn temperature anomaly"
                                              ))

# Define palette
paletteBlueRed <- c("blue3","white","red3")
plot_colors <-  colorRampPalette(paletteBlueRed)(6)

# Plot 
# SUPPLEMENTARY FIGURE 2 (Fig. S2)
sup_fig_5 <- ggplot(data = meanpcorr_per_species, aes(x=Species, y=mean_values, fill=predictors)) +
  geom_bar(position = position_dodge(preserve = "single"),
            stat="identity") + 
  geom_errorbar(aes(ymin=mean_values-SE_values, ymax=mean_values+SE_values),
                size=.3,    
                width=.5,
                position = position_dodge(width = 0.9, preserve = "single")) +
  labs(y = "Partial correlation coefficient") +
  coord_cartesian(ylim=c(-0.60,0.25)) +
  scale_fill_manual(values=plot_colors) +
  scale_y_continuous(breaks = c(-0.45,-0.30,-0.15,0,0.15)) +
  geom_hline(aes(yintercept=0), size=.3) +
  theme(aspect.ratio = 1, 
        axis.title.y=element_text(size=11, vjust=1),
        axis.text.y=element_text(size=7),
        axis.ticks.y=element_line(size=.3),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, size=7),
        axis.ticks.x=element_line(size=.3),
        legend.title=element_blank(),
        legend.text=element_text(size=7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black')
  ) 
sup_fig_5
ggsave(filename = "FigureS2.jpeg",
       device = "jpeg",
       width = 5.8*4, units = "cm",
       dpi = 600)


##----------------------------------------
## Linear regression
# DoY_off = autumn leaf senescence dates vs.
# cA_tot = photosynthesis activity accounting for water deficit
# with site (PEP_ID) as random effect
library(lme4)
library(MuMIn)
library(broom)

# Subset according ot species
all.species <- unique(drivers.df$Species)

output <- data.frame()
for (species in all.species){
  
  # Subset according to species
  sub.sp <- drivers.df %>% 
    select(Species,PEP_ID,DoY_off,cA_tot) %>% 
    filter(Species==species)
  
  # Fit the model
  f1 <- lmer(DoY_off~cA_tot+(1|PEP_ID), data=sub.sp)
  
  # Predict the data, removing the grouping value
  fix.pred1 <- predict(f1,re.form=NA)
  
  # To remove the random effects term and view the original spread of the data
  # with just the random noise added, add the residual error back 
  # onto the predicted values
  y.adj1 <- fix.pred1+resid(f1)
  
  # Store results
  sub.sp$cA_tot_random_effect <- y.adj1
  
  output <- rbind(output,sub.sp)
  print(species)
}

# Calculate r2s with site as random effect
fit.NPP <- lmer(DoY_off~cA_tot+(1|PEP_ID), data=output)
r2s <- as.data.frame(r.squaredGLMM(fit.NPP))
fit.CO2 <- lmer(DoY_off~ci+(1|PEP_ID), data=output)
r2_CO2 <- as.data.frame(r.squaredGLMM(fit.CO2))[1]
fit.tempaut <- lmer(DoY_off~temp_aut2+(1|PEP_ID), data=output)
r2_tempaut <- as.data.frame(r.squaredGLMM(fit.tempaut))[1]
fit.tempGS <- lmer(DoY_off~temp_GS+(1|PEP_ID), data=output)
r2_tempGS <- as.data.frame(r.squaredGLMM(fit.tempGS))[1]
fit.prec <- lmer(DoY_off~HRD+(1|PEP_ID), data=output)
r2_prec <- as.data.frame(r.squaredGLMM(fit.prec))[1]

# Plot 
# FIGURE 1A
# Define palettes
paletteBlueRed <- c("blue3","white","red3")

fig_1a <- ggplot(output, aes(x=cA_tot, y=cA_tot_random_effect)) +
  stat_bin2d(bins=300) +
  xlab(bquote('Photosynthesis rate (gC '~m^-2~year^-1*')')) +
  ylab("Autumn senescence (days)") +
  coord_cartesian(ylim=c(200,350))+
  stat_smooth(se=F, method="lm", colour="black", size=0.5) +
  scale_fill_gradientn(colours = paletteBlueRed,
                       limits = c(0,250),
                       breaks = c(5,125,250))+
  theme(aspect.ratio = 1,
        legend.position=c(0.15,0.2),
        legend.key.size=unit(0.35,"cm"),
        legend.title=element_blank(),
        axis.title.y=element_text(size=11),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=11),
        axis.text.x=element_text(size=9),
        strip.text=element_text(face="italic"),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black')
  )
fig_1a <- fig_1a + geom_text(data=r2s,
                             aes(x=4500,y=345,
                                 label=paste0("R2 = ",round(r2s[1,1],2))
                             ))
fig_1a

# Combine Fig. 1A and Fig. 1B
library(ggpubr)
fig_1ab <- ggarrange(fig_1a,fig_1b,
                     nrow=1,labels=c("a","b"),align="h")
fig_1ab
ggsave(filename = "Figure1AB.jpeg",
       device = "jpeg",
       width = 5.8*4, units = "cm",
       dpi = 600)


##----------------------------------------
## Path analysis 

# DoY_off = autumn leaf senescence dates 
# DoY_out = spring leaf-out dates (Fu et al. 2014)
# temp_GS = average temperature during the growing season (Xie et al. 2018; Liu et al. 2019)
# RD_summer = rainy days during summer (water supply during the driest season)
# temp_aut2 = average minimum temperature during autumn (2 months before leaf senescence)
# cA_tot = photosynthesis activity

# Use Structural Equation Modelling (SEM)
library(piecewiseSEM)
library(nlme)

# List structured equations
# Model including only climatic variables
climate_autumn_model <- piecewiseSEM::psem(
  
  lme(autumn_anomaly ~  spring_anomaly + tempaut2_anomaly + tempGS_anomaly + prec_anomaly, random = ~ 1 | PEP_ID, na.action = na.omit,
      data = drivers_sel_anomalies.df),
  
  lme(prec_anomaly ~ tempGS_anomaly, random = ~ 1 | timeseries, na.action = na.omit,
      data = drivers_sel_anomalies.df)
  
)

# Model including seasonal photosynthesis
# and climatic variables as indirect effects
full_autumn_model <- piecewiseSEM::psem(
  
  lme(NPP_anomaly ~ spring_anomaly + tempGS_anomaly + prec_anomaly + ci_anomaly, random = ~ 1 | PEP_ID, na.action = na.omit,
      data = drivers_sel_anomalies.df),
  
  lme(autumn_anomaly ~ NPP_anomaly + tempaut2_anomaly, random = ~ 1 | timeseries, na.action = na.omit,
      data = drivers_sel_anomalies.df)
  
)

# Summary of statistics
summary(climate_autumn_model, .progressBar = FALSE)
summary(full_autumn_model, .progressBar = FALSE)

# Compare models
anova(climate_autumn_model,full_autumn_model)

# Extract standardized parameter estimates
climate_stdest <- piecewiseSEM::coefs(climate_autumn_model)$Std.Estimate
full_stdest <- piecewiseSEM::coefs(full_autumn_model)$Std.Estimate

# Extract R2
climate_r2 <- piecewiseSEM::rsquared(climate_autumn_model)
full_r2 <- piecewiseSEM::rsquared(full_autumn_model)

# Plot
# FIGURE 1C
library(diagram)
jpeg("Figure1C.jpeg",width=5.8*6,height=5.8*3,units="cm",res=600)

# Plot flow-charts
par(mar=c(1,1,1,1), mfrow=c(1,2))

# Flow-chart with all significant environmental cues as predictors of leaf senescence
openplotmat()

# Add title-label
text(0.025,0.95,"c",cex=2.25,font=2)

# Define names of variables
names <- c("Spring\nAnomaly","Autumn\nTemperature","Summer\nTemperature","Summer\nPrecipitation",
           "Autumn\nAnomaly")

# Define coordinates
elpos <- coordinates(c(4,1))

# Define colors
Cols = c("blue3","blue3","red3","red3")

# Arrows
arrpos <- matrix(ncol=2, nrow=4)
for (i in 1:4) {
  arrpos[i,] <- straightarrow(from=elpos[i,],to=elpos[5,],lty=1,lcol=Cols[i],lwd=20*abs(climate_stdest)[i])
}

# Boxes
for(i in 1:5) {
  textrect(elpos[i,],0.12,0.05,lab=names[i],cex=1.25)
}

# Arrow-labels
for(i in 1:4) {
  text(arrpos[i,1]+0.05,arrpos[i,2],climate_stdest[i],col=Cols[i])
}

# Add r2
text(arrpos[i,1]+0.05,0.26,paste0("R2 = ",round(climate_r2[1,5],2)),cex=1.25)


# Flow-chart with seasonal photosynthesis and autumn temperature as predictors of leaf senescence 
openplotmat()

# Define coordinates
elpos <- coordinates(c(4,2,1))

# Define names of variables
names <- c("Spring\nAnomaly","Summer\nTemperature","Summer\nPrecipitation","CO2\nConcentration",
           "Photosynthesis","Autumn\nTemperature",
           "Autumn\nAnomaly")

# Define colors
Cols = c("red3","blue3","blue3","blue3",
         "red3","blue3")

# Arrows
arrpos <- matrix(ncol=2, nrow=6)
for (i in 1:4) {
  arrpos[i,] <- straightarrow(from=elpos[i,],to=elpos[5,],lty=1,lcol=Cols[i],lwd=20*abs(full_stdest)[i])
}
for (i in 5:6) {
  arrpos[i,] <- straightarrow(from=elpos[i,],to=elpos[7,],lty=1,lcol=Cols[i],lwd=20*abs(full_stdest)[i])
}

# Boxes
for(i in 1:7) {
  textrect(elpos[i,],0.12,0.05,lab=names[i],cex=1.25)
}

# Arrow-labels
for(i in 1:6) {
  text(arrpos[i,1]+0.05,arrpos[i,2],round(full_stdest[i],2),col=Cols[i])
}

# Add r2
text(arrpos[i,1]+0.12,0.18,paste0("R2 = ",round(full_r2[1,5],2)),cex=1.25)
dev.off()