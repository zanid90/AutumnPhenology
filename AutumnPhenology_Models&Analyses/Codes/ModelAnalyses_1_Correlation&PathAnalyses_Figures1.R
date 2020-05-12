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
variables <- fread("DataMeta_3_Drivers.csv")


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

## Photosynthesis activity
# cGSI = cumulative/annuals Growing Season Index (function of temperature, day-length and vapour pressure deficit)
# cA_tot-w = annual net photosynthetic activity
# cA_tot = accounting for water supply during the growing season

# Compute a correlation matrix
corr <- variables %>% 
  select(-c(timeseries,PEP_ID,LON,LAT,Species,YEAR,
            autumn_anomaly,spring_anomaly,temp_aut3)) %>% 
  cor() %>% # default method = "pearson"
  round(.,1)

# Visualize the correlation matrix
# EXTENDED DATA FIGURE 1
corr_plot <- ggcorrplot(corr,
                        type = "lower", # Get the lower triangle
                        colors = c("red3", "white", "blue3"), 
                        legend.title = "",
                        outline.col = "white",
                        lab = TRUE) # Add correlation coefficients
corr_plot


##----------------------------------------
## (Partial) Correlations
library(ppcor)

# DoY_off = autumn leaf senescence dates vs.

# DoY_out = spring leaf-out dates (Fu et al. 2014)
# temp_GS = average temperature during the growing season (Xie et al. 2018; Liu et al. 2019)
# RD_summer = rainy days during summer (water supply during the driest season)
# cGSI = cumulative/annuals Growing Season Index (function of temperature, day-length and vapour pressure deficit)
# cA_tot-w = annual net photosynthetic activity
# cA_tot = photosynthesis activity accounting for water deficit

# accounting for temp_aut2 = average minimum temperature during autumn (2 months before leaf senescence)

# Time-series level
pcorr.tm <- variables %>%
  group_by(timeseries) %>% 
  summarise(pcorr_cAtot = pcor.test(DoY_off, cA_tot, temp_aut2)[,1],
            `pcorr_cA_tot-w` = pcor.test(DoY_off, `cA_tot-w`, temp_aut2)[,1],
            pcorr_cGSI = pcor.test(DoY_off, cGSI, temp_aut2)[,1],
            pcorr_tempGS = pcor.test(DoY_off, temp_GS, temp_aut2)[,1],
            pcorr_RDsummer = pcor.test(DoY_off, RD_summer, temp_aut2)[,1],
            pcorr_DoYout = pcor.test(DoY_off, DoY_out, temp_aut2)[,1]
            )

# Average values across species
meanpcorr_across_species <- pcorr.tm %>% 
  summarise(`cA_tot` = mean(pcorr_cAtot),
            `cA_tot-w` = mean(`pcorr_cA_tot-w`),
            `cGSI` = mean(pcorr_cGSI),
            `Summer temperature` = mean(pcorr_tempGS),
            `Summer precipitation` = mean(pcorr_RDsummer),
            `Spring leaf-out` = mean(pcorr_DoYout)) %>% 
  add_column(Species = "Across species") %>% 
  pivot_longer(-Species,names_to="predictors",values_to = "mean_values")

# Species-specific time-series level
pcorr.sp <- variables %>%
  group_by(Species, timeseries) %>% 
  summarise(pcorr_cAtot = pcor.test(DoY_off, cA_tot, temp_aut2)[,1],
            `pcorr_cA_tot-w` = pcor.test(DoY_off, `cA_tot-w`, temp_aut2)[,1],
            pcorr_cGSI = pcor.test(DoY_off, cGSI, temp_aut2)[,1],
            pcorr_tempGS = pcor.test(DoY_off, temp_GS, temp_aut2)[,1],
            pcorr_RDsummer = pcor.test(DoY_off, RD_summer, temp_aut2)[,1],
            pcorr_DoYout = pcor.test(DoY_off, DoY_out, temp_aut2)[,1]
            )

# Average values per species
meanpcorr_per_species <- pcorr.sp %>%
  summarise(`cA_tot` = mean(pcorr_cAtot),
            `cA_tot-w` = mean(`pcorr_cA_tot-w`),
            `cGSI` = mean(pcorr_cGSI),
            `Summer temperature` = mean(pcorr_tempGS),
            `Summer precipitation` = mean(pcorr_RDsummer),
            `Spring leaf-out` = mean(pcorr_DoYout)) %>% 
  pivot_longer(cols=-Species,names_to="predictors",values_to = "mean_values") %>% 
  # Add average values across species
  rbind(meanpcorr_across_species,.) 

# Standard errors per species
SE_per_species <- pcorr.sp %>% 
  group_by(Species) %>%
  summarise(`cA_tot` = se(pcorr_cAtot),
            `cA_tot-w` = se(`pcorr_cA_tot-w`),
            `cGSI` = se(pcorr_cGSI),
            `Summer temperature` = se(pcorr_tempGS),
            `Summer precipitation` = se(pcorr_RDsummer),
            `Spring leaf-out` = se(pcorr_DoYout)) %>%
  pivot_longer(cols=2:7,names_to="predictors",values_to="SE_values")

# Standard errors across species
SE_across_species <- pcorr.tm %>% 
  summarise(`cA_tot` = se(pcorr_cAtot),
            `cA_tot-w` = se(`pcorr_cA_tot-w`),
            `cGSI` = se(pcorr_cGSI),
            `Summer temperature` = se(pcorr_tempGS),
            `Summer precipitation` = se(pcorr_RDsummer),
            `Spring leaf-out` = se(pcorr_DoYout)) %>% 
  add_column(Species = "Across species") %>% 
  pivot_longer(-Species,names_to="predictors",values_to = "SE_values")
SE <- rbind(SE_across_species,SE_per_species)
meanpcorr_per_species$SE_values <- SE$SE_values

# Set the order of predictors to plot
meanpcorr_per_species$predictors <- factor(meanpcorr_per_species$predictors,
                                           levels=c("cA_tot",
                                                    "cA_tot-w",
                                                    "cGSI",
                                                    "Summer temperature",
                                                    "Summer precipitation",
                                                    "Spring leaf-out"))

# Define palette
paletteSpectral <- rev(c("#3288BD", "#66C2A5", 
                         "#FFFFBF", "#FEE08B", "#FDAE61", 
                         "#F46D43", "#D53E4F"))

# Plot 
# EXTENDED DATA FIGURE 2A
ext_fig_2a <- ggplot(data = meanpcorr_per_species, aes(x=Species, y=mean_values, fill=predictors)) +
  geom_bar(position = position_dodge(preserve = "single"),
           stat="identity") + 
  geom_errorbar(aes(ymin=mean_values-SE_values, ymax=mean_values+SE_values),
                size=.3,    
                width=.5,
                position = position_dodge(width = 0.9, preserve = "single")) +
  labs(y = "Partial correlation coefficient") +
  coord_cartesian(ylim=c(-0.60,0.25)) +
  scale_fill_manual(values=rev(paletteSpectral[c(1,3:7)])) +
  scale_y_continuous(breaks = c(-0.45,-0.30,-0.15,0,0.15)) +
  geom_hline(aes(yintercept=0), size=.3) +
  theme(aspect.ratio = 1, 
        axis.title.y=element_text(size=8, vjust=1),
        axis.text.y=element_text(size=6),
        axis.ticks.y=element_line(size=.3),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, size=6),
        axis.ticks.x=element_line(size=.3),
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black')
  ) 
ext_fig_2a


## Pearson correlations
# DoY_off = autumn leaf senescence dates vs.

# DoY_out = spring leaf-out dates (Fu et al. 2014)
# temp_GS = average temperature during the growing season (Xie et al. 2018; Liu et al. 2019)
# temp_aut2 = average minimum temperature during autumn (2 months before leaf senescence)
# RD_summer = rainy days during summer (water supply during the driest season)
# cGSI = cumulative/annuals Growing Season Index (function of temperature, day-length and vapour pressure deficit)
# cA_tot-w = annual net photosynthetic activity
# cA_tot = photosynthesis activity accounting for water deficit

# Species-specific time-series level
corr.sp <- variables %>%
  group_by(Species, timeseries) %>% 
  summarise(corr_cAtot = cor(DoY_off, cA_tot),
            `corr_cA_tot-w` = cor(DoY_off, `cA_tot-w`),
            corr_cGSI = cor(DoY_off, cGSI),
            corr_tempGS = cor(DoY_off, temp_GS),
            corr_RDsummer = cor(DoY_off, RD_summer),
            corr_tempaut2 = cor(DoY_off, temp_aut2),
            corr_DoYout = cor(DoY_off, DoY_out)
  )

# Average values per species
meancorr_per_species <- corr.sp %>%
  summarise(`cA_tot` = mean(corr_cAtot),
            `cA_tot-w` = mean(`corr_cA_tot-w`),
            `cGSI` = mean(corr_cGSI),
            `Summer temperature` = mean(corr_tempGS),
            `Summer precipitation` = mean(corr_RDsummer),
            `Autumn temperature` = mean(corr_tempaut2),
            `Spring leaf-out` = mean(corr_DoYout)) %>% 
  pivot_longer(cols=-Species,names_to="predictors",values_to = "mean_values") 

# Standard errors across timeseries
SE_per_species <- corr.sp %>% 
  group_by(Species) %>%
  summarise(`cA_tot` = se(corr_cAtot),
            `cA_tot-w` = se(`corr_cA_tot-w`),
            `cGSI` = se(corr_cGSI),
            `Summer temperature` = se(corr_tempGS),
            `Summer precipitation` = se(corr_RDsummer),
            `Autumn temperature` = se(corr_tempaut2),
            `Spring leaf-out` = se(corr_DoYout)) %>%
  pivot_longer(cols=2:8,names_to="predictors",values_to="SE_values")
meancorr_per_species$SE_values <- SE_per_species$SE_values

# Set the order of predictors to plot
meancorr_per_species$predictors <- factor(meancorr_per_species$predictors,
                                                levels=c("cA_tot",
                                                         "cA_tot-w",
                                                         "cGSI",
                                                         "Summer temperature",
                                                         "Summer precipitation",
                                                         "Autumn temperature",
                                                         "Spring leaf-out"))

# Plot 
# EXTENDED DATA FIGURE 2B
ext_fig_2b <- ggplot(data=meancorr_per_species, aes(x=Species, y=mean_values, fill=predictors)) +
  geom_bar(position = position_dodge(preserve = "single"),
           stat="identity") + 
  geom_errorbar(aes(ymin=mean_values-SE_values, ymax=mean_values+SE_values),
                size=.3,    
                width=.5,
                position = position_dodge(width = 0.9, preserve = "single")) +
  labs(y = "Correlation coefficient") +
  coord_cartesian(ylim=c(-0.60,0.25)) +
  scale_fill_manual(values=rev(paletteSpectral)) +
  scale_y_continuous(breaks = c(-0.45,-0.30,-0.15,0,0.15),
                     position = 'right') +
  geom_hline(aes(yintercept=0), size=.3) +
  theme(aspect.ratio = 1, 
        axis.title.y=element_text(size=8, vjust=1),
        axis.text.y=element_text(size=6),
        axis.ticks.y=element_line(size=.3),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, size=6),
        axis.ticks.x=element_line(size=.3),
        legend.title=element_blank(),
        legend.text=element_text(size=6),
        legend.key.size=unit(0.3,"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black')
  )
ext_fig_2b


## Linear regression
# DoY_off = autumn leaf senescence dates vs.
# cA_tot = photosynthesis activity accounting for water deficit
# with site (PEP_ID) as random effect
library(lme4)
library(MuMIn)

# Subset according ot species
all.species <- unique(variables$Species)

output <- data.frame()
for (species in all.species){
  
  # Subset according to species
  sub.sp <- variables %>% 
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

# Calculate r2 with site as random effect
fit.all <- lmer(DoY_off~cA_tot+(1|PEP_ID), data=output)
r2s <- as.data.frame(r.squaredGLMM(fit.all))

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


## Pearson correlations
# DoY_off = autumn leaf senescence dates vs.

# DoY_out = spring leaf-out dates (Fu et al. 2014)
# temp_GS = average temperature during the growing season (Xie et al. 2018; Liu et al. 2019)
# temp_aut2 = average minimum temperature during autumn (2 months before leaf senescence)
# RD_summer = rainy days during summer (water supply during the driest season)
# cGSI = cumulative/annuals Growing Season Index (function of temperature, day-length and vapour pressure deficit)
# cA_tot = photosynthesis activity accounting for water deficit

# Time-series level
corr.tm <- variables %>%
  group_by(timeseries) %>% 
  summarise(corr_cAtot = cor(DoY_off, cA_tot),
            corr_cGSI = cor(DoY_off, cGSI),
            corr_tempGS = cor(DoY_off, temp_GS),
            corr_RDsummer = cor(DoY_off, RD_summer),
            corr_tempaut2 = cor(DoY_off, temp_aut2),
            corr_DoYout = cor(DoY_off, DoY_out)
            )

# Average values across timeseries
meancorr_across_timeseries <- corr.tm %>% 
  summarise(`cA_tot` = mean(corr_cAtot),
            `cGSI` = mean(corr_cGSI),
            `Summer temperature` = mean(corr_tempGS),
            `Summer precipitation` = mean(corr_RDsummer),
            `Autumn temperature` = mean(corr_tempaut2),
            `Spring leaf-out` = mean(corr_DoYout)) %>%
  add_column(timeseries = "Across timeseries") %>% 
  pivot_longer(cols=-timeseries,names_to="predictors",values_to="mean_values")

# Standard errors across timeseries
SE_across_timeseries <- corr.tm %>% 
  summarise(`cA_tot` = se(corr_cAtot),
            `cGSI` = se(corr_cGSI),
            `Summer temperature` = se(corr_tempGS),
            `Summer precipitation` = se(corr_RDsummer),
            `Autumn temperature` = se(corr_tempaut2),
            `Spring leaf-out` = se(corr_DoYout)) %>%
  pivot_longer(cols=1:6,names_to="predictors",values_to="SE_values")
meancorr_across_timeseries <- meancorr_across_timeseries %>% 
  left_join(SE_across_timeseries,by="predictors") 

# Set the order of predictors to plot
meancorr_across_timeseries$predictors <- factor(meancorr_across_timeseries$predictors,
                                            levels=c("cA_tot",
                                                     "cGSI",
                                                     "Summer temperature",
                                                     "Summer precipitation",
                                                     "Autumn temperature",
                                                     "Spring leaf-out"))

# Plot 
# FIGURE 1B
# Define palettes
paletteBlueRed <- c("blue3","white","red3")
plot_colors <-  colorRampPalette(paletteBlueRed)(6)

fig_1b <- ggplot(data=meancorr_across_timeseries, aes(x=predictors, y=mean_values)) +
  geom_bar(position=position_dodge(), stat="identity", width=.9,
           fill=plot_colors) +
  geom_errorbar(aes(ymin=mean_values-SE_values, ymax=mean_values+SE_values),
                size=.3,    
                width=.5,
                position=position_dodge(.9)) +
  labs(y = "Correlation coefficient") +
  coord_cartesian(ylim=c(-0.55,0.25)) +
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


## Path analysis 

# DoY_off = autumn leaf senescence dates 
# DoY_out = spring leaf-out dates (Fu et al. 2014)
# temp_GS = average temperature during the growing season (Xie et al. 2018; Liu et al. 2019)
# RD_summer = rainy days during summer (water supply during the driest season)
# temp_aut2 = average minimum temperature during autumn (2 months before leaf senescence)
# cA_tot = photosynthesis activity

# Use Structural Equation Modelling (SEM)
library(lavaan)
library(diagram)
library(lme4)
library(MuMIn)

# List structured equations for lavaan
climate_autumn_model = '
DoY_off ~ DoY_out + temp_summer + RD_summer + temp_aut2
'

full_autumn_model = '
DoY_off ~ cA_tot + temp_aut2

cA_tot ~  DoY_out + temp_summer + RD_summer 
'

# Fit vcov SEM
climate_autumn_model.fit = sem(climate_autumn_model, variables, estimator = "MLM")
full_autumn_model.fit = sem(full_autumn_model, variables, estimator = "MLM")

# Extract standardized parameter estimates
climate_est.std <- standardizedSolution(climate_autumn_model.fit, type="std.all")$est.std[1:4]
climate_est.std <- round(climate_est.std,2)
full_est.std <- standardizedSolution(full_autumn_model.fit, type="std.all")$est.std[c(3:5,1:2)]
full_est.std <- round(full_est.std,2)

# Calculate r2 with site as random effect
fit.1 <- lmer(DoY_off~DoY_out+temp_summer+RD_summer+temp_aut2+(1|PEP_ID),data=variables)
r2.1 <- r.squaredGLMM(fit.1) 
r2.1 <- round(r2.1[1,1],2) # 0.13

fit.2 <- lmer(DoY_off~cA_tot+temp_aut2+(1|PEP_ID),data=variables)
r2.2 <- r.squaredGLMM(fit.2) 
r2.2 <- round(r2.2[1,1],2) # 0.35

# Plot
# FIGURE 1C

# Plot flow-charts
par(mar=c(1,1,1,1), mfrow=c(1,2))

# Flow-chart with all significant environmental cues as predictors of leaf senescence
openplotmat()

# Add title-label
text(0.025,0.95,"c",cex=2.25,font=2)

# Define names of variables
names <- c("Spring\nLeaf-out","Summer\nTemperature","Summer\nPrecipitation","Autumn\nTemperature",
           "Autumn\nSenescence")

# Define coordinates
elpos <- coordinates(c(4,1))

# Define correlation coefficients (from SEM analysis)
climate_est.std <- c(0.06,-0.24,-0.13,0.22)

# Define colors
Cols = c("blue3","red3","red3","blue3")

# Arrows
arrpos <- matrix(ncol=2, nrow=4)
for (i in 1:4) {
  arrpos[i,] <- straightarrow(from=elpos[i,],to=elpos[5,],lty=1,lcol=Cols[i],lwd=20*abs(climate_est.std)[i])
}

# Boxes
for(i in 1:5) {
  textrect(elpos[i,],0.12,0.05,lab=names[i],cex=1.25)
}

# Arrow-labels
for(i in 1:4) {
  text(arrpos[i,1]+0.05,arrpos[i,2],climate_est.std[i],col=Cols[i])
}

# Add r2
text(arrpos[i,1]+0.05,0.26,bquote(~R^2*' = '~0.13),cex=1.25)


# Flow-chart with seasonal photosynthesis (cAtot,water) and autumn temperature as predictors of leaf senescence 
openplotmat()

# Define coordinates
elpos <- coordinates(c(3,2,1))

# Define names of variables
names <- c("Spring\nLeaf-out","Summer\nTemperature","Summer\nPrecipitation",
           "Photosynthesis","Autumn\nTemperature",
           "Autumn\nSenescence")

# Define correlation coefficients (from SEM analysis)
full_est.std <- c(-0.17,0.44,0.43,-0.44,0.17)

# Define colors
Cols = c("red3","blue3","blue3",
         "red3","blue3")

# Arrows
arrpos <- matrix(ncol=2, nrow=5)
for (i in 1:3) {
  arrpos[i,] <- straightarrow(from=elpos[i,],to=elpos[4,],lty=1,lcol=Cols[i],lwd=20*abs(full_est.std)[i])
}
for (i in 4:5) {
  arrpos[i,] <- straightarrow(from=elpos[i,],to=elpos[6,],lty=1,lcol=Cols[i],lwd=20*abs(full_est.std)[i])
}

# Boxes
for(i in 1:6) {
  textrect(elpos[i,],0.12,0.05,lab=names[i],cex=1.25)
}

# Arrow-labels
for(i in 1:5) {
  text(arrpos[i,1]+0.05,arrpos[i,2],full_est.std[i],col=Cols[i])
}

# Add r2
text(arrpos[i,1]+0.12,0.18,bquote(~R^2*' = '~0.35),cex=1.25)