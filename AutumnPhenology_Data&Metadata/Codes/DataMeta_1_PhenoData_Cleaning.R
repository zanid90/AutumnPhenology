################
## Data cleaning

# From the text:
# Following ref(43), we removed 
# (i) individual time-series with less than 15 years of leaf-out and leaf senescence observations for the same year, 
# (ii) dates deviating from an individual's median more than 3 times the median absolute deviation, 
# (iii) time-series for which the standard deviation of phenological observations across years was higher than 15 for leaf flushing and 20 for leaf senescence. 
# The thresholds differed because the mean absolute deviation for leaf-out time-series (8.8) was lower than that for leaf senescence (10.0). 

# Define directory path
setwd(".../AutumnPhenology/AutumnPhenology_Data&Metadata/Data/")

# Load libraries
library(data.table)
library(tidyverse)

##----------------------------------------
# Import data
# Acknowledgement: Data were provided by the members of the PEP725 project.
# http://www.pep725.eu/index.php
# PEP725 Pan European Phenology Data. Data set accessed 2019-02-08.
pheno_raw.df <- fread("PEP725_RawData.csv")
# Contains:
# Species ... the scientific plant name
# PEP_ID ... station identification code (see "PhenologyObs_PEP725_Stations.csv")
# BBCH ... leaf phenology BBCH Code (as defined by PEP725_BBCH.csv in PEP725 Datataset; see below) 
# YEAR ... year of the phenological observation
# DAY ... date of the year for a phenological observation

# Eliminate species with less than 3000 observations
table(pheno_raw.df$Species)
pheno_raw.df <- pheno_raw.df[!(pheno_raw.df$Species=="Acer pseudoplatanus"), ]

# Keep only records on leaf flushing and leaf senescence dates
# From PEP725_BBCH.csv (with definition of the BBCH Code) from the PEP725 Dataset
# Spring phenology codes
# 10	First leaves seperated (mouse ear)
# 11	Leaf unfolding (first visible leaf stalk)
# Autumn phenology codes
# 94	Autumnal coloring of leaves (50%)
# 95	50 % autumnal leaf fall
pheno_out.df <- pheno_raw.df[pheno_raw.df$BBCH==10 | pheno_raw.df$BBCH==11,] 
pheno_off.df <- pheno_raw.df[pheno_raw.df$BBCH==94 | pheno_raw.df$BBCH==95,]

# Create leaf phenology identifier
pheno_out.df$phenology <- "leaf.out"
pheno_off.df$phenology <- "leaf.off"

# Create timeseries identifier (siteXspecies)
pheno_out.df$timeseries <- paste0(pheno_out.df$PEP_ID,"_",pheno_out.df$Species)
pheno_off.df$timeseries <- paste0(pheno_off.df$PEP_ID,"_",pheno_off.df$Species)

# Keep only timeseries with observations for both spring and autumn leaf phenology 
ts_common <- intersect(pheno_out.df$timeseries, pheno_off.df$timeseries)
pheno_out.df <- pheno_out.df[pheno_out.df$timeseries %in% ts_common,] 
pheno_out.df <- pheno_off.df[pheno_off.df$timeseries %in% ts_common,]


##----------------------------------------
# Filter data

# Function to calculate MAD_ratio for each phenological observation
madratio_fx <- function(x) {
  abs(median(x)-x) / mad(x)
}

# Create data.frame to store output
pheno_clean.df <- data.frame()

for(ts in ts_common) {
  
  # Subset per timeseries
  pheno_out.sub <- pheno_out.df[pheno_out.df$timeseries==ts, ]
  pheno_off.sub <- pheno_off.df[pheno_off.df$timeseries==ts, ]
  
  # Eliminate duplicate phenological observations per year
  pheno_out.sub <- pheno_out.sub %>% 
    distinct(YEAR, .keep_all = TRUE)
  pheno_off.sub <- pheno_off.sub %>% 
    distinct(YEAR, .keep_all = TRUE)
  
  # Remove observations with MAD_ratio >=3
  pheno_out.sub <- pheno_out.sub %>% 
    mutate(MAD_ratio = madratio_fx(.$DAY)) %>% 
    filter(MAD_ratio<3)
  pheno_off.sub <- pheno_off.sub %>% 
    mutate(MAD_ratio = madratio_fx(.$DAY)) %>% 
    filter(MAD_ratio<3)
  
  # Keep only common observations (i.e. for the same years)
  yr_common <- intersect(pheno_out.sub$YEAR, pheno_off.sub$YEAR)
  pheno_out.sub <- pheno_out.sub[pheno_out.sub$YEAR %in% yr_common,]
  pheno_off.sub <- pheno_off.sub[pheno_off.sub$YEAR %in% yr_common,]
  
  # Error check
  if(nrow(pheno_out.sub)!=nrow(pheno_off.sub)) {
    print("ERROR! Unequal timeseries length")
    break
  }
  
  # Keep timeseries with more than 15 years of common observations
  if(nrow(pheno_out.sub)>=15) {
    
    # Calculate mean and standard deviation
    pheno_out.sub$mean_DoY <- mean(pheno_out.sub$DAY)
    pheno_out.sub$SD_DoY <- sd(pheno_out.sub$DAY)
    pheno_off.sub$mean_DoY <- mean(pheno_off.sub$DAY)
    pheno_off.sub$SD_DoY <- sd(pheno_off.sub$DAY)
    
    # Remove timeseries with standard deviation >15 for leaf flushing 
    # and >20 for leaf senescence
    if(sd(pheno_out.sub$DAY)<15 && sd(pheno_off.sub$DAY)<20) {
      
      # Timeseries length
      pheno_clean.sub <- rbind(pheno_out.sub,pheno_off.sub)
      pheno_clean.sub$timeseries_length <- nrow(pheno_clean.sub)
      
      # Bind to output dataset
      pheno_clean.df <- rbind(pheno_clean.df,pheno_clean.sub)
      print(paste0("Timeseries ",ts," included"))
    } else {
      print(paste0("Timeseries ",ts," excluded"))
    }
  } else {
    print(paste0("Timeseries ",ts," excluded"))
  }
}

##----------------------------------------
# Data wrangling 
pheno_clean.df <- pheno_clean.df %>% 
  rename(DoY = DAY) %>% 
  select(timeseries, everything()) %>% 
  arrange(timeseries,YEAR)


##----------------------------------------
# Summary of phenological observations

# Observations count (per species) 
(nrow(pheno_clean.df)/2)
pheno_clean.df %>% 
  group_by(Species) %>% 
  tally(name="Observations") %>% 
  mutate(Observations=Observations/2)

# Timeseries count
(length(unique(pheno_clean.df$timeseries)))

# Sites count (per species)
(length(unique(pheno_clean.df$PEP_ID)))
pheno_clean.df %>% 
  count(timeseries,Species) %>% 
  count(Species,name="Sites")


##----------------------------------------
# Export dataset
write.table(pheno_clean.df,"DataMeta_2_PhenologyObs_PEP725_CleanData.csv",sep=";",row.names=F)