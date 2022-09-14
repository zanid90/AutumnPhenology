#---
#title: Increased growing-season productivity drives earlier autumn leaf senescence in temperate trees (Zani et al. 2020 Science)
#author: Constantin Zohner, Deborah Zani
#date: "last updated May 11, 2022"

#R code to download, format, and clean PEP725 data
 

 
## 1. Load packages and set directories
  

require(devtools)
require(phenor)
require(data.table)
require(reshape2)
require(dplyr)
require(ggplot2)
require(data.table)
require(maptools)
require(maps)

# set the working directory
wd  = "/Users/consti/Desktop/PhD/Publication_material/17_Autumn_phenology_tier2/PEP_Analysis/Analysis/PEP_data"
wd2 = "/Users/consti/Desktop/PhD/Publication_material/000_Zani_et_al_reanalysis/2_Combined_analysis/Analysis_input"

# paths
tmp_path = paste(wd, "PEPzip", sep="/")
PEPcsv_path = paste(wd, "PEPcsv", sep="/")


  
## 2. Download and merge PEP data	
  

## Download PEP data from website

# load or download necessary data
# [create a proper pep725login.txt file first]
#check manually for errors. These are 1) umlauts and 2) missing phenology file info
#download_pep725(credentials = "PEPcredentials/PEPcredentials.txt",
#                species= check_pep725_species(list = TRUE)$number,#download everything
#                internal = F,
#                path = tmp_path)

#merge PEP files into one csv (check errors first)
tidy_pep_data = merge_pep725(path=tmp_path)


  
## 3. Format and clean PEP data
  

#keep only leaf-out and leaf senescence observations
tidy_pep_data = tidy_pep_data %>%
  filter(bbch %in% c("10","11","13","94","95")) %>% 
  mutate(day=as.numeric(day))

#check observed phenophases
barplot(table(tidy_pep_data$bbch), las=2)

#reshape table to short format (each bbch gets separate column)
tidy_pep_data = dcast(tidy_pep_data, pep_id+year+country+species+lat+lon+alt ~ bbch, value.var = "day")

####################
# Table operations #
####################

tidy_pep_data = tidy_pep_data %>%
  #rename bbch columns
  rename(leaf_out_10 = "10", leaf_out_11 = "11", leaf_out_13 = "13", 
         leaf_off_94 = "94", leaf_off_95 = "95") %>%
  
  #keep only individuals that have leaf-out and leaf-off information
  filter(complete.cases(leaf_off_94 | leaf_off_95)) %>%
  filter(complete.cases(leaf_out_10 | leaf_out_11 | leaf_out_13)) %>%
  
  #keep only years after 1947
  filter(year>1947) %>% 
  
  #group by species and site
  group_by(species, pep_id) %>% 
  
  #create leaf-out and leaf-off columns: which are the most common spring/autumn bbch codes within each time series?
  mutate(leaf_out_bbch = ifelse(sum(!is.na(leaf_out_11)) < sum(!is.na(leaf_out_13)), "use_bbch_13", "use_bbch_11"),
         leaf_off_bbch = ifelse(sum(!is.na(leaf_off_94)) < sum(!is.na(leaf_off_95)), "use_bbch_95", "use_bbch_94"),
         leaf_out = ifelse(species == "Larix decidua", leaf_out_10, #bbch 10 for larix
                           ifelse(leaf_out_bbch == "use_bbch_13", leaf_out_13, leaf_out_11)),
         leaf_off = ifelse(leaf_out_bbch == "use_bbch_95", leaf_off_95, leaf_off_94)) %>% 
  
  #delete NAs
  filter(!is.na(leaf_off),
         !is.na(leaf_out)) %>% 
  
  ##############
  #Data cleaning
  ##############

  #delete dates deviating from median more than 3 times MAD
  filter(!(abs(leaf_off-median(leaf_off))>3*mad(leaf_off) | 
           abs(leaf_out-median(leaf_out, na.rm=T))>3*mad(leaf_out, na.rm=T))) %>% 
  #delete time-series with standard deviation of leaf-off dates > 20
  filter(!(sd(leaf_off, na.rm=T)>20 | sd(leaf_out, na.rm=T)>15)) %>% 
  #delete time series with fewer than 15 years
  filter(n() >= 15) %>% 
  #delete time series with leaf-out > solstice
  filter(!max(leaf_out)>171) %>%
  
  mutate(
    #add mean leaf-out date
    leaf_out_mean = round(mean(leaf_out)),
    #add mean leaf-off date
    leaf_off_mean = round(mean(leaf_off))) %>%
  ungroup() %>% 
  
  mutate(
    #rename species
    species = case_when(species == "Aesculus hippocastanum" ~ "Aesculus",
                        species == "Betula(Betula pendula_(B._verrucosa|_B._alba))" ~ "Betula",
                        species == "Fagus(Fagus sylvatica)" ~ "Fagus",
                        species == "Larix decidua" ~ "Larix",
                        species == "Quercus robur_(Q.peduncula)" ~ "Quercus",
                        species == "Sorbus aucuparia" ~ "Sorbus"),
    #add timeseries identifier
    timeseries = paste0(pep_id, "_", species)) %>%
  
  #keep Aesculus, Betula, Fagus, Larix, Quercus, and Sorbus
  filter(species %in% c('Aesculus', 'Betula', 'Fagus', 'Larix', 'Quercus', 'Sorbus')) %>%
  
  #delete columns
  select(-c(leaf_out_10, leaf_out_11, leaf_out_13, leaf_off_94, leaf_off_95, leaf_out_bbch, leaf_off_bbch)) #delete columns

#Safe table
write.csv(tidy_pep_data, paste(wd2, "DataMeta_2_PhenologyObs_PEP725_CleanData.csv", sep="/"))


  
## 4. Sample sizes
  

#total observations
nrow(tidy_pep_data)

#how many species have data?
length(unique(tidy_pep_data$species))

#how many sites in total?
length(unique(tidy_pep_data$pep_id))

#time span 
range(tidy_pep_data$year)
hist(tidy_pep_data$year, xlab="Year", main="Temporal distribution of data")

#latitudinal gradient 
range(tidy_pep_data$lat)
hist(tidy_pep_data$lat, xlab="Latitude", main="Latitudinal gradient")

#elevational gradient 
range(tidy_pep_data$alt)
hist(tidy_pep_data$alt, xlab="Elevation (m)", main="Elevational gradient")

#leaf-out data
mean(tidy_pep_data$leaf_out)
sd(tidy_pep_data$leaf_out)
range(tidy_pep_data$leaf_out)
hist(tidy_pep_data$leaf_out, xlab="Leaf-out date", main="Leaf-out gradient")

#leaf-off data
mean(tidy_pep_data$leaf_off)
sd(tidy_pep_data$leaf_off)
range(tidy_pep_data$leaf_off)
hist(tidy_pep_data$leaf_off, xlab="Senescence date", main="Senescence gradient")

#Create summary dataframe
sample.size = tidy_pep_data %>%
  group_by(species) %>%
  summarize(n.time.series = length(unique(pep_id)),
            n.observations = n())

#Sample sizes within species
data.frame(sample.size)

#Total amount of times series
sum(sample.size$n.time.series)

#Bar Plot 
sample.size = melt(sample.size, id.vars = c("species"), measure.vars = c("n.time.series","n.observations"))
ggplot(sample.size, aes(x=species, y=value)) + 
  geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~variable, nrow=2, scales="free_y")

#Map the observations
mp <- NULL
mapWorld <- borders("world", colour="gray60", fill="gray60") # create a layer of borders
mp <- ggplot() + mapWorld

#Now Layer the stations on top
mp <- mp + geom_point(data=tidy_pep_data[!duplicated(tidy_pep_data[ , c("lat", "lon")]), ], 
                      aes(x=lon, y=lat) ,color="blue", size=.3) +
  coord_cartesian(ylim = c(43, 69), xlim = c(-10, 31))
mp


  
## 5. Reproducibility
  

## date time
Sys.time()

## session info
sessionInfo()


