######################################
## Estimate future leaf flushing dates

# From the text:
# Future projections of leaf-out dates were obtained using 
# three of the best performing spring phenology models from the literature 
# AT: Cannell, M. G. R. & Smith, R. I. Thermal time, chill days and prediction of bud burst in Picea sitchensis. J. Appl. Ecol. 20, 951-963 (1983).
# M1: Hänninen, H. Modelling bud dormancy release in trees from cool and temperate regions. Acta For. Fenn. 213, 1-47 (1990).
# PM1: Blümel, K. & Chmielewski, F. M. Shortcomings of classical phenological forcing models and a way to overcome them. Agric. For. Meteorol. 164, 10-19 (2020). 

# Load libraries
library(data.table)
library(dplyr)
library(parallel)
library(phenor)
library(data.table)
library(raster)
library(stringr)
source("phenoModels/format_cmip5LM.R")

# Import data
pep725_data = readRDS(paste("phenoModels/pep725.rds", sep="/"), refhook = NULL)

# Load the ID and coordinates
lookupTable = fread("phenoModels/Reference_Coords_ID_DF_for_AutaumPhenology.csv")[,-1]

# Get the difference location IDs
IDs = unique(lookupTable$DiffCoord)

# Function to predict future spring phenology
FutureProjection <- function(id) {

	# Subset according to location (lon, lat in id)
	LocalCoord = lookupTable[lookupTable$DiffCoord==id,]
	
	# Get the location vector of that ID 
	CoordIn = c(unique(LocalCoord$lat),unique(LocalCoord$lon))

	subsetList = lapply(pep725_data,function(x) if(all(x$location==CoordIn)){return (x)})
	
	# kick out the null elements in the filtered list
	subsetList = Filter(Negate(is.null), subsetList)
	
	# get the name list in that id
  sites = names(subsetList)
  sites = LocalCoord$pep_id
  
  for (site in sites) 
  {
	    # subset the data frame for each site of each species
	    perDataFrame = pep725_data[[site]]
	    # get the coordinates of each site
	    Coordinates = perDataFrame$location
      # compare with the lookup table 
      SubsetLookup = lookupTable[lookupTable$lat==Coordinates[1]&lookupTable$lon==Coordinates[2],]
	    # get the reference ID
      ReferenceID = id
      # model parameter optimization
	    # as each mode has different parameter range we have to use the parameter reference table
	    ParameterRef = fread("phenoModels/parameter_ranges.txt",sep=",")
	    # models names vector
	    Models = c("M1", "AT", "PM1")
	    # make an empty data frame to save the data 
	    RowsDF = data.frame()
	    
	    for (mdl in Models)
	    {
		    # get the model lower and uper range
		    SubsetParamters = as.data.frame(ParameterRef[ParameterRef$model==mdl,][,-c("model","boundary")])
		    RetainedPar = as.matrix(SubsetParamters[,-which(apply(SubsetParamters,2,function(x)all(is.na(x))))])
		    print(RetainedPar)
		    set.seed(1000)
		    # cleaning the NA values in the data frame 
		    optim.par = optimize_parameters(par = NULL,
                        data = pep725_data[[site]],
                        cost = rmse,
                        model = mdl,
                        method = "GenSA",
						lower = as.numeric(RetainedPar[1,]),
                        upper = as.numeric(RetainedPar[2,]))
		    
		    # loop for rcp senarios and every future year
		    for (rcp in c(85)) {   
		      # init a vector
			    EmptyVec <- vector()
			    for (yr in c(2016:2100)) {
				    # read the rds data for each in each senario 
				    FormatedCmipData   <- readRDS(paste("phenoModels/FormatedCMIP/",id,"/",id,"_Formated_cmip_RCP",rcp,"_",yr,".rds",sep=""))
            # execute the projection
				    ProjectedRaster    <- estimate_phenology(data = FormatedCmipData,
				                                             par = optim.par$par,
				                                             model=mdl)
				    # as the result from the last step is a one cell raster, we need to transfer the value to a numeric value
				    ProjectedValue     <- ProjectedRaster@data@values
				    # allocate this value to a vector
				    EmptyVec <- append(EmptyVec,ProjectedValue)
				    }
			    
			    OneRowDF <- data.frame(site,ReferenceID,t(Coordinates),mdl,RCP=paste("RCP",rcp,sep=""),t(EmptyVec))
			    RowsDF <- rbind(RowsDF,OneRowDF) 
			    print(paste(mdl,"_",yr,"_",rcp,"_",ReferenceID,sep=""))  
		    }
	    }
	    
	    # allocat names to the data frame
	    names(RowsDF) <- c("Site","ReferenceID","Lat","Lon","Model","Senarios",paste("Year_",2016:2100,sep=""))  
	    
	    # return the caclulation data frame
	    write.csv(RowsDF,paste("phenoModels/SplitedIDResult/",id,"_",site,"_Formated_cmip_RCP",rcp,"_Projection.csv",sep=""))	
    }
	# unlink(paste("phenoModels/FormatedCMIP/",id,"/",sep=""), recursive = T)	
}

system.time(mclapply(IDs,FutureProjection,mc.preschedule=F,mc.cores=12))

# Merge all the tables in the country tables folder
tableList = list.files(path ="phenoModels/SplitedIDResult/",pattern=".csv")

# Bind all those tables
rbindedTable = rbindlist(lapply(paste("phenoModels/SplitedIDResult/",tableList,sep=""),fread))

# Export dataset
write.csv(rbindedTable,"DataMeta_5_FutureSpringProjections_3Models.csv")