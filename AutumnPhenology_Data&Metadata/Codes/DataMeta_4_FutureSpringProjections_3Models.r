library(parallel)
library(phenor)
library(data.table)
library(raster)
library(stringr)
# setwd("/Volumes/CrowtherLabRAID/Lidong_Mo/PhenologyFutureHPC/")
source("phenoModels/format_cmip5LM.R")

# pep725_dataOld = readRDS(paste("phenoModels/pep725_for_PhenologyFuture.rds", sep="/"), refhook = NULL)
pep725_data = readRDS(paste("phenoModels/pep725.rds", sep="/"), refhook = NULL)

# # load the ID and coordinates look up table
# lookupTable = fread("phenoModels/Reference_Coords_ID_DF_for_AutaumPhenology.csv")[,-1]

# lookupTable$DiffCoord = paste(lookupTable$lat,lookupTable$lon,sep="_")
# # get the different coord combinates into a vector
# diffCoordVector = unique(lookupTable$DiffCoord)

# for (coo in diffCoordVector)
# {
# 	# get the order of the coordinate comibination in the coordinates vector
# 	orderInfo = which(diffCoordVector==coo)
# 	lookupTable[lookupTable==coo]  = paste("Location_",str_pad(orderInfo,5, pad = "0"),sep="")
# }

# write.csv(lookupTable,"phenoModels/Reference_Coords_ID_DF_for_AutaumPhenology_Updated.csv")


lookupTable = fread("phenoModels/Reference_Coords_ID_DF_for_AutaumPhenology_Updated.csv")[,-1]
# get the difference location IDs
IDs = unique(lookupTable$DiffCoord)

# 
# set.seed(1000)
# IDs                      <- sample(IDs,size=420)
# write the function for one ID
FutureProjection           <- function(id)
{
	# if(!dir.exists(paste("phenoModels/FormatedCMIP/",id,sep="")))
	# {
	# 	# check the existence of the folder, if not create a new one
	# 	dir.create(paste("phenoModels/FormatedCMIP/",id,sep=""))
	# }

	# subset the lat lon information by the id
	LocalCoord = lookupTable[lookupTable$DiffCoord==id,]
	# get the location vector of that ID 
	CoordIn = c(unique(LocalCoord$lat),unique(LocalCoord$lon))
	# format future data
	# for (rcp in c(85))
	# {
	# 	for (yr in 2016:2100)
	# 	{
	# 		FormatedCmipData = format_cmip5LM(path = "phenoModels/FutureClimateData",
 #                                              year = yr,
 #                                              offset = 264,
 #                                              model = "MIROC5",
 #                                              scenario = paste("rcp",rcp,sep=""),
 #                                              extent = c(-40.5,75.5,25,75.5),
 #                                              internal = TRUE,
	# 					                      coord= CoordIn)
	# 		# save the formated future data into the local folder which has the name 
	# 		saveRDS(FormatedCmipData, paste("phenoModels/FormatedCMIP/",id,"/",id,"_Formated_cmip_RCP",rcp,"_",yr,".rds",sep=""))
	# 		print(paste("--- Data formtting for ",id," of RCP",rcp," has been done for Year ",yr,"---", sep=""))
	# 	}
	# }

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
	        for (rcp in c(85))
	        {   
			    # init a vector
			    EmptyVec               <- vector()
			    for (yr in c(2016:2100))
	            {
				    # read the rds data for each in each senario 
				    FormatedCmipData   <- readRDS(paste("phenoModels/FormatedCMIP/",id,"/",id,"_Formated_cmip_RCP",rcp,"_",yr,".rds",sep=""))
                    # execute the projection
				    ProjectedRaster    <- estimate_phenology(data = FormatedCmipData,
                                                         par = optim.par$par,
							                             model=mdl)
                    # as the result from the last step is a one cell raster, we need to transfer the value to a numeric value
				    ProjectedValue     <- ProjectedRaster@data@values
				    # allocate this value to a vector
				    EmptyVec           <- append(EmptyVec,ProjectedValue)
	            }
			    OneRowDF               <- data.frame(site,ReferenceID,t(Coordinates),mdl,RCP=paste("RCP",rcp,sep=""),t(EmptyVec))
			    RowsDF                 <- rbind(RowsDF,OneRowDF) 
			    print(paste(mdl,"_",yr,"_",rcp,"_",ReferenceID,sep=""))  
		    }
	    }
	    # allocat names to the data frame
	    names(RowsDF)                  <- c("Site","ReferenceID","Lat","Lon","Model","Senarios",paste("Year_",2016:2100,sep=""))       
	    # return the caclulation data frame
	    write.csv(RowsDF,paste("phenoModels/SplitedIDResult/",id,"_",site,"_Formated_cmip_RCP",rcp,"_Projection.csv",sep=""))	
    }
	# unlink(paste("phenoModels/FormatedCMIP/",id,"/",sep=""), recursive = T)	
}

system.time(mclapply(IDs,FutureProjection,mc.preschedule=F,mc.cores=12))


# Befeore we start the merge of the data we need to cd to the directory we stroed our data in HPC
#  cd  /nfs/nas22.ethz.ch/fs2201/usys_ibz_cr_lab/Lidong/PhenologyFutureHPC
library(data.table)
library(dplyr)


# merge all the tables in the country tables folder
tableList = list.files(path ="phenoModels/SplitedIDResult/",pattern=".csv")
# rbind all those tables
rbindedTable = rbindlist(lapply(paste("phenoModels/SplitedIDResult/",tableList,sep=""),fread))
# write the data frame into the local folder
write.csv(rbindedTable,"Future_Projection_Table_For_Autaum_Phenology.csv")