# sets working directory
setwd("/home/jim/Dropbox/REM/tasks/scal_lengths")

# sources mseRtools.r, needed for lisread function()
source("mseRtools.r")

# sources df_to_dat.r, needed for df_to_dat function
source("df_to_dat.r")

# The following names should be used consistently throughout the various data files,
# ADMB, and R code: 
# Year RV_4VWX HS LL_NAFO3_Obs LL_NAFO4_Obs OT_NAFO3_Obs OT_NAFO4_Obs LL_PS OT_PS
# each of the three sample types has it's own lengthComp file with suffix _m, _f,
# or _c
allNames <- c("RV_4VWX", "HS", "LL_NAFO3_Obs", "LL_NAFO4_Obs", "OT_NAFO3_Obs", "OT_NAFO4_Obs")

# read in the catch data
scaHalCatch_Index_lb5_Sept152014 <- lisread("scaHalCatch_Index_lb5_ Sept152014 .dat.txt")
# pull out the landed catch matrix
landCatchMatrix <- scaHalCatch_Index_lb5_Sept152014[[2]]
landCatchMatrix <- as.data.frame(landCatchMatrix)

# set the names for the landed catch matrix
names_cat <- "TimeStep Year X3_long.line X3_other X3_otter.trawl X4_long.line X4_other X4_otter.trawl rv4VWX rv3NPO.1 rv3NPO.2 rv3NPO.3 hal.surv ci"
names_cat <- strsplit(names_cat," ")[[1]]

# create a copy of names_cat, will relabel the names with the names from the sample data for consistency
names_sam <- names_cat
for(i in seq(from=1,to=length(names_cat),by=1)){
    if(names_cat[i] == "X3_long.line"){names_sam[i] <- "LL_NAFO3_Obs"}
    else if(names_cat[i] == "X3_otter.trawl"){names_sam[i] <- "OT_NAFO3_Obs"}
    else if(names_cat[i] == "X4_long.line"){names_sam[i] <- "LL_NAFO4_Obs"}
    else if(names_cat[i] == "X4_otter.trawl"){names_sam[i] <- "OT_NAFO4_Obs"}
    else if(names_cat[i] == "rv4VWX"){names_sam[i] <- "RV_4VWX"}
    else if(names_cat[i] == "hal.surv"){names_sam[i] <- "HS"}
}

# set new names on landCatchMatrix to be in line with the names in the sample data
names(landCatchMatrix) <- names_sam

# set missing data to NA
landCatchMatrix[landCatchMatrix == -1] <- NA

# lump in OT with other
landCatchMatrix$OT_NAFO3_Obs <- landCatchMatrix$OT_NAFO3_Obs + landCatchMatrix$X3_other
landCatchMatrix$OT_NAFO4_Obs <- landCatchMatrix$OT_NAFO4_Obs + landCatchMatrix$X4_other

landCatchMatrix <- landCatchMatrix[,c("TimeStep","Year","LL_NAFO3_Obs", "LL_NAFO4_Obs", "OT_NAFO3_Obs", "OT_NAFO4_Obs","RV_4VWX", "HS")]

cf <- 1
for(cf in c(1,2,4)){
	if(cf == 1){
		landCatchMatrix_cf <- landCatchMatrix
		landCatchMatrix_cf$Commercial_fishery <- landCatchMatrix_cf$LL_NAFO3_Obs + landCatchMatrix_cf$LL_NAFO4_Obs + landCatchMatrix_cf$OT_NAFO3_Obs + landCatchMatrix_cf$OT_NAFO4_Obs
		landCatchMatrix_cf <- landCatchMatrix_cf[,c("TimeStep","Year","Commercial_fishery","RV_4VWX","HS")]
	} else if (cf == 2){
		landCatchMatrix_cf <- landCatchMatrix
		landCatchMatrix_cf$LL <- landCatchMatrix_cf$LL_NAFO3_Obs + landCatchMatrix_cf$LL_NAFO4_Obs
		landCatchMatrix_cf$OT <- landCatchMatrix_cf$OT_NAFO3_Obs + landCatchMatrix_cf$OT_NAFO4_Obs
		landCatchMatrix_cf <- landCatchMatrix_cf[,c("TimeStep","Year","LL","OT","RV_4VWX","HS")]
	} else if (cf == 4){
		landCatchMatrix_cf <- landCatchMatrix
		landCatchMatrix_cf <- landCatchMatrix_cf[,c("TimeStep","Year","LL_NAFO3_Obs", "LL_NAFO4_Obs", "OT_NAFO3_Obs", "OT_NAFO4_Obs","RV_4VWX", "HS")]
	}

	head(landCatchMatrix_cf)
	landCatchMatrix_cf[is.na(landCatchMatrix_cf) == TRUE] <- 1

	char <- paste("## scal_catch_cf1.dat
			##",df_to_dat(names(landCatchMatrix_cf)),
			"# landCatchMatrix(tonnes) \n",
			df_to_dat(landCatchMatrix_cf),sep="")

	# writes char to scal_lengths_cf_mod.dat
	write(char,file=paste("scal_catch_cf",cf,".dat",sep=""))	

}