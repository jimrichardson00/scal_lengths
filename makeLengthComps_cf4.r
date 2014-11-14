# sets working directory
# setwd("/home/jim/Dropbox/REM/tasks/scal_lengths")

# sources mseRtools.r, needed for lisread function()
source("mseRtools.r")
# findInterval(x,vec) returns the index of vec containing the closest match to x

# sources df_to_dat.r, needed for df_to_dat function
source("df_to_dat.r")

# The following names should be used consistently throughout the various data files,
# ADMB, and R code: 
# Year RV_4VWX HS LL_NAFO3_Obs LL_NAFO4_Obs OT_NAFO3_Obs OT_NAFO4_Obs LL_PS OT_PS
# each of the three sample types has it's own lengthComp file with suffix _m, _f,
# or _c
allNames <- c("RV_4VWX", "HS", "LL_NAFO3_Obs", "LL_NAFO4_Obs", "OT_NAFO3_Obs", "OT_NAFO4_Obs")

# Read the sample data
samples_m <- read.table("NumberMaleMeasured.txt",    header=T)
samples_f <- read.table("NumberFemaleMeasured.txt",  header=T)
samples_c <- read.table("NumberCombinedMeasured.txt",header=T)

# Column names will be used later to grab data
names(samples_m) <- c("Year",allNames)
names(samples_f) <- c("Year",allNames)
names(samples_c) <- c("Year",allNames)

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

# The original length comps all had different dimensions, so these give
# the corresponding start rows in the final 54 x 44 matrix
firstRow <- c(2,5,1,1,1,1)
lastRow  <- c(36,42,54,54,54,54)

# creates empty list, scal_length, will write this to .dat file
scal_lengths <- list()

# sexes; male, female and combined
sexes <- c("m","f","c")

# fisheries; commercial fishery, RV survey, halibut survey
# 1 - LL - 03 - Commercial fishery
# 2 - LL - 04 - Commercial fishery
# 3 - OT - 03 - Commercial fishery
# 4 - OT - 04 - Commercial fishery
# 5 - RV
# 6 - HS
fisheries <- c(1,2,3,4,5,6)

# commercial fisheries
commercial_fisheries <- c("LL","OT")

# define expanded grid (each fishery, sex combination)
grid_FS <- expand.grid(fisheries,sexes)
names(grid_FS) <- c("fishery","sex")

# check the order is the same in each lengthComps file
for(s in sexes){
    lengths <- lisread(paste("lengthComps_",s,".txt",sep=""))
    for(name in names(lengths)){
        print(paste("sex = ",s,". name = ",name))
    }
}

sf <- 1
for(sf in seq(from=1,to=nrow(grid_FS),by=1)){

    s <- grid_FS[sf,"sex"]
    f <- grid_FS[sf,"fishery"]

    name_sf <- paste("lenObsProp_",s,f,sep="")

    lengths <- lisread(paste("lengthComps_",s,".txt",sep=""))

    # set sample according to sex
    if(s == "m"){
        samples <- samples_m
    } else if (s == "f"){
        samples <- samples_f
    } else if (s == "c"){
        samples <- samples_c
    }

    # cycles through fisheries
    if(f == 5){
        scal_lengths[[name_sf]] <- lengths[[1]]
    } else if (f == 6){
        scal_lengths[[name_sf]] <- lengths[[2]]
    } else if (f == 1 | f == 2 | f == 3 | f == 4){

        # Pool length comps across fisheries:
        #   i. first re-scale back up to original samples, pal*numberMeasured 
        #   ii. sum numbers measured in each length class over fisheries
        #  iii. get the total sample size
        #   iv. re-normalize length comps back to proportions-at-length

        # LL & 03
        if(f == 1){
            idxFisheries <- c(3)
        # LL & 04
        } else if (f == 2){
            idxFisheries <- c(4)
        } else if (f == 3){
            idxFisheries <- c(5)
        } else if (f == 4){
            idxFisheries <- c(6)
        }

        # The ADMB format is m for all fisheries, f for all fisheries, 
        # and c for all fisheries
        idxNames <- allNames(idxFisheries)
        nlens    <- matrix(0,nrow=54,ncol=nrow(samples))

        for( idxF in idxFisheries ){
            # get number measured for all years in this fishery
            nMeas <- samples[ ,allNames[idxF] ]

            # get catch data for all years in this fishery
            nCatch <- landCatchMatrix[,allNames[idxF]]

            # get the pal for this fishery: SET LENGTH DATA!!
            lens  <- lengths[[idxF]]
            # get the dimensions
            fRow  <- firstRow[idxF]
            lRow  <- lastRow[idxF]

            # find the row/col indices of missing values
            idxNA <- which( lens==-1, arr.ind=T )
            # convert -1 to NA
            lens[ idxNA ] <- NA
            #browser()
            # multiply props and nMeas to get totals measured
            # by length and accumulate over fisheries
            for( i in 1:ncol(lens) ){
                #if(f==4 & nMeas[i]>0) browser()
                nlens[fRow:lRow,i] <- nlens[fRow:lRow,i] + nCatch[i]*nMeas[i]*lens[,i]
            }
        }

        # adds up colums so we can renormalize
        totalSampled <- colSums(nlens,na.rm=T)

        # renormalize back to proportions
        plens        <- matrix(NA,nrow(nlens),ncol(nlens))
        for( i in 1:ncol(nlens) ){
            plens[,i] <- nlens[,i]/totalSampled[i]
        }
        # replaces missing data with -1
        plens[ is.na(plens)] <- -1

        # adds table to list file
        scal_lengths[[name_sf]] <- as.data.frame(plens)
    }
}


# calculated the proportion female data set for each fishery and stores it in the list
f <- 1
for(f in seq(from=1,to=length(fisheries),by=1)){

    # names the proportion at length, year for males in data set
    lenObsProp_m <- scal_lengths[[paste("lenObsProp_m",f,sep="")]]
    lenObsProp_m[lenObsProp_m == -1] <- NA
    # names the proportion at length, year for females in data set
    lenObsProp_f <- scal_lengths[[paste("lenObsProp_f",f,sep="")]]
    lenObsProp_f[lenObsProp_f == -1] <- NA

    # creates a new table called PropFemale, which measures the proportion of females in a given length class, year, fishery
    PropFemale <- lenObsProp_f/(lenObsProp_m + lenObsProp_f)
    # replaces NA with -1, these are entries where prop_m & prop_f are zero
    # (no fish in this length class, year, so data is missing)
    PropFemale[is.na(PropFemale)] <- -1
    scal_lengths[[paste("PropFemale",f,sep="")]] <- PropFemale
}

# write new scal_lengths_mod file to .dat file
char <- "## scal_lengths.dat using only 1 fishery and 1 survey
# nLenSeries
6
# nLenBins
54
# firstBin
1 1 1 1 2 5
# lastBin
54 54 54 54 36 42
# binSize
5
# sizeLimit
81 81 81 81 81 81
# lenIndex
1 2 3 4 5 6
# lenLikeWeight
1 1 1 1 1 1
# minLen
11 11 11 11 11 11
# maxLen
276 276 276 276 276 176
# lenFirstYear
19 19 18 18 1 29
# lenLastYear
44 44 44 44 44 44 \n"

# adds proportion for length class & year by sex, fishery to char for writing
# cycles through each sex-fishery combination
for(sf in seq(from=1,to=nrow(grid_FS),by=1)){

    # sets the sex s and fishery f from grid
    s <- grid_FS[sf,"sex"]
    f <- grid_FS[sf,"fishery"]

    # sets name to obtain table for sex, fishery combo
    name_sf <- paste("lenObsProp_",s,f,sep="")

    # pastes table for sex, fishery combo into char
    char <- paste(char,"# ",name_sf,"\n",df_to_dat(scal_lengths[[name_sf]]),sep="")

}

# adds proportion female for length class and year, by fishery to char for writing
# cycles through fisheries
for(f in seq(from=1,to=length(fisheries),by=1)){

    # sets name to obtain prop female table for fishery
    name_f <- paste("PropFemale",f,sep="")

    # pasters prop female table into char
    char <- paste(char,"# ",name_f,"\n",df_to_dat(scal_lengths[[name_f]]),sep="")

}

# writes char to file, scal_lengths2.dat
write(char,file="scal_lengths_cf4.dat")