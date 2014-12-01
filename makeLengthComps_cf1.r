# -----------------------------------------------------------
# Code to create aggregate proportion at length data from 4 fisheries into 1 (weighed by catch and sample size), and output data in a formate that can be read by ADMB code.
# Author: Jim Richardson

# Set the working directory, (set this to a local folder on your computer)
setwd("/home/jim/Dropbox/REM/Tasks/scal_lengths")

# sources mseRtools.r, needed for lisread function()
source("mseRtools.r")

# sources df_to_dat.r, needed for df_to_dat function
source("df_to_dat.r")

# The following names should be used consistently throughout the various data files:
allNames <- c("RV_4VWX", "HS", "LL_NAFO3_Obs", "LL_NAFO4_Obs", "OT_NAFO3_Obs", "OT_NAFO4_Obs")


# Read the sample data
samples_m <- read.table("NumberMaleMeasured.txt",    header=T)
samples_f <- read.table("NumberFemaleMeasured.txt",  header=T)
samples_c <- read.table("NumberCombinedMeasured.txt",header=T)

# Read in the catch data
landCatchMatrix <- read.table("landCatchMatrix.txt")

# set missing data to NA
landCatchMatrix[landCatchMatrix == -1] <- NA

# lump in OT with other
landCatchMatrix$OT_NAFO3_Obs <- landCatchMatrix$OT_NAFO3_Obs + landCatchMatrix$X3_other
landCatchMatrix$OT_NAFO4_Obs <- landCatchMatrix$OT_NAFO4_Obs + landCatchMatrix$X4_other

# Set the start and end rows in each length prop data file. For example, the full range of length classes is 0-5cm, 5-10cm,... 265-270cm (54 rows). But the data may only include 20cm-210cm (37 rows) since there were no fish measured outside this range.
# The start and end rows are different for each fishery, set them here
firstRow <- rep(NA,6) ; names(firstRow) <- allNames
lastRow <- rep(NA,6) ; names(lastRow) <- allNames

firstRow["RV_4VWX"] <- 2 ; lastRow["RV_4VWX"]  <- 36
firstRow["HS"] <- 5 ; lastRow["HS"] <- 42
firstRow["LL_NAFO3_Obs"] <- 1 ; lastRow["LL_NAFO3_Obs"] <- 54
firstRow["LL_NAFO4_Obs"] <- 1 ; lastRow["LL_NAFO4_Obs"] <- 54
firstRow["OT_NAFO3_Obs"] <- 1 ; lastRow["OT_NAFO3_Obs"] <- 54
firstRow["OT_NAFO4_Obs"] <- 1 ; lastRow["OT_NAFO4_Obs"] <- 54

# Creates empty list, scal_length, will write this to .dat file
scal_lengths <- list()

# Sexes; male, female and combined
sexes <- c("m","f","c")

# Fisheries; commercial fishery, RV survey, halibut survey
# 1 - Commercial fishery
# 2 - RV
# 3 - HS
fisheries <- c("Cf","RV_4VWX","HS")

# define expanded grid (each fishery, sex combination)
grid_FS <- expand.grid(fisheries,sexes)
names(grid_FS) <- c("fishery","sex")

# Check the order is the same in each lengthComps file
for(s in sexes){
    lengths <- lisread(paste("lengthComps_",s,".txt",sep=""))
    for(name in names(lengths)){
        print(paste("sex = ",s,". name = ",name))
    }
}

# Cycles through each fishery, sex combination. If fishery = Cf, aggregates 4 commercial fisheries together
sf <- 1
for(sf in seq(from=1,to=nrow(grid_FS),by=1)){

    # Set sex and fishery, s and f
    s <- grid_FS[sf,"sex"]
    f <- grid_FS[sf,"fishery"]

    # Set name to be written to .dat file
    name_sf <- paste("lenObsProp_",s,f,sep="")

    # Read in prop at length file for sex s.
    lengths <- lisread(paste("lengthComps_",s,".txt",sep=""))

    # Set sample according to sex
    if(s == "m"){
        samples <- samples_m
    } else if (s == "f"){
        samples <- samples_f
    } else if (s == "c"){
        samples <- samples_c
    }

    # cycles through fisheries
    # If RV or HS, no need to aggregate, writes file as is to scal_lengths list
    if(f == "RV_4VWX"){
        scal_lengths[[name_sf]] <- lengths[[paste("RV_4VWX_",s,sep='')]]
    } else if (f == "HS"){
        scal_lengths[[name_sf]] <- lengths[[paste("HS_",s,sep='')]]
    } 
    # if f = Cf, aggregates 4 commercial fisheries into 1.
    else if (f == "Cf"){

        # Settings to format and output aggregated data
        # indices of the fisheries to sum over in aggregating length comps
        commercial_fisheries <- c("LL_NAFO3_Obs", "LL_NAFO4_Obs", "OT_NAFO3_Obs", "OT_NAFO4_Obs")

        # Sets empty prop at length matrix to be filled later
        nlens    <- matrix(0,nrow=54,ncol=nrow(samples))

        # cycles through fisheries that we will aggregate over
        for( comm_f in commercial_fisheries )
        {
            # Get number measured for all years in this fishery
            nMeas <- samples[, comm_f]

            # Get catch data for all years in this fishery
            nCatch <- landCatchMatrix[, comm_f]

            # Get the prop at length data for this fishery
            lens  <- lengths[[paste(comm_f,"_",s,sep='')]]

            # Get the dimensions of the prop at length data (start length class, end length class for this fishery)
            fRow  <- firstRow[comm_f]
            lRow  <- lastRow[comm_f]

            # Set missing variables as NA
            idxNA <- which( lens==-1, arr.ind=T )
            lens[ idxNA ] <- NA

            # multiply props and nMeas to get totals measured
            # by length and accumulate over fisheries
            for( i in 1:ncol(lens) )
            {
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

gi        name_sf

        # adds table to list file
        scal_lengths[[name_sf]] <- as.data.frame(plens)
    }
}

# calculated the proportion female data set for each fishery and stores it in the list
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
3
# nLenBins
54
# firstBin
1 2 5
# lastBin
54 36 42
# binSize
5
# sizeLimit
81 81 81
# lenIndex
1 2 3
# lenLikeWeight
1 1 1
# minLen
11 11 11
# maxLen
276 276 276
# lenFirstYear
19 1 29
# lenLastYear
44 44 44 \n"

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

# writes char to file, scal_lengths_cf1.dat
write(char,file="scal_lengths_cf1.dat")