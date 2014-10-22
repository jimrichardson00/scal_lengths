scal_lengths / makeLengthComps.r
jimrichardson00jimrichardson00 37 minutes ago new branch, modifying makeLengthComps.r to include table that contain…
1 contributor
231 lines (181 sloc)  6.453 kb RawBlameHistory  
# sets working directory
setwd("/home/jim/Dropbox/REM/tasks/scal_lengths")

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

# The original length comps all had different dimensions, so these give
# the corresponding start rows in the final 54 x 44 matrix
firstRow <- c(2,5,1,1,1,1)
lastRow  <- c(36,42,54,54,54,54)

# Read the numbers measured
samples_m <- read.table("NumberMaleMeasured.txt",    header=T)
samples_f <- read.table("NumberFemaleMeasured.txt",  header=T)
samples_c <- read.table("NumberCombinedMeasured.txt",header=T)

# Column names will be used later to grab data
names(samples_m) <- c("Year",allNames)
names(samples_f) <- c("Year",allNames)
names(samples_c) <- c("Year",allNames)

# The original length comps all had different dimensions, so these give
# the corresponding start rows in the final 54 x 44 matrix
firstRow <- c(2,5,1,1,1,1)
lastRow  <- c(36,42,54,54,54,54)

# creates empty list, scal_length, will write this to .dat file
scal_lengths <- list()

# sexes; male, female and combined
sexes <- c("m","f","c")

# fisheries; all
fisheries <- allNames

# years: from 1 to 44 (1970-2013)
years <- seq(from=1,to=44,by=1)

# length classes: from 1-54 (varies by fishery)
lengths <- seq(from=1,to=54,by=1)

grid <- expand.grid(fisheries,sexes,years,lengths)
grid <- as.data.frame(grid)
names(grid) <- c("fisheries","sexes","years","lengths")

grid$Prop_fsyl <- rep(NA,nrow(grid))
grid$Total_fsy <- rep(NA,nrow(grid))

nrow(grid)

v <- 1
l <- 36
for(v in seq(from=1,to=nrow(grid),by=1)){

    f <- as.character(grid[v,"fisheries"])
    s <- as.character(grid[v,"sexes"])
    y <- grid[v,"years"]
    l <- grid[v,"lengths"]

    # set sample
    if(s == "m"){
        samples <- samples_m
    } else if (s == "f"){
        samples <- samples_f
    } else if (s == "c"){
        samples <- samples_c
    }

    f_name <- grep(f,names(lisread(paste("lengthComps_",s,".txt",sep=""))),value=TRUE)

    lenComps_fs <- lisread(paste("lengthComps_",s,".txt",sep=""))[[f_name]]


    if(l > nrow(lenComps_fs) | y > ncol(lenComps_fs)){
        grid[v,"Prop_fsyl"] <- NA
        grid[v,"Total_fsy"] <- NA
    } else {
        grid[v,"Prop_fsyl"] <- lenComps_fs[l,y]
        grid[v,"Total_fsy"] <- samples[y,f]
    }

}

# /////////////////////////////////////////////////


# define expanded grid (each fishery, sex combination)
grid_FS <- expand.grid(fisheries,sexes)
names(grid_FS) <- c("fishery","sex")

for(sf in seq(from=1,to=nrow(grid_FS),by=1)){

    s <- grid_FS[sf,"sex"]
    f <- grid_FS[sf,"fishery"]

    name_sf <- paste("lenObsProp_",s,f,sep="")

    lengths <- lisread(paste("lengthComps_",s,".txt",sep=""))

    # set sample
    if(s == "m"){
        samples <- samples_m
    }
    else if (s == "f"){
        samples <- samples_f
    }
    else if (s == "c"){
        samples <- samples_c
    }

    # cycles through fisheries
    if(f == 2){
        scal_lengths[[name_sf]] <- lengths[[1]]
    }
    else if (f == 3){
        scal_lengths[[name_sf]] <- lengths[[2]]
    }
    else if (f == 1){

        # Pool length comps across fisheries:
        #   i. first re-scale back up to original samples, pal*numberMeasured 
        #   ii. sum numbers measured in each length class over fisheries
        #  iii. get the total sample size
        #   iv. re-normalize length comps back to proportions-at-length

        # Settings to format and output aggregated data
        # indices of the fisheries to sum over in aggregating length comps
        idxFisheries <- c(3,4,6)

        # The ADMB format is m for all fisheries, f for all fisheries, 
        # and c for all fisheries
        idxNames <- allNames(idxFisheries)
        nlens    <- matrix(0,nrow=54,ncol=nrow(samples))
        for( idxF in idxFisheries )
        {
            # get number measured for all years in this fishery
            nMeas <- samples[ ,allNames[idxF] ]
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
            for( i in 1:ncol(lens) )
            {
              #if(f==4 & nMeas[i]>0) browser()
              nlens[fRow:lRow,i] <- nlens[fRow:lRow,i] + nMeas[i]*lens[,i]
            }
        }

        totalSampled <- colSums(nlens,na.rm=T)

        # renormalize back to proportions
        plens        <- matrix(NA,nrow(nlens),ncol(nlens))
        for( i in 1:ncol(nlens) ){
            plens[,i] <- nlens[,i]/totalSampled[i]
        }
        plens[ is.na(plens) ] <- -1

        scal_lengths[[name_sf]] <- as.data.frame(plens)

        scal _lengths[[paste("nlens_",s,)]] nlens

    }
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

for(sf in seq(from=1,to=nrow(grid_FS),by=1)){

    s <- grid_FS[sf,"sex"]
    f <- grid_FS[sf,"fishery"]

    name_sf <- paste("lenObsProp_",s,f,sep="")

    char <- paste(char,"# ",name_sf,"\n",df_to_dat(scal_lengths[[name_sf]]),sep="")

}

write(char,file="scal_lengths2.dat")

scal_lengths <- lisread("scal_lengths.dat")
scal_lengths2 <- lisread("scal_lengths2.dat")

# check that the code is doing the right thing
for(name in names(scal_lengths)){
    print(paste(name," difference = ",norm(as.matrix(scal_lengths[[name]]) - as.matrix(scal_lengths2[[name]]),type="f"),sep=""))
}