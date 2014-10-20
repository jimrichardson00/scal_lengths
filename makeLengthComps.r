source("mseRtools.r")
# findInterval(x,vec) returns the index of vec containing the closest match to x

# The following names should be used consistently throughout the various data files,
# ADMB, and R code: 
# Year RV_4VWX HS LL_NAFO3_Obs LL_NAFO4_Obs OT_NAFO3_Obs OT_NAFO4_Obs LL_PS OT_PS
# each of the three sample types has it's own lengthComp file with suffix _m, _f,
# or _c
allNames <- c("RV_4VWX", "HS", "LL_NAFO3_Obs", "LL_NAFO4_Obs", "OT_NAFO3_Obs", "OT_NAFO4_Obs")

# The original length comps all had different dimensions, so these give
# the corresponding start rows in the final 54 x 44 matrix. The final matrix
# rows (1:54) correspond to 5-cm length bins seq(from=11,to=276,by=5).
firstRow <- c(2,5,1,1,1,1)
lastRow  <- c(36,42,54,54,54,54)

# Pool length comps across fisheries:
#   i. first re-scale back up to original samples, pal*numberMeasured 
#   ii. sum numbers measured in each length class over fisheries
#  iii. get the total sample size
#   iv. re-normalize length comps back to proportions-at-length

# Read the length comp files
lengths_m <- lisread("lengthComps_m.txt",quiet=T)
lengths_f <- lisread("lengthComps_f.txt",quiet=T)
lengths_c <- lisread("lengthComps_c.txt",quiet=T)

# Read the numbers measured
samples_m <- read.table("NumberMaleMeasured.txt",    header=T)
samples_f <- read.table("NumberFemaleMeasured.txt",  header=T)
samples_c <- read.table("NumberCombinedMeasured.txt",header=T)

# Column names will be used later to grab data
names(samples_m) <- c("Year",allNames)
names(samples_f) <- c("Year",allNames)
names(samples_c) <- c("Year",allNames)

# Settings to format and output aggregated data
# indices of the fisheries to sum over in aggregating length comps
idxFisheries <- c(3,4,6)
samples      <- samples_c      # SET SAMPLES
outFile      <- "plens_c.txt"  # SET OUTPUT FILE

# The ADMB format is m for all fisheries, f for all fisheries, 
# and c for all fisheries
idxNames <- allNames(idxFisheries)
nlens    <- matrix(0,nrow=54,ncol=nrow(samples))
for( f in idxFisheries )
{
    # get number measured for all years in this fishery
    nMeas <- samples[ ,allNames[f] ]
    # get the pal for this fishery: SET LENGTH DATA!!
    lens  <- lengths_c[[f]]
    # get the dimensions
    fRow  <- firstRow[f]
    lRow  <- lastRow[f]

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
# get totals sampled and compared to nMeas - these should match
# if the above code worked

totalSampled <- colSums(nlens,na.rm=T)

# renormalize back to proportions
plens        <- matrix(NA,nrow(nlens),ncol(nlens))
for( i in 1:ncol(nlens) )
  plens[,i] <- nlens[,i]/totalSampled[i]
plens[ is.na(plens) ] <- -1

write( x=t(plens), file=outFile, ncolumns=ncol(nlens) )




