setwd("/home/jim/Dropbox/REM/tasks/scal_lengths")

source("df_to_dat.r")

source("mseRtools.r")

NumberMaleMeasured <- read.table("NumberMaleMeasured.txt",header=TRUE)
NumberFemaleMeasured <- read.table("NumberFemaleMeasured.txt",header=TRUE)
NumberCombinedMeasured <- read.table("NumberCombinedMeasured.txt",header=TRUE)

NumberMaleMeasured_copy <- read.table("NumberMaleMeasured copy.txt",header=TRUE)
NumberFemaleMeasured_copy <- read.table("NumberFemaleMeasured copy.txt",header=TRUE)
NumberCombinedMeasured_copy <- read.table("NumberCombinedMeasured copy.txt",header=TRUE)

NumberMaleMeasured[,2] <- NumberMaleMeasured_copy[,2]
NumberFemaleMeasured[,2] <- NumberFemaleMeasured_copy[,2]
NumberCombinedMeasured[,2] <- NumberCombinedMeasured_copy[,2]

allNames <- c("RV_4VWX", "HS", "LL_NAFO3_Obs", "LL_NAFO4_Obs", "OT_NAFO3_Obs", "OT_NAFO4_Obs")
names(NumberMaleMeasured) <- allNames
names(NumberFemaleMeasured) <- allNames
names(NumberCombinedMeasured) <- allNames

write(paste(df_to_dat(c("Year",allNames)),df_to_dat(NumberMaleMeasured),sep=""),'NumberMaleMeasured.txt')
write(paste(df_to_dat(c("Year",allNames)),df_to_dat(NumberFemaleMeasured),sep=""),'NumberFemaleMeasured.txt')
write(paste(df_to_dat(c("Year",allNames)),df_to_dat(NumberCombinedMeasured),sep=""),'NumberCombinedMeasured.txt')
