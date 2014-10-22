# sets working directory
setwd("/home/jim/Dropbox/REM/tasks/scal_lengths")

# imports df_to_dat.r, needed for writing tables
source("df_to_dat.r")

# reads in Number(sex)Measured.txt
NumberMaleMeasured <- read.table("NumberMaleMeasured.txt",header=TRUE)
NumberFemaleMeasured <- read.table("NumberFemaleMeasured.txt",header=TRUE)
NumberCombinedMeasured <- read.table("NumberCombinedMeasured.txt",header=TRUE)

# pulls out names from NumberCombinedMeasured, only one with names for some reason
names <- names(NumberCombinedMeasured)

# reads in Number(sex)Measured_copy.txt
NumberMaleMeasured_copy <- read.table("NumberMaleMeasured copy.txt",header=TRUE)
NumberFemaleMeasured_copy <- read.table("NumberFemaleMeasured copy.txt",header=TRUE)
NumberCombinedMeasured_copy <- read.table("NumberCombinedMeasured copy.txt",header=TRUE)

# replaces RV_4VWX column with column from copied data set
NumberMaleMeasured[,2] <- NumberMaleMeasured_copy[,2]
NumberFemaleMeasured[,2] <- NumberFemaleMeasured_copy[,2]
NumberCombinedMeasured[,2] <- NumberCombinedMeasured_copy[,2]

# sets names for each data set
names(NumberMaleMeasured) <- names
names(NumberFemaleMeasured) <- names
names(NumberCombinedMeasured) <- names

# writes new data sets to .txt file
write(paste(df_to_dat(names),df_to_dat(NumberMaleMeasured),sep=""),'NumberMaleMeasured.txt')
write(paste(df_to_dat(names),df_to_dat(NumberFemaleMeasured),sep=""),'NumberFemaleMeasured.txt')
write(paste(df_to_dat(names),df_to_dat(NumberCombinedMeasured),sep=""),'NumberCombinedMeasured.txt')
