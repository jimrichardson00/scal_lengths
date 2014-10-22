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

write.table(NumberMaleMeasured,'NumberMaleMeasured.txt')
write.table(NumberFemaleMeasured,'NumberFemaleMeasured.txt')
write.table(NumberCombinedMeasured,'NumberCombinedMeasured.txt')




