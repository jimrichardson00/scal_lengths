# data cleaning code

# sets working directory
# setwd("/home/jim/Dropbox/REM/tasks/scal_lengths")

# replaced RV_4VWX column in Number(gender)Measured.txt, with RV_4VWX column in Number(gender)Measured_copy.txt
source("replace_RV_4VWX_col.r")

# take the propotion at length data (lengthComps_m etc.) and weigh it by the number measured data (NumberMaleMeasured etc.) for the commercial fisheries (RV and HS are unchanged) and output this into a dat file, scal_lengths2.dat
# outputs to scal_lenghts2.dat

# cf1
source("makeLengthComps_cf1.r")

# cf2
source("makeLengthComps_cf2.r")

# cf4
source("makeLengthComps_cf4.r")

# take output from makeLengthComps data and for each year and fishery, aggregate the top and bottom length classes so they are at least 1%
# output the modified file in each case as: scal_lengths_cf_mod.dat
source("lengthComp_cumu.r")