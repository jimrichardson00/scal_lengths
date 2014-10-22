# data cleaning code

# replaced RV_4VWX column in Number(gender)Measured.txt, with RV_4VWX column in Number(gender)Measured_copy.txt
source("replace_RVWX_col.r")

# take the propotion at length data (lengthComps_m etc.) and weigh it by the number measured data (NumberMaleMeasured etc.) for the commercial fisheries (RV and HS are unchanged) and output this into a dat file, scal_lengths.data
source("makeLengthComps.r")

# take scal_lengths.dat and for each year and fishery, aggregate the top and bottom length classes so they are at least 1%
# output the modified file as: scal_lengths_mod.dat
source("lengthComps_cumu.r")