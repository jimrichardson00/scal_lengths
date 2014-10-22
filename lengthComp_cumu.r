setwd('/home/jim/Dropbox/REM/tasks/scal_lengths')

# sources mseRtools.r, needed for lisread() function
source("mseRtools.r")

# reads in scal_lengths2.dat
scal_lengths2 <- lisread("scal_lengths2.dat")

# defined copy of scal_lengths2 for modification
# will be the same, except will add nObsBins, nstartBins, nendBins, and will change each length prop table
scal_lengths_mod <- scal_lengths2

# sources df_to_dat.r, needed for df_to_dat function
source("df_to_dat.r")

# sets theshhold for which values are considered too small
# will sum up values approaching from the top and bottom until threshhold is reached
# current threshhold value is 1% (0.01)
threshhold <- 0.01

# sexes; male, female and combined
sexes <- c("m","f","c")

# fisheries; commercial fishery, RV survey, halibut survey
fisheries <- c(1,2,3)

# define expanded grid (each fishery, sex combination)
grid_FS <- expand.grid(fisheries,sexes)

# empty startBin data.frame
startBins <- matrix(-1,nrow=9,ncol=44)
startBins <- as.data.frame(startBins)

# empy endBin data.frame
endBins <- matrix(-1,nrow=9,ncol=44)
endBins <- as.data.frame(endBins)

# empty nObsBins data.frame
nObsBins <- matrix(-1,nrow=9,ncol=44)
nObsBins <- as.data.frame(nObsBins)

for(sf in seq(from=1,to=nrow(grid_FS),by=1)){

	# current fishery
	f <- grid_FS[sf,1]
	# current sex
	s <- grid_FS[sf,2]

	# pulls data set from .dat file for each fishery, sex combination
	scal_lengths_sf <- scal_lengths2[paste("lenObsProp_",s,f,sep="")]
	scal_lengths_sf <- as.data.frame(scal_lengths_sf)

	# number of years for current fishery, sex combination (always 44)
	nyears <- length(scal_lengths_sf)

	# number of length classes for current fishery, sex combination (varies by fishery)
	nlength_classes <- nrow(scal_lengths_sf)

	# sets modified scal_lengths_sf data.frame, that we fill in with cumu values
	scal_lengths_sf_mod <- matrix(-1,ncol=nyears,nrow=nlength_classes)
	scal_lengths_sf_mod <- as.data.frame(scal_lengths_sf_mod)

	# firstBin and lastBin
	firstBin <- scal_lengths2$firstBin[f]
	lastBin <- scal_lengths2$lastBin[f]

	# lenFirstYear and lenLastYear
	lenFirstYear <- scal_lengths_sf$lenFirstYear[f]
	lenLastYear <- scal_lengths_sf$lenLastYear[f]

	for(y in seq(from=1,to=nyears,by=1)){

		# sets startBin as first bin to get over theshhold (currently 1%) in cumsum
		# -1 if bins never get over 1% in cumsum
		startBin <- sum(cumsum(scal_lengths_sf[,y]) < threshhold) + 1
		startBin <- ifelse(startBin > nlength_classes,-1,startBin)
		# sets value at startBin, i.e. cumu value in first n bins
		# -1 if bins never get over theshhold in cumsum
		startBin_value <- cumsum(scal_lengths_sf[,y])[startBin]
		startBin_value <- ifelse(startBin > nlength_classes,-1,startBin_value)

		#sets endBin as first bin to get over theshhold (currently 1%) in cumsum (working backwards)
		endBin <- nlength_classes - sum(cumsum(scal_lengths_sf[,y][seq(from=nlength_classes,to=1,by=-1)]) < threshhold)
		endBin <- ifelse(endBin == 0 | endBin > nlength_classes,-1,endBin)
		# sets value at endBin, i.e. cumu value in first n bins
		endBin_value <- cumsum(scal_lengths_sf[,y][seq(from=nlength_classes,to=1,by=-1)])[nlength_classes - endBin + 1]
		endBin_value <- ifelse(endBin > nlength_classes,-1,endBin_value)

		# sets nObsBin as the number of bins with observations in them
		# if either startBin or endBin is -1, sets noObsBin as -1
		nObsBin <- endBin - startBin + 1
		nObsBin <- ifelse(endBin == -1 | startBin == -1, -1,nObsBin)

		# for each year & fishery-sex combination, writes startBin value to matrix startBins
		startBins[sf,y] <- startBin

		# for each year & fishery-sex combination, writes endBin value to matrix endBin
		endBins[sf,y] <- endBin

		# for each year & fishery-sex combination, writes nObsBin value to matrix endBin
		nObsBins[sf,y] <- nObsBin

		# writes new scal_lengths_sf_mod file, with values summed up on top and bottom
		# cycles through rows (currently cycling through year and fishery-sex data set)
		for(l in seq(from=1,to=nlength_classes,by=1)){
			# if line is startBin, sets as the sum of the first 1-startBin values
			if(l == startBin){
				scal_lengths_sf_mod[l,y] <- startBin_value
			} 
			# if line is endBin, sets as sum of the last nlength_classes - endBin values
			else if (l == endBin){
				scal_lengths_sf_mod[l,y] <- endBin_value
			}
			# if startBin < l < endBin, sets as value from .dat file
			else if (l > startBin & l < endBin){
				scal_lengths_sf_mod[l,y] <- scal_lengths_sf[l,y]
			}
			# otherwise does not modify scal_lengths_sf_mod, which is all -1
		}
	}

	# updates the lensProp tables with cumu values
	scal_lengths_mod[[paste("lenObsProp_",s,f,sep="")]] <- scal_lengths_sf_mod

}

# adds the nObsBins, startBins and endBins tables to the modified list
scal_lengths_mod$nObsBins <- nObsBins
scal_lengths_mod$startBins <- startBins
scal_lengths_mod$endBins <- endBins

# names of new scal_lengths2 file written out in order
names <- c("nLenSeries", 
"nLenBins", 
"firstBin", 
"lastBin", 
"binSize", 
"sizeLimit", 
"lenIndex", 
"lenLikeWeight", 
"minLen", 
"maxLen", 
"lenFirstYear", 
"lenLastYear", 

# adding startBins
"startBins",
# adding endBins
"endBins",
# adding nObsBins
"nObsBins",

"lenObsProp_m1", 
"lenObsProp_m2", 
"lenObsProp_m3",
"lenObsProp_f1", 
"lenObsProp_f2", 
"lenObsProp_f3", 
"lenObsProp_c1", 
"lenObsProp_c2", 
"lenObsProp_c3")

# write new scal_lengths_mod file to .dat file
char <- "## scal_lengths2.dat using only 1 fishery and 1 survey \n"

# add data set from each fishery-sex combo
for(name in names){
	char <- paste(char,"# ",name,"\n",df_to_dat(scal_lengths_mod[[name]]),sep="")
}

# adds proportion female for length class and year, by fishery to char for writing
# cycles through fisheries
for(f in seq(from=1,to=3,by=1)){

    # sets name to obtain prop female table for fishery
    name_f <- paste("PropFemale",f,sep="")

    # pasters prop female table into char
    char <- paste(char,"# ",name_f,"\n",df_to_dat(scal_lengths2[[name_f]]),sep="")

}

# writes char to scal_lengths_mod.dat
write(char,file="scal_lengths_mod.dat")

# read in scal_lengths2 files
# scal_lengths <- lisread("scal_lengths.dat")
scal_lengths_mod <- lisread("scal_lengths_mod.dat")

# write.csv(scal_lengths2$lenObsProp_f2,"lenObsProp_f2.csv")
# write.csv(scal_lengths_mod$lenObsProp_f2,"lenObsProp_f2_mod.csv")

# write.csv(nObsBins,"nObsBins.csv")
# write.csv(startBins,"startBins.csv")
# write.csv(endBins,"endBins.csv")

# write.csv(grid_FS,"grid_FS.csv")