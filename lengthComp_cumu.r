# sets working directory
setwd('/home/jim/Dropbox/REM/tasks/scal_lengths')

# sources mseRtools.r, needed for lisread() function
source("mseRtools.r")

# sources df_to_dat.r, needed for df_to_dat function
source("df_to_dat.r")

# sets theshhold for which values are considered too small
# will sum up values approaching from the top and bottom until threshhold is reached
# current threshhold value is 1% (0.01)
threshhold <- 0.01


# # check that each .dat file has all the relevant matricies
# for(cf in c(1,2,4)){

# 	scal_lengths_cf <- lisread(paste("scal_lengths_cf",cf,".dat",sep=""))

# 	print("----------------------------------")
# 	print(paste("cf = ",cf,sep=""))
# 	print("")

# 	for(name in names(scal_lengths_cf)){
# 		print(paste(name))
# 		print(paste("number rows ",nrow(scal_lengths_cf[[name]])))
# 		print(paste("number cols ",ncol(scal_lengths_cf[[name]])))
# 		print(paste("length ",length(scal_lengths_cf[[name]])))
# 		print("")
# 	}
# }

for(cf in c(1,2,4)){

	# reads in scal_lengths_cf.dat
	scal_lengths_cf <- lisread(paste("scal_lengths_cf",cf,".dat",sep=""))

	# defined copy of scal_lengths for modification
	# will be the same, except will add nObsBins, nstartBins, nendBins, and will change each length prop table
	scal_lengths_cf_mod <- scal_lengths_cf

	# sexes; male, female and combined
	sexes <- c("m","f","c")

	# fisheries; commercial fishery, RV survey, halibut survey
	if (cf == 1) {
		# 1 - Commercial fishery
		# 2 - RV
		# 3 - HS
		fisheries <- c(1,2,3)
	} else if (cf == 2) {
		# fisheries; commercial fishery, RV survey, halibut survey
		# 1 - LL - Commercial fishery
		# 2 - OT - Commercial fishery
		# 3 - RV
		# 4 - HS
		fisheries <- c(1,2,3,4)
	} else if (cf == 4) {
		# 1 - LL - 03 - Commercial fishery
		# 2 - LL - 04 - Commercial fishery
		# 3 - OT - 03 - Commercial fishery
		# 4 - OT - 04 - Commercial fishery
		# 5 - RV
		# 6 - HS
		fisheries <- c(1,2,3,4,5,6)
	}

	# define expanded grid (each fishery, sex combination)
	grid_FS <- expand.grid(fisheries,sexes)
	names(grid_FS) <- c("fishery","sex")

	# empty startBin data.frame
	startBins <- matrix(-1,nrow=nrow(grid_FS),ncol=44)
	startBins <- as.data.frame(startBins)

	# empy endBin data.frame
	endBins <- matrix(-1,nrow=nrow(grid_FS),ncol=44)
	endBins <- as.data.frame(endBins)

	# empty nObsBins data.frame
	nObsBins <- matrix(-1,nrow=nrow(grid_FS),ncol=44)
	nObsBins <- as.data.frame(nObsBins)

	for(sf in seq(from=1,to=nrow(grid_FS),by=1)){

		# current sex
		s <- grid_FS[sf,"sex"]
		# current fishery
		f <- grid_FS[sf,"fishery"]

		# pulls data set from .dat file for each fishery, sex combination
		scal_lengths_sf <- scal_lengths_cf[paste("lenObsProp_",s,f,sep="")]
		scal_lengths_sf <- as.data.frame(scal_lengths_sf)

		# number of years for current fishery, sex combination (always 44)
		nyears <- length(scal_lengths_sf)

		# number of length classes for current fishery, sex combination (varies by fishery)
		nlength_classes <- nrow(scal_lengths_sf)

		# sets modified scal_lengths_sf data.frame, that we fill in with cumu values
		scal_lengths_sf_mod <- matrix(-1,ncol=nyears,nrow=nlength_classes)
		scal_lengths_sf_mod <- as.data.frame(scal_lengths_sf_mod)

		# firstBin and lastBin
		firstBin <- scal_lengths_cf$firstBin[f]
		lastBin <- scal_lengths_cf$lastBin[f]

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
		scal_lengths_cf_mod[[paste("lenObsProp_",s,f,sep="")]] <- scal_lengths_sf_mod
	}

	# adds the nObsBins, startBins and endBins tables to the modified list
	scal_lengths_cf_mod$nObsBins <- nObsBins
	scal_lengths_cf_mod$startBins <- startBins
	scal_lengths_cf_mod$endBins <- endBins

	# names of new scal_lengths_cf file written out in order
	names <- c(names(scal_lengths_cf)[seq(1,12,1)],
		c("startBins","endBins","nObsBins"),
		names(scal_lengths_cf)[seq(13,length(fisheries)*3 + 12,1)])

	# write new scal_lengths_cf_mod file to .dat file
	char <- "## scal_lengths_cf.dat using only 1 fishery and 1 survey \n"

	# add data set from each fishery-sex combo
	for(name in names){
		char <- paste(char,"# ",name,"\n",df_to_dat(scal_lengths_cf_mod[[name]]),sep="")
	}

	# adds proportion female for length class and year, by fishery to char for writing
	# cycles through fisheries
	for(f in seq(from=1,to=length(fisheries),by=1)){

	    # sets name to obtain prop female table for fishery
	    name_f <- paste("PropFemale",f,sep="")

	    # pasters prop female table into char
	    char <- paste(char,"# ",name_f,"\n",df_to_dat(scal_lengths_cf[[name_f]]),sep="")

	}

	# writes char to scal_lengths_cf_mod.dat
	write(char,file=paste("scal_lengths_cf",cf,"_mod.dat",sep=""))

}
