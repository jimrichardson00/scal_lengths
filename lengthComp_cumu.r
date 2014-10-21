setwd('/home/jim/Dropbox/REM/tasks/scal_lengths')

# sources mseRtools.r, needed for lisread() function
source("mseRtools.r")

# reads in scal_lengths.dat
scal_lengths <- lisread("scal_lengths.dat")

# defined copy of scal_lengths for modification
# will be the same, except will add nObsBins, nstartBins, nendBins, and will change each length prop table
scal_lengths_mod <- scal_lengths

# data.frame to .dat file function
df_to_dat <- function(data){
	
	# start with empty string
	char <- ""

	# for each row, spit out a character with all row elements separated by spaces
	ncols <- length(as.data.frame(data))
	nrows <- nrow(as.data.frame(data))

	# if string is single element, convert it to character
	if(ncols == 1 & nrows == 1){

		char <- paste(as.character(data),"\n",sep="")

	# if string is vector, cycle through elements of vector and output a single line
	} else if (ncols == 1) {

		char_row <- ""
		for(j in seq(from=1,to=nrows,by=1)){
			char_row <- paste(char_row,data[j]," ",sep="")
		}
		char <- paste(char_row,"\n",sep="")

	# if string is data frame, cycle through rows and columns and output a line for each row
	} else if (ncols > 1 & nrows > 1) {

		data <- as.data.frame(data)

		for(i in seq(from=1,to=nrows,by=1)){
			char_row <- ""
			for(j in seq(from=1,to=ncols,by=1)){
				char_row <- paste(char_row,data[i,j]," ",sep="")
			}
			char <- paste(char,char_row,"\n",sep="")
		}
	}

	return(char)
}


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
	scal_lengths_sf <- scal_lengths[paste("lenObsProp_",s,f,sep="")]
	scal_lengths_sf <- as.data.frame(scal_lengths_sf)

	# number of years for current fishery, sex combination (always 44)
	nyears <- length(scal_lengths_sf)

	# number of length classes for current fishery, sex combination (varies by fishery)
	nlength_classes <- nrow(scal_lengths_sf)

	# sets modified scal_lengths_sf data.frame, that we fill in with cumu values
	scal_lengths_sf_mod <- matrix(-1,ncol=nyears,nrow=nlength_classes)
	scal_lengths_sf_mod <- as.data.frame(scal_lengths_sf_mod)

	# firstBin and lastBin
	firstBin <- scal_lengths$firstBin[f]
	lastBin <- scal_lengths$lastBin[f]

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

# names of new scal_lengths file written out in order
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
char <- "## scal_lengths.dat using only 1 fishery and 1 survey \n"
for(name in names){
	char <- paste(char,"# ",name,"\n",df_to_dat(scal_lengths_mod[[name]]),sep="")
}
write(char,file="scal_lengths_mod.dat")

# read in scal_lengths files
# scal_lengths <- lisread("scal_lengths.dat")
# scal_lengths_mod <- lisread("scal_lengths_mod.dat")

# write.csv(scal_lengths$lenObsProp_m2,"lenObsProp_m2.csv")
# write.csv(scal_lengths_mod$lenObsProp_m2,"lenObsProp_m2_mod.csv")

# write.csv(nObsBins,"nObsBins.csv")
# write.csv(startBins,"startBins.csv")
# write.csv(endBins,"endBins.csv")

# write.csv(grid_FS,"grid_FS.csv")