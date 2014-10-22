# data.frame to .dat file function
df_to_dat <- function(data){

	# convert data to data.frame, if something else (matrix?)
	data <- as.data.frame(data)
	
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
