setwd("/home/jim/Dropbox/REM/tasks/scal_lengths")

source("mseRtools.r")
source("scal_globals.r")
#source("mcmcpairs.r")
require("MCMCpack")

# This rep file loading can be improved by (i) adding some time/date
# info to the first line of the rep file, (ii) reading the first line
# using scan(), (iii) letting the scal_plots user know when the rep 
# file was created and perhaps whether it is a converged fit.
repFileName <- "scal.rep"
if( file.exists(repFileName) )
{
  repObj <- lisread(repFileName)
  cat( paste(repFileName," successfully imported...","\n",sep="") )  
}
if( !file.exists(repFileName) )
{
  cat("No .rep file to import...","\n")  
}

# read in scal_lengths_mods
scal_lengths_mod <- lisread("scal_lengths_mod.dat")

# plotLenFit
# Plots observed length-frequency (bars) and fitted model length-frequency
# for each year. A separate page is generated for each fishery/survey
plotLenFit_PropFemale <- function( obj, modelCheck=F, doType="f", resGear=1, png=F )
{
  # the model uses pl_m(g)(lengths)(t)
  # data are lenObsProp_m(g)

  firstBin <- obj$firstBin
  lastBin  <- obj$lastBin
  minLen   <- obj$minLen
  maxLen   <- obj$maxLen
  nGear    <- 3
  nLens    <- 54 #nrow( pl_m[[1]] )
  nT       <- 44 #ncol( pl_m[[1]] )
  lenFirstYear <- obj$lenFirstYear
  lenLastYear  <- obj$lenLastYear

  # Model: pl_m,f,c are lists of arrays with dimensions( nLens,nT ).
  pFemale_g      <- vector("list")

  pFemale_g[[1]] <- matrix(NA, nrow=nLens,ncol=nT)
  # pFemale_g[[1]] <- switch(doType,
  #                   m = obj$pl_m1,
  #                   f = obj$pl_f1,
  #                   c = obj$pl_c1
  #                   )
  pFemale_g[[1]] <- repObj$pFemale_1

  pFemale_g[[2]] <- matrix(NA, nrow=nLens,ncol=nT)
  # pFemale_g[[2]] <- switch(doType,
  #                   m = obj$pl_m2,
  #                   f = obj$pl_f2,
  #                   c = obj$pl_c2
  #                   )
  pFemale_g[[2]] <- repObj$pFemale_2

  pFemale_g[[3]] <- matrix(NA, nrow=nLens,ncol=nT)
  # pFemale_g[[3]] <- switch(doType,
  #                   m = obj$pl_m3,
  #                   f = obj$pl_f3,
  #                   c = obj$pl_c3
  #                   )
  pFemale_g[[3]] <- repObj$pFemale_3
  
  lenClasses <- vector("list")
  lenClasses[[1]] <- seq( minLen[1], maxLen[1], by=obj$binSize )
  lenClasses[[2]] <- seq( minLen[2], maxLen[2], by=obj$binSize )
  lenClasses[[3]] <- seq( minLen[3], maxLen[3], by=obj$binSize )

  # These are the lengths.
  xAxisLengths    <- obj$lfBins
  xLim <- range( xAxisLengths )
  
  # OBS: lenObsProp is a 3-dim array with dimensions( nGear,nLens,nT ).
  lenObsProp      <- array( 0,dim=c(nGear,nLens,nT) )
  # lenObsProp[1,firstBin[1]:lastBin[1],] <- switch(doType,
  #                                                 m = obj$lenObsProp_m1,
  #                                                 f = obj$lenObsProp_f1,
  #                                                 c = obj$lenObsProp_c1
  #                                                 )
  lenObsProp[1,firstBin[1]:lastBin[1],] <- scal_lengths_mod$PropFemale1

  # lenObsProp[2,firstBin[2]:lastBin[2],] <- switch(doType,
  #                                                 m = obj$lenObsProp_m2,
  #                                                 f = obj$lenObsProp_f2,
  #                                                 c = obj$lenObsProp_c2
  #                                                 )
  lenObsProp[2,firstBin[2]:lastBin[2],] <- scal_lengths_mod$PropFemale2

  # lenObsProp[3,firstBin[3]:lastBin[3],] <- switch(doType,
  #                                                 m = obj$lenObsProp_m3,
  #                                                 f = obj$lenObsProp_f3,
  #                                                 c = obj$lenObsProp_c3
  #                                                 )
  lenObsProp[3,firstBin[3]:lastBin[3],] <- scal_lengths_mod$PropFemale3

  # Starting Bins: only using length bins where lenObsProp > 0.01. All 
  # entries < 0.01 are accumulated into start and end bin classes
  startBins <- matrix( 0,nrow=nGear,ncol=nT )
  startBins <- switch(doType,
                      m = obj$startBin_m,
                      f = obj$startBin_f,
                      c = obj$startBin_c
                      )
  # Ending Bins: only using length bins where lenObsProp > 0.01. All 
  # entries < 0.01 are accumulated into start and end bin classes
  endBins <- matrix( 0,nrow=nGear,ncol=nT )
  endBins <- switch(doType,
                    m = obj$endBin_m,
                    f = obj$endBin_f,
                    c = obj$endBin_c
                    )
  # Starting Bins: only using length bins where lenObsProp > 0.01. All 
  # entries < 0.01 are accumulated into start and end bin classes
  nObsBins <- matrix( 0,nrow=nGear,ncol=nT )
  nObsBins <- switch(doType,
                      m = obj$nObsBins_m,
                      f = obj$nObsBins_f,
                      c = obj$nObsBins_c
                      )
  # Get those -1 values outta there for plotting...convert to NA
  setMissing <- function( vec )
  {
    n <- length(vec)
    vec[ vec<0 ] <- NA
    return(vec)
  }
  for( i in 1:dim(lenObsProp)[1] )
    lenObsProp[i,,] <- apply( X=lenObsProp[i,,], MARGIN=2,FUN=setMissing )
 
 
  # These are the fitted len proportions.
  yLim <- NULL
  if ( is.null(yLim) )
    yLim <- c( 0,1.1 )

  for ( g in 1:nGear )
  {   

    if ( g > 1 )
    {
      dev.new()
    }
    
    yLim    <- c( 0, 1.1 )
    sumProp <- colSums( lenObsProp[g,,], na.rm=TRUE )      
    nPanels <- obj$lenLastYear[g] - obj$lenFirstYear[g] + 1

    myMar   <- c( 1.5,1.25,0.5,0.5 )
    par( oma=.OMA, mar=myMar, mfcol=c(5,5) )

    # Years where the sum is not 0, i.e., there are fitted lens.
    idxPlot <- c(1:nT)[sumProp>0.0]
    idxPlot <- seq(obj$lenFirstYear[g],obj$lenLastYear[g],by=1) 
    idxPlot <- idxPlot[sumProp[idxPlot]!=0.0]
    
    # List to hold the residuals
    if( g==resGear )
    {
      resids      <- vector("list",length=length(idxPlot)+3)
      cumulObs    <- rep(0,length(lenClasses[[g]]))
      cumulPreds  <- rep(0,length(lenClasses[[g]]))
      cumulResids <- rep(0,length(lenClasses[[g]]))
    }
    delta <- 1
    iPlot <- 0
    for ( i in idxPlot )
    {
      if( g==2 & i==26 )
      {
        dev.new()
        myMar   <- c( 1.5,1.25,0.5,0.5 )
        par( oma=.OMA, mar=myMar, mfcol=c(5,5) )
      }

      plot( xLim, yLim, type="n", axes=FALSE, xlab="", ylab="" )
      
      abline( h=0 )

      # The observed length proportions.
      rect( lenClasses[[g]]-delta, 0, lenClasses[[g]]+delta, 
            lenObsProp[g,,i], col="white",lwd=0.5 )
      
      # The fitted length proportions.      
      lines( lenClasses[[g]], pFemale_g[[g]][,i], lty=1, lwd=2, col="red" )
      points( lenClasses[[g]], pFemale_g[[g]][,i], bg="white", 
              col="black", cex=0.5, lwd=1, pch=21 )
      
      if( g==resGear )
      {
        iPlot <- iPlot + 1
        resids[[iPlot]]$points    <- vector("numeric",length(lenClasses[[g]]))  
        resids[[iPlot]]$stdResids <- vector("numeric",length(lenClasses[[g]]))  
        resids[[iPlot]]$std       <- vector("numeric",length(lenClasses[[g]]))  
        
        resids[[iPlot]]$points    <- lenObsProp[g,,i] - pFemale_g[[g]][,i]
        resids[[iPlot]]$std       <- sd( resids[[iPlot]]$points, na.rm=TRUE )  
        resids[[iPlot]]$stdResids <- resids[[iPlot]]$points/resids[[iPlot]]$std
        resids[[iPlot]]$stdResids[is.na(resids[[iPlot]]$stdResids)] <- 0
        resids[[iPlot]]$year    <- paste( i+.INITYEAR-1,",",i,sep="") 

        tmp <- lenObsProp[g,,i]
        tmp[is.na(tmp)] <- 0.
        cumulObs <- cumulObs + tmp
        
        tmp <- pFemale_g[[g]][,i]
        tmp[is.na(tmp)] <- 0.
        cumulPreds <- cumulPreds + tmp

        tmp <- resids[[iPlot]]$points
        tmp[ is.na(tmp) ] <- 0.
        cumulResids <- cumulResids + tmp
        #browser()       
      }
      resids[[iPlot+1]]$cumulObs   <- cumulObs
      resids[[iPlot+2]]$cumulPreds <- cumulPreds
      resids[[iPlot+3]]$cumulResids<- cumulResids
  
      if ( .USEYEAR )
      { 
        xPos <- seq( .INITYEAR,.INITYEAR+nT-1, 5 )
        panLab( 0.5, 0.9, cex=1.0, paste( i+.INITYEAR-1,",",i,sep="") )
      }
      else      
        panLab( 0.5, 0.9, cex=1.0, paste( "Year",i ) )
      
      mfg <- par( "mfg" )
      
      axis( side=1, cex.axis=.CEXAXIS, pos=0 )
      
      if ( mfg[2]==1 )
        axis( side=2, cex.axis=.CEXAXIS, las=.YAXISLAS )
      else
        axis (side=2, labels=FALSE )
        
    }     # i=1,..,nT
    mtext( side=3, line=0, cex=.CEXLAB, outer=TRUE, 
              switch(g,
                     "Commercial",
                     "RV_4VWX",
                     "Halibut Survey" 
                     )
          ) 
    mtext( side=1, line=.OUTLINE1, cex=.CEXLAB, outer=TRUE, "Length class (cm)" )
    mtext( side=2, line=2, cex=.CEXLAB, outer=TRUE, "Proportion-Female-at-length" )

  }
  return( resids )
}