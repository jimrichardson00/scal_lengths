#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#
#-- Plots for scal Atlantic Halibut models                                   --#
#-- Most of these functions are based on mseR plotting conventions
# Notes:
#  1. Need to plot the selectivity functions by fishery and year. 
#     Need a 3-panel figure, one plot for each fishery. Put
#     multiple years on same graph, but fisheries on different graphs.
#  2. Need to output Bayesian MCMC samples in scal.tpl, then make a 
#     function here to read them and produce summary tables. There is
#     some code already here to do much of that, but will need
#     customizing from crab to halibut.
#-----------------------------------------------------------------------#
rm( list=ls() )
source("mseRtools.r")
source("scal_globals.r")
#source("mcmcpairs.r")
require(MCMCpack)

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

#------------------------------------------------------------------------------#
#-- Usage
#   These plots all require a list object argument ("obj") defined
#      by the ADMB rep file. For gbcod, the list will be "scal.rep", which
#      can be input and formatted to a list via the R command line:
#      > repObj <- lisread( "scal.rep" )
#   The above list in "repFileName" can then be used as the arguments to the 
#   function.
#      An example for plotting Mt output is:
#      > plotAll_M( obj=repObj )
#   IMPORTANT: if the ADMB output is changed, then repObj needs to be re-created.
#      Otherwise, the output plots will use the old repObj stored in the current
#      R environment.
# catch-equation w discarding
# > f<-expression(s*F*p*B*(1-exp(-(s*F*(p+(1-p)*D)+M)))/(s*F*(p+(1-p)*D)+M))
# derivative wrt F
# > D(f,"F")
# (s * p * B * (1 - exp(-(s * F * (p + (1 - p) * D) + M))) + s * 
#     F * p * B * (exp(-(s * F * (p + (1 - p) * D) + M)) * (s * 
#     (p + (1 - p) * D))))/(s * F * (p + (1 - p) * D) + M) - s * 
#     F * p * B * (1 - exp(-(s * F * (p + (1 - p) * D) + M))) * 
#     (s * (p + (1 - p) * D))/(s * F * (p + (1 - p) * D) + M)^2
# > 

plotLenResiduals <- function( repObj )
{

plot( x=NULL, y=NULL, xlim=c(1,nT), ylim=c(1,nAges),
        xlab="Year",ylab="Age", main="residPropAge" )
  resCol <- vector("numeric")
  for( t in 1:nT )
  {
    resid <- 5.*(result$obsPropAge[t,]-result$predPropAge[t,])
    resCol[resid > 0] <- "blue"
    resCol[resid <= 0] <- "red"
    points( rep(t,nAges),c(1:nAges), cex=resid, col=resCol )  
  }


}
  
# plotBiomass
#   plot Total, Legal, and SSB time-series.
#   Could add a vector specifying which series
#   to plot in case user doesn't want them all.
plotBiomass <- function( obj, scaler=1.e3 )
{
  nT       <- obj$nT
  lenAge_f <- obj$lenAge_f
  lenAge_m <- obj$lenAge_m

  Bta_m <- obj$Bta_m
  Bta_f <- obj$Bta_f
  totalBt <- rowSums(Bta_m + Bta_f)

  # Use length-at-age to determine whether each element of
  # the biomass-at-age is legal and return legal biomass
  # only
  setLegal <- function( bvec,type="m",sizeLimit=15 )
  {
    legal <- rep(0.,length(bvec))
    lenAge <- switch(type,
                     m = lenAge_m,
                     f = lenAge_f
                     )
    legal[ lenAge >= sizeLimit ] <- 1.
    # return the vector of legal biomass-at-age
    return( bvec*legal )
  }
  # subset the biomass-at-age matrices for the years in which
  # the size limit was active
  idxSizeLimitYear <- obj$sizeLimYear-.INITYEAR+1
  tmpB_m  <- Bta_m[idxSizeLimitYear:nT,]
  
  # compute legal biomass-at-age by year
  tmpB_m  <- apply( X=tmpB_m, MARGIN=1, FUN=setLegal, 
                    type="m", sizeLimit=obj$sizeLimit[1] )
  
  # replace original biomass-at-age with legal amounts
  Bta_m[idxSizeLimitYear:nT,] <- t(tmpB_m)
  # repeat for females
  tmpB_f  <- Bta_f[idxSizeLimitYear:nT,]
  tmpB_f   <- apply( X=tmpB_f, MARGIN=1, FUN=setLegal, 
                    type="f", sizeLimit=obj$sizeLimit[1] )
  Bta_f[idxSizeLimitYear:nT,] <- t(tmpB_f)
  legalBt <- rowSums( Bta_m + Bta_f )

  # reset graphical parameters.
  graphics.off()
  xLim   <- c( .INITYEAR, .LASTYEAR )
  yLimit <- c( 0, max(legalBt) )/scaler
  plot( xLim, yLimit, type="n", axes=F, xlab="Year", ylab="" )

  # Add three biomass time-series
  lines( c( xLim[1]:xLim[2] ), legalBt/scaler, col="black", lty="solid", lwd=2 )
  lines( c( xLim[1]:xLim[2] ), totalBt/scaler, col="black", lty="dashed", lwd=2 )
  lines( c( xLim[1]:xLim[2] ), obj$SSBt/scaler, col="black", lty="dotted", lwd=2 )
  
  xPos <- seq( .INITYEAR,.LASTYEAR, 2 )
  xLabs <- paste( xPos )
  axis( side=1, cex.axis=.CEXAXIS2, at=xPos, labels=xLabs )
  axis( side=2, cex.axis=.CEXAXIS2, las=.YAXISLAS )
  panLegend( 0,1,legTxt="Total",lty="solid",lwd=2, bty="n" )
  panLegend( 0,.95,legTxt="Legal",lty="dashed",lwd=2, bty="n" )
  panLegend( 0,.90,legTxt="Spawning",lty="dotted",lwd=2, bty="n" )

  mtext( side=2, line=2.5, text="Biomass (tonnes)", cex=.CEXAXIS2)
  box(bty="l")
   
} # end plotBiomass

# plotSelectivity
#   Plot the selectivity-at-age curves by fishery
#   for a given type: m=male, f=female
plotSelectivity <- function( obj, type="m" )
{

  nT <- 1
  ageClasses <- obj$ages
  nSelBlocks <- obj$nT
  lType      <- c("solid", "dashed", "dotted" )
  lCol       <- c("black","green","blue","red")

  sel_1       <- switch(type,
                      m = obj$sel_gta_m1,
                      f = obj$sel_gta_f1)
  sel_2       <- switch(type,
                      m = obj$sel_gta_m2,
                      f = obj$sel_gta_f2)
  sel_3       <- switch(type,
                      m = obj$sel_gta_m3,
                      f = obj$sel_gta_f3)
  nSelLines <- 0
  for( j in 1:nT )
  {
    nSelLines <- nSelLines+1
    if( j==1 )
    {
      #browser()
      plot( ageClasses, sel_1[1,], type="l", lty=lType[1], 
            lwd=2, bty="n", axes=F, xlab="",
            ylab="", col=lCol[1], ylim=c(0,1) )
      axis(side=1)
      axis(side=2,las=1)
    }
    else
    {
      sel       <- sel_1[j,]
      lines( ageClasses, sel, lty=lType[1],col=lCol[j], lwd=2 )  
    }
    panLegend( 0.5,1.-0.05*nSelLines,
               legTxt=paste("FISHERY_",1970+j-1,sep=""), bty="n", 
               lty=lType[1], col=lCol[j],cex=0.6 )
  }
  for( j in 1:nT )
  {
    nSelLines <- nSelLines + 1
    sel       <- sel_2[j,]
    lines( ageClasses, sel, lty=lType[2], lwd=2,col=j )
    panLegend( 0.5,1.-0.05*nSelLines,
               legTxt=paste("RV_",1970+j-1,sep=""), 
               bty="n", lty=lType[2], col=j,cex=0.6 )
  }

  for( j in 1:nT )
  {
    nSelLines <- nSelLines + 1
    sel       <- sel_3[j,]
    lines( ageClasses, sel, lty=lType[3], col=lCol[j], lwd=2 )
    panLegend( 0.5,1.-0.05*nSelLines,
               legTxt=paste("HS_",1970+j-1,sep=""), 
               bty="n", lty=lType[3], col=lCol[j],cex=0.6 )
  }

  mtext( side=1, line=2, outer=FALSE, text="Age", cex=1)
  mtext( side=2, line=3, outer=FALSE, text="Selectivity",cex=1)
#browser()
} # end plotSelectivity

plotWtAge <- function( obj, label=NULL,
                   gfx=list( annotate=TRUE, xLim=NULL, yLim=NULL ) )
{

  ages   <- obj$ages
  Wta_m    <- obj$wtAge_m*1.e3 # wtAge is in tonnes
  Wta_f    <- obj$wtAge_f*1.e3 # wtAge is in tonnes
  maxWt  <- max(max(Wta_m,Wta_f))
  maxAge <- max(ages)
  xRange <- c( 0,maxAge )

  matYr <- matrix(NA, nrow=35, ncol=35, byrow=T )
  matWt <- matrix(NA, nrow=35, ncol=35, byrow=T )
  for( cohRow in 1:35 )
  {
    maxCol <- min( cohRow + maxAge -1, ncol(matWt) )
    a <- 0
    for( j in cohRow:maxCol )
    {
      a <- a + 1
      matYr[ cohRow, j ] <- (.INITYEAR + cohRow -1) + a - 1
      matWt[ cohRow, j ] <- Wta[ cohRow, a ]
    }
  }
  xLim <- c( .INITYEAR,.LASTYEAR )
  yLim <- c( 0,max(Wta) )
  
  # Weight at age.
  plot( xLim,yLim,type="n",axes=FALSE,xlab="",ylab="" )
  for( t in 1:nrow(matWt) )
  {
    lines( matYr[t,], matWt[t,], lty=1, lwd=1 )
    points( matYr[t,], matWt[t,], bg="white", col="black", cex=1, pch=21 )
  }
  for( w in seq(2,14,by=2) )
    abline( h=w, lty="dotted", lwd=0.5, col="gray" )
  
  axis( side=1, cex.axis=.CEXAXIS )
  axis( side=2, cex.axis=.CEXAXIS, las=.YAXISLAS )
  box()
    
  mtext( side=1, line=2.5, cex=.CEXLAB, outer=FALSE, "Cohort/Year" )
  mtext( side=2, line=2, cex=.CEXLAB, outer=FALSE, "Weight (kg)" )
}     # .plotWgtAtAge function

# plotHt
# Plots exploitation rates by M group.
# These need to be computed in scal.tpl
plotHt <- function( obj )
{
  graphics.off()
  par( mar=c(6,6,2,2) )
  Ht        <- obj$subHt
  splitAge <- obj$splitAge
  nT        <- obj$nT
  katch     <- obj$landCatchMatrix[,3]
  
  xLim  <- c( .INITYEAR, .LASTYEAR )
  xVals <- c( xLim[1]:xLim[2] )
  xPos  <- seq( .INITYEAR,.INITYEAR+nT-1, 2 )
  xLabs <- paste( xPos )
  delta <- 0.2
  
  # Catch plot
  layout(matrix(c(1,1,2,2,2,2),3,2,byrow=T))  
  yLimit <- c( 0, max(katch) )
  plot( xLim, yLimit, type="n", axes=F, xlab="", ylab="" )
  
  # Catch bars.
  rect( xVals-delta, 0, xVals+delta, katch, col="gray88" )
  axis( side=1, cex.axis=.CEXAXIS2, at=xPos, labels=xLabs )
  axis( side=2, cex.axis=.CEXAXIS2, las=.YAXISLAS )

  mtext( side=1, line=3, cex=.CEXLAB, outer=TRUE, "Year" )
  mtext( side=2, line=4, cex=.CEXLAB, outer=FALSE, "Catch" )
  
  
  # Ht plot
  yLimit <- c( 0, max(Ht) )
  plot( xLim, yLimit, type="n", axes=F, xlab="", ylab="" )

  lines( xVals, Ht[1,],  col="black", lty="solid", lwd=1 )
  lines( xVals, Ht[2,], col="black", lty="dashed", lwd=1 )

  points( xVals, Ht[1,], bg="black", cex=1.2, pch=21 )
  points( xVals, Ht[2,], bg="white", cex=1.2, pch=21 )
  
  axis( side=1, cex.axis=.CEXAXIS2, at=xPos, labels=xLabs )
  axis( side=2, cex.axis=.CEXAXIS2, las=.YAXISLAS )

  mtext( side=1, line=3, cex=.CEXLAB, outer=FALSE, "Year" )
  mtext( side=2, line=4, cex=.CEXLAB, outer=FALSE, "Exploitation rate" )
  
  box()
  Hnames1 <- paste( "Ht_1", splitAge[2]-1, sep="" )
  Hnames2 <- paste( "Ht_",splitAge[2],"+", sep="")
  
  panLegend( 0.05,0.95, legTxt=c(Hnames1,Hnames2),
             col="black", lty=c(.LTYGROUP1,.LTYGROUP2),
             lwd=2, bty="n",
             pch=c(19,21) )

}     # plotHt function

# plotLenFit
# Plots observed length-frequency (bars) and fitted model length-frequency
# for each year. A separate page is generated for each fishery/survey
plotLenFit <- function( obj, modelCheck=F, doType="m", resGear=1, png=F )
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
  pl      <- vector("list")
  pl[[1]] <- matrix(NA, nrow=nLens,ncol=nT)
  pl[[1]] <- switch(doType,
                    m = obj$pl_m1,
                    f = obj$pl_f1,
                    c = obj$pl_c1
                    )
  pl[[2]] <- matrix(NA, nrow=nLens,ncol=nT)
  pl[[2]] <- switch(doType,
                    m = obj$pl_m2,
                    f = obj$pl_f2,
                    c = obj$pl_c2
                    )
  pl[[3]] <- matrix(NA, nrow=nLens,ncol=nT)
  pl[[3]] <- switch(doType,
                    m = obj$pl_m3,
                    f = obj$pl_f3,
                    c = obj$pl_c3
                    )
  
  lenClasses <- vector("list")
  lenClasses[[1]] <- seq( minLen[1], maxLen[1], by=obj$binSize )
  lenClasses[[2]] <- seq( minLen[2], maxLen[2], by=obj$binSize )
  lenClasses[[3]] <- seq( minLen[3], maxLen[3], by=obj$binSize )

  # These are the lengths.
  xAxisLengths    <- obj$lfBins
  xLim <- range( xAxisLengths )
  
  # OBS: lenObsProp is a 3-dim array with dimensions( nGear,nLens,nT ).
  lenObsProp      <- array( 0,dim=c(nGear,nLens,nT) )
  lenObsProp[1,firstBin[1]:lastBin[1],] <- switch(doType,
                                                  m = obj$lenObsProp_m1,
                                                  f = obj$lenObsProp_f1,
                                                  c = obj$lenObsProp_c1
                                                  )
  lenObsProp[2,firstBin[2]:lastBin[2],] <- switch(doType,
                                                  m = obj$lenObsProp_m2,
                                                  f = obj$lenObsProp_f2,
                                                  c = obj$lenObsProp_c2
                                                  )
  lenObsProp[3,firstBin[3]:lastBin[3],] <- switch(doType,
                                                  m = obj$lenObsProp_m3,
                                                  f = obj$lenObsProp_f3,
                                                  c = obj$lenObsProp_c3
                                                  )
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
    yLim <- c( 0,1 )

  for ( g in 1:nGear )
  {   

    if ( g > 1 )
    {
      dev.new()
    }
    
    yLim    <- c( 0, 0.25 )
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
      lines( lenClasses[[g]], pl[[g]][,i], lty=1, lwd=2, col="red" )
      points( lenClasses[[g]], pl[[g]][,i], bg="white", 
              col="black", cex=0.5, lwd=1, pch=21 )
      
      if( g==resGear )
      {
        iPlot <- iPlot + 1
        resids[[iPlot]]$points    <- vector("numeric",length(lenClasses[[g]]))  
        resids[[iPlot]]$stdResids <- vector("numeric",length(lenClasses[[g]]))  
        resids[[iPlot]]$std       <- vector("numeric",length(lenClasses[[g]]))  
        
        resids[[iPlot]]$points    <- lenObsProp[g,,i] - pl[[g]][,i]
        resids[[iPlot]]$std       <- sd( resids[[iPlot]]$points, na.rm=TRUE )  
        resids[[iPlot]]$stdResids <- resids[[iPlot]]$points/resids[[iPlot]]$std
        resids[[iPlot]]$stdResids[is.na(resids[[iPlot]]$stdResids)] <- 0
        resids[[iPlot]]$year    <- paste( i+.INITYEAR-1,",",i,sep="") 

        tmp <- lenObsProp[g,,i]
        tmp[is.na(tmp)] <- 0.
        cumulObs <- cumulObs + tmp
        
        tmp <- pl[[g]][,i]
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
    mtext( side=2, line=2, cex=.CEXLAB, outer=TRUE, "Proportion-at-length" )

  }
  return( resids )
}     # .plotLenFit function

# plotIndexFit
# Plots observed total survey catch and fitted exploitable biomass - 
# biomass available via survey-specific selectivity. All data is 
# re-scaled to biomass units via division by "q"
# All surveys are shown on same page
plotIndexFit <- function( obj, ssb=FALSE )
{
  #png( filename="Index.png" )
  expBtg <- obj$expBtg
  Itg    <- obj$Itg
  Itg[ Itg < 0 ] <- NA
  nT     <- obj$nT
  nIndex <- nrow(expBtg)
  xLim <- c( .INITYEAR, .LASTYEAR )
  scaler <- c(1.e6,1.e3)
  
  myMar <- c( 3,5,0.5,0.5 )
  par( oma=c(1,3,1,1), mar=myMar, mfcol=c(nIndex,1) )
  legs <- c("(A)","(B)")
  yMax1  <- max( max(expBtg[1,],na.rm=TRUE), max(Itg[1,],na.rm=TRUE)  )
  yMax2  <- max( max(expBtg[2,],na.rm=TRUE), max(Itg[2,],na.rm=TRUE)  ) 
  yMax   <- c( yMax1/scaler[1], yMax2/scaler[2] )
  for( i in 1:nIndex )
  {
    yLimit <- c( 0, yMax[i] )
    yLab <- c("Number of halibut (millions)",
              "Halibut biomass (1000s tonnes)")
    plot( xLim, yLimit, type="n", axes=F, xlab="", 
          ylab=yLab[i] )

    # Model biomass
    lines( c( xLim[1]:xLim[2] ), expBtg[i,]/scaler[i], col="black", 
           lty="solid", lwd=2 )       
  
    # Biomass index data.
    points( c(xLim[1]:xLim[2] ), Itg[i,]/scaler[i], bg="white", cex=1.8, pch=21 )

    xPos  <- seq( .INITYEAR,.INITYEAR+nT-1, 2 )
    xLabs <- paste( xPos )
    axis( side=1, cex.axis=.CEXAXIS2, at=xPos, labels=xLabs )
    axis( side=2, cex.axis=.CEXAXIS2, las=.YAXISLAS )  
  }
  mtext( side=1, line=3, text="Year" )

  #dev.off()
}     # plotIndexFit function

# plotAll_R
# Plots estimated total age-1 recruitment (top) and recruitment deviations
# from random walk (bottom). The parameters controlling recruitment are 
# priorSD_R (prior standard deviation on recruitment residuals) and 
# gamma_R, which is the estimated autocorrelation in recruitment. 
plotAll_R <- function( obj )
{
  par( oma=c(3.5,4,1,1), mar=c(2,2,1,1), mfrow=c(2,1 ) )
  plotRecruitment( obj=obj )
  plotRecDevs( obj=obj )
}

plotRecruitment <- function( obj, ssb=FALSE, newPlot=F )
{

  # This will need to use an Rt variable.
  Rt     <- 2.*obj$Nta_m[,1]/1.e6
  nT     <- obj$nT
  if( !newPlot) 
  {
    par( oma=.OMA)
    cexY <- .CEXAXIS
  }
  else
    cexY <- 0.8
    
  xLim <- c( .INITYEAR, .LASTYEAR )
  yLimit <- c( 0, max(Rt) )
  plot( xLim, yLimit, type="n", axes=F, xlab="Year", ylab="" )

  # Recruits
  lines( c( xLim[1]:xLim[2] ), Rt, col="black", lty="solid", lwd=2 )
  
  if( exists(obj$avgR) )
    avgR <- obj$avgR/1.e6
  else
    avgR <- mean( obj$Rt )
  abline( h=avgR, lty="dashed", col="gray" )
  
  xPos <- seq( .INITYEAR,.LASTYEAR, 2 )
  xLabs <- paste( xPos )
  axis( side=1, cex.axis=.CEXAXIS2, at=xPos, labels=xLabs )
  axis( side=2, cex.axis=cexY, las=.YAXISLAS )
  mtext( side=2, line=3, text="Recruitment (millions)", cex=cexY)
  box(bty="l")
   
}     # plotRecruitment function

plotRecDevs <- function( obj, label=NULL, gfx=list( annotate=TRUE, bygears=FALSE,
                            doLegend=TRUE, xLim=NULL, yLim=NULL ) )
{
  # Recruitment deviations.

  recDevs <- obj$recDevs
  priorSD_R <- round( obj$priorSD_R, digits=2)
  gamma_R   <- round( obj$gamma_R, digits=2)
  nT      <- length( recDevs )
  
  xLim <- gfx$xLim
  yLim <- gfx$yLim

  # X-axis limits.
  if ( is.null(xLim) )
    xLim <- c( 1,nT )

  # Y-axis limits.
  if ( is.null(yLim) )
    yLim <- range( recDevs,na.rm=TRUE )
    
  plot( xLim, yLim, type="n", axes=FALSE, xlab="", ylab="" )

  abline( h=0, lty=3 )
      
  lines( c(1:nT), recDevs, col="black", lty=1, lwd=1 )
  points( c(1:nT), recDevs, bg="black", cex=1.2, pch=21 )

  
  axis( side=1, cex.axis=.CEXAXIS )
  axis( side=2, cex.axis=.CEXAXIS, las=.YAXISLAS )
  box()
  
  mtext( side=1, line=.OUTLINE1, cex=.CEXLAB, outer=TRUE, "Year" )
  mtext( side=2, line=.OUTLINE2, cex=.CEXLAB, outer=TRUE, "Recruitment Deviation" )
  
  parNames    <- vector(length=2)
  parNames[1] <- paste( "priorSD_R = ", priorSD_R, sep="" )
  parNames[2] <- paste( "gamma_R = ", gamma_R, sep="" )
  
  panLegend( 0.05,0.35, legTxt=parNames,
             col="black", lty=c(.LTYGROUP1,.LTYGROUP2),
             pch=c(19,21),
             lwd=2, bty="n" )

  return( invisible() )
}     # .plotMLErecDevs


#------------------------------------------------------------------------------#
# MCMC Plotting Functions                                                      #
#------------------------------------------------------------------------------#
getMCMC <- function( parName )
{ 
  parsHdr<-c("MSY","FMSY","initF","B0")
  if( any( parsHdr==parName ) )
  {
    fName <- paste("mc",parName,".dat",sep="")
    pars  <- as.matrix( read.table(file="mcPar.dat", header=TRUE,sep="") )
    mcPar <- pars[ , parName ]
  }
  else
  {
    fName <- paste("mc",parName,".dat",sep="")
    mcPar <- as.matrix( read.table(file=fName, header=F, sep="") )
  }
   
  return( mcPar )
}

# drawPoly
#   This is a generic function meant to draw a polygon on whatever
#   list object is given.
drawPoly <- function( obj, new=FALSE )
{
  # obj is a list. If length(obj)>1 polys will overlay
  # obj$x       - x-axis values
  # obj$yMatrix - y values where quantile will be used on columns.
  #             - MLE assumed in 1st row
  # obj$color   - vector of rgb colors
  if( new )
    par( oma=c(4,4,1,1), mar=c(3,2,1,1) )
  for( i in 1:length(obj) )
  {
    x       <- obj[[i]]$x
    rev.x   <- rev(x) 
    yMatrix <- obj[[i]]$yMatrix
    cols    <- obj[[i]]$color
    quants  <- apply( X=yMatrix, FUN=quantile, MARGIN=2, prob=c(0.025, 0.5, 0.975) )
    yMean   <- apply( X=yMatrix, FUN=mean, MARGIN=2 )
    if( i==1 )
    {
      xLim    <- range(x)
      yLim    <- c(0,max(quants[3,]))
      plot( xLim, yLim, type="n", xlab="",ylab="",axes=F )
    }
    polygon( c(x, rev.x), c(quants[1,], rev(quants[3,])), 
             col=rgb(cols[1],cols[2],cols[3],alpha=0.1), 
             lty=NULL )
      lines( x, quants[3,], lwd=1.5 )
      lines( x, quants[2,], lwd=2 )
      lines( x, quants[1,], lwd=1.5 )
  
  }
  axis( side=1 )
  axis( side=2, las=1 ) 
  
  box( bty="l" )
} 
# adapt this to handle the 3 biomass
# time-series plotted in plotBiomass
# function. Fix mcout header
# to match what getMCMC is looking
# for, e.g., "Bio_t"
makeBioList <- function()
{
  x <- getMCMC("Bio")
  idx <- seq(1,nrow(x)-1,by=2)
  x.pre <- x[idx,]
  idx <- seq(2,nrow(x),by=2)
  x.post <- x[idx,]
  
  obj <- vector("list",length=2)
  obj[[1]]$x <- 1991:2012
  obj[[2]]$x <- 1991:2012

  obj[[1]]$yMatrix <- x.pre
  obj[[2]]$yMatrix <- x.post

  obj[[1]]$color <- c(0,0,1)
  obj[[2]]$color <- c(1,0,0)
  
  return( obj )
}
plotPriorPost <- function( parName, writeCSV="new" )
{
  
  z.rep <- lisread("lbmcrab.rep")
  if( parName=="Par" )
  {
    parsHdr<-c("avgR","M1","M2","beta1","beta2","sigG")
    prNames <- c("pr_log_AvgR", "pr_log_M","pr_log_M","pr_beta1", "pr_beta2", "pr_sigG")
    xlabs   <- c( "Average male recruitment (millions)",
                  "Non-moulting mortality (/yr)",
                  "Moulting mortality (/yr)",
                  "Absolute moult increment (mm)",
                  "Relative moult increment",
                  "Growth standard deviation (mm)" 
                 )
                 
    lb      <- c( 0, 0.2, 0.2, 16, 0.6, 1 )
    ub      <- c( 10, 0.8, 0.8, 32, 1.2, 5)
    nPar   <- length(parsHdr)
    zNames <- vector(mode="character",length=nPar)
    zSumm  <- matrix(NA, nrow=nPar, ncol=4)
    zRow <- 0
    par( mfcol=c(3,2), oma=c(3,3,1,1), mar=c(3,2,2,2) )
    for( i in 1:length(parsHdr) )
    {       
      pName <- parsHdr[i]
      z.mc  <- getMCMC( pName )
      zVals <- z.mc 
    #if(i==5) browser()
      # create breaks based on prior distribution
      repName <- paste( "z.rep$", prNames[i], sep="" )
      prVals <- eval(parse(text=repName))
      if( i==1 )
      {
        mu   <- exp(prVals[1])/1.e6
        zVals <- zVals/1.e6
      }
      else
      {
        mu   <- prVals[1]
      }
      sd   <- prVals[2]
      
      brks <- "Sturges" #seq( from=lb[i], to=ub[i], length=30 )
      hst <- hist( x=zVals, axes=F, breaks=brks, main=" ")
      brks <- hst$breaks
      axis( side=1 )
      mtext( side=1, line=2.5, text=xlabs[i], cex=0.8 )

      # draw prior
      xVals  <- seq(min(brks),max(brks),length=100)
      if( i !=6 )
        pr_pdf <- dnorm( x=xVals, mean=mu, sd=sd )
      else
        pr_pdf <- dinvgamma( x=xVals, shape=mu, scale=sd )
      lines( xVals, max(hst$counts)*pr_pdf/max(pr_pdf) )
      
      # compute posterior summary and save to summary matrix
      zRow           <- zRow + 1
      zNames[zRow]   <- pName
      quants         <- quantile(x=zVals, probs=c(0.025,0.5,0.975))
      means          <- mean(zVals)
      zSumm[zRow,]   <- round( c( quants[1:2],means,quants[3] ),2 )
    }
  }
  if( parName != "Par" )
  {
    xVals <- seq(from=92.5,to=207.5,by=2.5)
    z.mc  <- getMCMC( parName )  
  }
  if( parName=="Rec" )
  { 
    par( oma=c(4,4,1,1), mar=c(2,2,1,1) )
    plotHST.Rec( hst=z.mc, rep=z.rep )
    
    nVals <- ncol(z.mc)
    tVals <- seq(1991,2012,by=1)
    zNames <- vector(mode="character",length=nVals)
    zSumm  <- matrix(NA, nrow=nVals, ncol=4)
    zRow   <- 0
    for( i in 1:nVals )
    {
      zNames[i] <- paste( "Rec_",tVals[i],sep="" )
    }
    quants <- apply(X=z.mc,MARGIN=2,FUN=quantile,probs=c(0.025,0.5,0.975))
    means  <- colMeans(x=z.mc)
    zSumm[1:nVals,1]  <- round( quants[1,],2 )/1.e6
    zSumm[1:nVals,2]  <- round( quants[2,],2 )/1.e6
    zSumm[1:nVals,3]  <- round( means,2 )/1.e6
    zSumm[1:nVals,4]  <- round( quants[3,],2 )/1.e6
    
  }
  if( parName=="rho" )
  { 
    par( oma=c(4,4,1,1), mar=c(2,2,1,1) )
    plotHST.rho( hst=z.mc, rep=z.rep )
    
    nVals <- ncol(z.mc)
    lVals <- z.rep$cutPoints[1:5]
    zNames <- vector(mode="character",length=nVals)
    zSumm  <- matrix(NA, nrow=nVals, ncol=4)
    zRow   <- 0
    for( i in 1:nVals )
    {
      zNames[i] <- paste( "rho_",lVals[i],sep="" )
    }
    quants <- apply(X=z.mc,MARGIN=2,FUN=quantile,probs=c(0.025,0.5,0.975))
    means  <- colMeans(x=z.mc)
    zSumm[1:nVals,1]  <- round( quants[1,],2 )
    zSumm[1:nVals,2]  <- round( quants[2,],2 )
    zSumm[1:nVals,3]  <- round( means,2 )
    zSumm[1:nVals,4]  <- round( quants[3,],2 )
    
  }
  if( parName=="N1l" )
  { 
    par( oma=c(4,4,1,1), mar=c(2,2,1,1) )
    plotHST.N1l( hst=z.mc, rep=z.rep )
    
    nVals <- ncol(z.mc)
    lVals <- z.rep$cutPoints[6:z.rep$max_initN_class]
    zNames <- vector(mode="character",length=nVals)
    zSumm  <- matrix(NA, nrow=nVals, ncol=4)
    zRow   <- 0
    for( i in 1:nVals )
    {
      zNames[i] <- paste( "init_N_",lVals[i],sep="" )
    }
    quants <- apply(X=z.mc,MARGIN=2,FUN=quantile,probs=c(0.025,0.5,0.975))
    means  <- colMeans(x=z.mc)
    zSumm[1:nVals,1]  <- round( quants[1,],2 )
    zSumm[1:nVals,2]  <- round( quants[2,],2 )
    zSumm[1:nVals,3]  <- round( means,2 )
    zSumm[1:nVals,4]  <- round( quants[3,],2 )
    
  }
  if( parName=="PMoult" )
  {   
    ind   <- which( z.mc > max(xVals), arr.ind=TRUE )
    if( nrow(ind) > 0 )
      z.mc[ ind[,1], ind[,2] ] <- max(xVals)-1
    nPar  <- z.rep$nPBlocks
    par( mfcol=c(nPar,1), oma=c(3,3,1,1), mar=c(1,2,1,1) )
    p95Col <- "black"
    p50Col <- "gray"
    zNames <- vector(mode="character",length=2*nPar)
    zSumm  <- matrix(NA, nrow=2*nPar, ncol=4)
    zRow <- -1
    for( i in 1:nPar )
    {
      zVals <- z.mc[,c(i,i+nPar)]
      z.hst <- apply(X=zVals,MARGIN=2,FUN=hist,plot=F, breaks=xVals) 
      plotHST.pMoult( hst=z.hst[[1]], rep=z.rep, addHST=FALSE, bCol=p95Col ) 
      plotHST.pMoult( hst=z.hst[[2]], rep=z.rep, addHST=TRUE, bCol=p50Col )
      x      <- paste("z.rep$tBlock_P_",i,sep="")
      tBlock <- eval(parse(text=x))    
      parName1 <- paste("P95_",1990+tBlock,sep="")
      parName2 <- paste("P50_",1990+tBlock,sep="")
      panLegend( 0.,1.,legTxt=parName1, 
                 bty="n", lty="solid", col=p95Col,cex=0.6 )
      panLegend( 0.,.85,legTxt=parName2, 
                 bty="n", lty="solid", col=p50Col,cex=0.6 )
      axis( side=1 )

      # compute posterior summary and save
      zRow           <- zRow + 2
      print(zRow)
      zNames[zRow]   <- parName1
      zNames[zRow+1] <- parName2
      quants         <- apply(X=zVals,MARGIN=2,FUN=quantile,probs=c(0.025,0.5,0.975))
      means          <- colMeans(x=zVals)
      zSumm[zRow,]   <- round( c( quants[1:2,1],means[1],quants[3,1] ),2 )
      zSumm[zRow+1,] <- round( c( quants[1:2,2],means[2],quants[3,2] ), 2 )
    }
    mtext( side=1, line=2, text="Length (mm)", outer=T )
  }
  if( grepl(x=parName,pattern="Sel") )
  {   
    ind   <- which( z.mc > max(xVals), arr.ind=TRUE )
    if( nrow(ind) > 0 )
      z.mc[ ind[,1], ind[,2] ] <- max(xVals)-1
    gdx  <- as.integer( strsplit(x=parName, split="l")[[1]][2] )
    gear <- switch( gdx, "Fishery","DFO_PRE","DFO_POST" )
    
    nPar  <- z.rep$nSelBlocks[gdx]
    par( mfcol=c(2,nPar), oma=c(3,3,1,1), mar=c(2,2,1,1) )
    p50Col <- "black"
    p95Col <- "gray"
    zNames <- vector(mode="character",length=2*nPar)
    zSumm  <- matrix(NA, nrow=2*nPar, ncol=4)
    zRow <- -1
    for( i in 1:nPar )
    {  
      zVals <- z.mc[,c(i,i+nPar)]
      #browser()
      z.hst <- apply(X=zVals,MARGIN=2,FUN=hist,plot=F,breaks=xVals) 
      plotHST.Sel( hst=z.hst[[1]], rep=z.rep, addHST=FALSE, bCol=p50Col )
      plotHST.Sel( hst=z.hst[[2]], rep=z.rep, addHST=TRUE,  bCol=p95Col )
      x      <- paste("z.rep$tBlock_S_",gdx,i,sep="")
      tBlock <- eval(parse(text=x))
      parName1 <- paste(gear,"_S50_",1990+tBlock,sep="")
      parName2 <- paste(gear,"_S95_",1990+tBlock,sep="")
      panLegend( 0.5,.92,legTxt=parName1, 
                 bty="n", lty="solid", col=p50Col,cex=0.5 )
      panLegend( 0.5,.83,legTxt=parName2, 
                 bty="n", lty="solid", col=p95Col,cex=0.5 )
      axis( side=1 )
      
      # calc selectivity function for each row of pars
      sel   <- matrix(NA,nrow=nrow(zVals), ncol=100 )
      lVals <- seq( 100, 200, length=ncol(sel) )
      for( i in 1:nrow(zVals) )
      {         
        sel[i,] <- 1./( 1. + exp( -log(19.)*(lVals-zVals[i,1])/(zVals[i,2]-zVals[i,1]) ) )
      }
      obj              <- vector("list",length=1)
      obj[[1]]$x       <- lVals
      obj[[1]]$yMatrix <- sel
      obj[[1]]$color   <- c(0,0,1)
      drawPoly(obj, new=FALSE)
      mtext( side=2, line=3, text="Selectivity", outer=F, cex=0.8)

      # compute posterior summary and save
      zRow           <- zRow + 2
      print(zRow)
      zNames[zRow]   <- parName1
      zNames[zRow+1] <- parName2
      quants         <- apply(X=zVals,MARGIN=2,FUN=quantile,probs=c(0.025,0.5,0.975))
      means          <- colMeans(x=zVals)
      zSumm[zRow,]   <- round( c( quants[1:2,1],means[1],quants[3,1] ), 2 )
      zSumm[zRow+1,] <- round( c( quants[1:2,2],means[2],quants[3,2] ), 2 )
    }
    # add xlab for whole page
    mtext( side=1, line=2, text="Length (mm)", outer=T )
    
  }
  
  # assemble data frame and write param summary to csv table
  z.df <- data.frame( nams=zNames, sums=zSumm ) 
  names(z.df) <- c("Parameter", "2.5%", "Median", "Mean", "97.5%")
  print(z.df)
  
  if( !is.null( writeCSV ) )
  {
    if( writeCSV=="new" )
    {
      write.table(x=z.df,file="summPars.csv", sep=",", col.names=TRUE, row.names=FALSE, append=F)
    }
    if( writeCSV=="app" )
    {
      write.table(x=z.df,file="summPars.csv", sep=",", col.names=F, row.names=FALSE, append=T)
    }  
  }


}
plotHST.Prod <- function( hst, rep, writeCSV=NULL )
{   
  
  hst <- getMCMC(parName="Prod")
  par( oma=c(4,4,1,1), mar=c(2,3,1,1) )

  minT     <- 1992
  maxT     <- 2012
  tClasses <- seq(minT,maxT,by=1)
  
  # rect spacing
  delta  <- 0.25
  quants <- apply(X=hst,MARGIN=2,FUN=quantile,probs=c(0.025,0.5,0.975),breaks="Sturges")
  means  <- colMeans(x=hst)

  xLim <- c(1988,2014)
  yLim <- c(0, 1.2*max(hst))
  plot( xLim, yLim, type="n", axes=FALSE, xlab="", ylab="" )
  for( i in 1:length(tClasses) )
  { 
    jitX <- runif( n=nrow(hst),min=tClasses[i]-delta, max=tClasses[i]+delta )
    points( jitX, hst[,i], pch=19, cex=0.05, col="gray", bg="white" )
    
    points( tClasses[i], means[i], pch=19, col="black", cex=0.8 )
    
    rect( tClasses[i]-delta, quants[1,i], 
          tClasses[i]+delta, quants[3,i], 
          bg="transparent" )
    segments( x0=tClasses[i]-delta, x1=tClasses[i]+delta, 
              y0=quants[2,i], y1=quants[2,i],
              col="black",
              lwd=1.5 )
  }
  
  axis( side=1 )
  axis( side=2,las=1 )
  mtext( side=1, line=3, outer=F, text="Year", cex=1)
  mtext( side=2, line=3, outer=F, text="Legal crab production (tonnes)",cex=1)
  
  zNames <- vector(mode="character",length=21)
  zSumm  <- matrix(NA,nrow=21,ncol=4)
  for( i in 1:21 )
  {
    zNames[i] <- paste("Prod_",1991+i,sep="")
  }
  quants <- apply(X=hst,MARGIN=2,FUN=quantile,probs=c(0.025,0.5,0.975))
  means  <- colMeans(x=hst)
  zSumm[1:21,1]  <- round( quants[1,],2 )
  zSumm[1:21,2]  <- round( quants[2,],2 )
  zSumm[1:21,3]  <- round( means,2 )
  zSumm[1:21,4]  <- round( quants[3,],2 )
  # assemble data frame and write param summary to csv table
  z.df <- data.frame( nams=zNames, sums=zSumm ) 
  names(z.df) <- c("Parameter", "2.5%", "Median", "Mean", "97.5%")
  print(z.df)
  
  if( !is.null(writeCSV) )
  {
    write.table(x=z.df,file="summProd.csv", sep=",", col.names=TRUE, row.names=FALSE, append=F)
  }
  
  
}


plotHST.Rec <- function( hst, rep )
{   
  # vector of length classes
  minT     <- 1991
  maxT     <- 2012
  tClasses <- seq(minT,maxT,by=1)
  
  # rect spacing
  delta  <- 0.25
  quants <- apply(X=hst,MARGIN=2,FUN=quantile,probs=c(0.025,0.5,0.975),breaks="Sturges")
  means  <- colMeans(x=hst)

  # prior bands
  mu <- exp( rep$pr_log_AvgR[1] )
  sd <- rep$pr_log_AvgR[2]*1e6
  xLim <- c(1988,2014)
  yLim <- c(1e5, 1e8)
  plot( xLim, yLim, type="n", log="y", axes=FALSE, xlab="", ylab="" )
  for( i in 1:length(tClasses) )
  { 
    jitX <- runif( n=nrow(hst),min=tClasses[i]-delta, max=tClasses[i]+delta )
    points( jitX, hst[,i], pch=19, cex=0.05, col="gray", bg="white" )
    
    points( tClasses[i], means[i], pch=19, col="black", cex=0.8 )
    
    rect( tClasses[i]-delta, quants[1,i], 
          tClasses[i]+delta, quants[3,i], 
          bg="transparent" )
    segments( x0=tClasses[i]-delta, x1=tClasses[i]+delta, 
              y0=quants[2,i], y1=quants[2,i],
              col="black",
              lwd=1.5 )
  }
  
  axis( side=1 )
  axis( side=2,las=1 )
  mtext( side=1, line=1, outer=TRUE, text="Year", cex=1)
  mtext( side=2, line=2, outer=TRUE, text="Recruitment",cex=1)
  
}

plotHST.rho <- function( hst, rep )
{   
  # vector of length classes
  minLen     <- 1
  maxLen     <- 5
  lenClasses <- rep$cutPoints[minLen:maxLen]
  
  # rect spacing
  delta  <- 1
  quants <- apply(X=hst,MARGIN=2,FUN=quantile,probs=c(0.025,0.5,0.975), breaks="Sturges")
  means  <- colMeans(x=hst)
  
  xLim <- c(97.5,127.5)
  yLim <- c(0,1)  
  plot( xLim, yLim, type="n", axes=FALSE, xlab="", ylab="" )
  for( i in 1:length(lenClasses) )
  { 
    jitX <- runif( n=nrow(hst),min=lenClasses[i]-delta, max=lenClasses[i]+delta )
    points( jitX, hst[,i], pch=19, cex=0.05, col="gray", bg="white" )
    
    points( lenClasses[i], means[i], pch=19, col="black", cex=0.8 )
    
    rect( lenClasses[i]-delta, quants[1,i], 
          lenClasses[i]+delta, quants[3,i], 
          bg="transparent" )
    segments( x0=lenClasses[i]-delta, x1=lenClasses[i]+delta, 
              y0=quants[2,i], y1=quants[2,i],
              col="black",
              lwd=1.5 )
  }
  axis( side=1 )
  axis( side=2,las=1 )
  mtext( side=1, line=1, outer=TRUE, text="Carapace width (mm)", cex=1)
  mtext( side=2, line=2, outer=TRUE, text="Proportion recruiting",cex=1)
  
}
plotHST.N1l <- function( hst, rep )
{   
  # vector of length classes
  minLen     <- 6
  maxLen     <- rep$max_initN_class
  lenClasses <- rep$cutPoints[minLen:maxLen]
  
  # rect spacing
  delta  <- 1
  quants <- apply(X=hst,MARGIN=2,FUN=quantile,probs=c(0.025,0.5,0.975), breaks="Sturges")
  means  <- colMeans(x=hst)
  
  xLim <- range(lenClasses)
  yLim <- c(100,100*max(quants))  
  plot( xLim, yLim, type="n", log="y", axes=FALSE, xlab="", ylab="" )
  for( i in 1:length(lenClasses) )
  { 
    jitX <- runif( n=nrow(hst),min=lenClasses[i]-delta, max=lenClasses[i]+delta )
    points( jitX, hst[,i], pch=19, cex=0.05, col="gray", bg="white" )
    
    points( lenClasses[i], means[i], pch=19, col="black", cex=0.8 )
    
    rect( lenClasses[i]-delta, quants[1,i], 
          lenClasses[i]+delta, quants[3,i], 
          bg="transparent" )
    segments( x0=lenClasses[i]-delta, x1=lenClasses[i]+delta, 
              y0=quants[2,i], y1=quants[2,i],
              col="black",
              lwd=1.5 )
  }
  axis( side=1 )
  axis( side=2,las=1 )
  mtext( side=1, line=1, outer=TRUE, text="Carapace width (mm)", cex=1)
  mtext( side=2, line=2, outer=TRUE, text="Male crab abundance in 1991",cex=1)
  
}
plotHST.Sel <- function( hst, rep, addHST=FALSE, bCol="red" )
{

  # Plot hist frequencies as segments
  freq      <- hst$counts/max(hst$counts)
  cutPoints <- hst$breaks[1:length(freq)]
  #browser()
  # create new plot if necessary
  if( addHST==FALSE )
  {
    #browser()
    plot( cutPoints, freq, type="n", 
          xlim=range(cutPoints), 
          ylim=c(0,1),
          axes=F,
          xlab="")
  }
  # draw hist bars individually
  for( l in 1:length(cutPoints) )
  {
    # add segments for histogram
    segCol <- bCol
    delta  <- ifelse( addHST, 0.5, -0.5 )
        x0 <- cutPoints[l] + delta
        x1 <- cutPoints[l] + delta
        y0 <- rep(0,length(x0))
        y1 <- freq[l]

    segments( x0=x0,x1=x1,y0=y0,y1=y1,
             col=segCol,
              lty="solid",
              lwd=3
            )
  }
  # add prior distribution
  xVals  <- seq( from=min(cutPoints), to=max(cutPoints), length=1000 )
  if( addHST==FALSE )
  {
    sd <- rep$pr_S50[2]
    mu <- rep$pr_S50[1]
  }
  else
  {
    sd <- getSD( df=c(1,exp(2.*rep$pr_delatS[1])), sig=c(rep$pr_S50[2],rep$pr_deltaS[2]) )
    mu <- rep$pr_S50[1] + exp(rep$pr_deltaS[1])
  }
  pr_pdf <- dnorm( x=xVals, mean=mu, sd=sd )
  # re-scale to hist y scale
  pr_pdf <- pr_pdf/max(pr_pdf) 
  # add line
  lines( xVals, pr_pdf, lty="solid", col=bCol, lwd=1.25 )
}

# obj    - a histogram object
# addHST - add a hist to existing plot?
plotHST.pMoult <- function( hst, rep, addHST=FALSE, bCol="red" )
{
  # Plot hist frequencies as segments
  freq      <- hst$counts/max(hst$counts)
  cutPoints <- hst$breaks[1:length(freq)]
  # create new plot if necessary
  if( addHST==FALSE )
  {
    #browser()
    plot( cutPoints, freq, type="n", 
          xlim=range(cutPoints), 
          ylim=c(0,1),
          axes=F,
          xlab="")
  }
  # draw hist bars individually
  for( l in 1:length(cutPoints) )
  {
    # add segments for histogram
    segCol <- bCol
    delta  <- ifelse( addHST, 0.5, -0.5 )
        x0 <- cutPoints[l] + delta
        x1 <- cutPoints[l] + delta
        y0 <- rep(0,length(x0))
        y1 <- freq[l]

    segments( x0=x0,x1=x1,y0=y0,y1=y1,
             col=segCol,
              lty="solid",
              lwd=3
            )
  }
  # add prior distribution
  xVals  <- seq( from=min(cutPoints), to=max(cutPoints), length=1000 )
  if( addHST==FALSE )
  {
    sd <- rep$pr_P95[2]
    mu <- rep$pr_P95[1]
  }
  else
  {
    sd <- getSD( df=c(1,exp(2.*rep$pr_delatP[1])), sig=c(rep$pr_P95[2],rep$pr_deltaP[2]) )
    mu <- rep$pr_P95[1] + exp(rep$pr_deltaP[1])
  }
  pr_pdf <- dnorm( x=xVals, mean=mu, sd=sd )
  # re-scale to hist y scale
  pr_pdf <- pr_pdf/max(pr_pdf) 
  # add line
  lines( xVals, pr_pdf, lty="solid", col=bCol, lwd=1.25 )
}

summaryBio <- function( writeCSV=NULL )
{
  # get biomass posterior summary and append to posterior.csv file.
  obj    <- makeBioList()
  z.mc   <- obj[[1]]$yMatrix
  zNames <- vector(mode="character",length=22)
  zSumm  <- matrix(NA,nrow=22,ncol=4)
  for( i in 1:22 )
  {
    zNames[i] <- paste("Legal_B_",obj[[1]]$x[i],sep="")
  }
  quants <- apply(X=z.mc,MARGIN=2,FUN=quantile,probs=c(0.025,0.5,0.975))
  means  <- colMeans(x=z.mc)
  zSumm[1:22,1]  <- round( quants[1,],2 )
  zSumm[1:22,2]  <- round( quants[2,],2 )
  zSumm[1:22,3]  <- round( means,2 )
  zSumm[1:22,4]  <- round( quants[3,],2 )
  # assemble data frame and write param summary to csv table
  z.df <- data.frame( nams=zNames, sums=zSumm ) 
  names(z.df) <- c("Parameter", "2.5%", "Median", "Mean", "97.5%")
  print(z.df)
  
  if( !is.null(writeCSV) )
  {
    write.table(x=z.df,file="summBio.csv", sep=",", col.names=TRUE, row.names=FALSE, append=F)
  }

  drawPoly( obj )
  mtext( side=1, line=3, text="Year" )
  mtext( side=2, line=3, text="Legal biomass (tonnes)" )
  
}

# Delta approximation to sd of functions
# df - vector of first derivs
# sig2 <- vector of variances
getSD <- function( df, sig )
{                       
  d <- sum( df*df*sig*sig )
  return( sqrt(d) )
}

summaryAll <- function( repObj )
{              
  png( file="Par.png" )
    plotPriorPost( parName="Par", writeCSV="new" )
  dev.off()

  png( file="PMoult.png" )
  plotPriorPost( parName="PMoult", writeCSV="app" )
  dev.off()

  png( file="Sel1.png" )
  plotPriorPost( parName="Sel1", writeCSV="app" )
  dev.off()

  png( file="Sel2.png" )
  plotPriorPost( parName="Sel2", writeCSV="app" )
  dev.off()

  png( file="Sel3.png" )
  plotPriorPost( parName="Sel3", writeCSV="app" )
  dev.off()

  png( file="Rec.png" )
  plotPriorPost( parName="Rec", writeCSV="app" )
  dev.off()

  png( file="rho.png" )
  plotPriorPost( parName="rho", writeCSV="app" )
  dev.off()

  png( file="N1l.png" )
  plotPriorPost( parName="N1l", writeCSV="app" )
  dev.off()

  png( file="LegalB.png" )
    summaryBio( writeCSV="new" )
  dev.off()
  
  png( file="Prod.png" )
    plotHST.Prod( hst=NULL, rep=repObj, writeCSV="new" )
  dev.off()

  png( file="mcPairs-Pars.png" )
    plotPairs( parName="Pars" )
  dev.off()

  #png( file="mcPairs-PMoult.png" )
  #  plotPairs( parName="PMoult" )
  #dev.off()

  png( file="mcPairs-Sel.png" )
    plotPairs( parName="Sel" )
  dev.off()
  

}

plotPairs <- function( parName="Pars" )
{ 

  if( parName=="Pars" )
  {
    mcObj <- matrix(NA,nrow=1000,ncol=6)
    parsHdr<-c("avgR","M1","M2","beta1","beta2","sigG")
    for( i in 1:length(parsHdr) )
    {       
      pName <- parsHdr[i]
      mcObj[,i]  <- getMCMC( pName )
      if( i==1 ) mcObj[,i] <- mcObj[,i]/1.e6
    }
    colnames(mcObj) <- parsHdr
  }
  if( parName=="PMoult")
  {
    mcObj <- matrix(NA,nrow=1000,ncol=6)
    pm <- getMCMC(parName="PMoult")
    iCol <- 0
    for( j in 1:3 )
    { 
      iCol <- iCol + 1
      mcObj[,iCol]<- pm[,j]
      iCol <- iCol + 1
      mcObj[,iCol] <- pm[,j+ncol(pm)/2]
    }
    colnames( mcObj ) <- c("P95_1","P50_1","P95_2","P50_2","P95_3","P50_3")
  }
  if( parName=="Sel")
  {  
    repObj <- lisread("lbmcrab.rep")
    nCol   <- sum( repObj$nSelBlocks ) - 1
    mcObj  <- matrix(NA,nrow=1000,ncol=nCol*2)
    mcObj[,1:2] <- getMCMC(parName="Sel1")
    cNames      <- c("S50_F","S95_F")
    sp          <- getMCMC(parName="Sel2")
    iCol <- 2  
    for( j in 1:repObj$nSelBlocks[2] )
    { 
      iCol <- iCol + 1
      mcObj[,iCol]<- sp[,j]
      iCol <- iCol + 1
      mcObj[,iCol] <- sp[,j+ncol(sp)/2]
      cNames <- c( cNames, paste("S50_DFO_",j,sep=""),paste("S95_DFO_",j,sep="") )
    }
    colnames( mcObj ) <- cNames
  }
  .mcmcPairs( mcObj )
}




