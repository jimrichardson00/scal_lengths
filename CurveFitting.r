# setwd("/home/jim/Dropbox/REM/tasks/scal_lengths")

# define age frequency vector
Age <- rexp(1000,rate=0.115)
Age <- Age[Age >= 1 & Age <= 30]

CurveFitting <- function(L_inf=205,k=0.1,t0=0.45,cv=0.1,Age=Age){

breaks <- 30

# define von bert function
Von_Bert <- function(t){
	return(L_inf*(1-exp(-k*(t-t0))))
}

# need to defined this selectivity function
selectivity <- function(t){
	1/(1 + 1*exp(-0.5*(t - 81)))
}
plot(seq(0,2*81,1),selectivity(seq(0,2*81,1)),type='l')

# set fitted x values and y values
Xfit <- seq(min(Age),max(Age),length=breaks)
Yfit <- Von_Bert(Xfit)

# set y limits
Ylim <- range(Xfit,max(Von_Bert(Xfit) + 4*cv*Von_Bert(Xfit)))
Xlim <- range(Xfit)

# age histogram
h <- hist(Age,breaks=c(0,Xfit),plot=FALSE)

# calculate theoretical length distribution
TheoreticalLengthDist <- vector()
for(i in seq(1,length(Xfit),1)){
	X <- Xfit[i]
	TheoreticalLengthDist <- c(TheoreticalLengthDist,rnorm(h$counts[i]*100 + 1,mean=Von_Bert(X),sd=cv*Von_Bert(X)))
}
dens_t <- density(TheoreticalLengthDist)

# calculate observed length
ObservedLengthDist <- vector()
for(i in seq(1,length(Xfit),1)){
	X <- Xfit[i]
	LenDis_i <- rnorm(h$counts[i]*10 + 1,mean=Von_Bert(X),sd=cv*Von_Bert(X))
	LenSel_i <- vector()
	for(l in LenDis_i){
		LenSel_i <- c(LenSel_i,rep(l,floor(1000*selectivity(l))))
	}
	ObservedLengthDist <- c(ObservedLengthDist,LenSel_i)
}
dens_o <- density(ObservedLengthDist)

# #####################################################################################

curve <- function(){

	set.seed(5431)

	plot(range(Xfit),range(Yfit),type='n',ylim=Ylim,xlim=Xlim,
		xaxt="n",yaxt="n",xlab="",ylab="")
	lines(Xfit,Yfit,type='l',ylim=Ylim,xlim=Xlim,
		xaxt="n",yaxt="n",xlab="",ylab="")

	Xpon <- seq(from=5,to=25,length=5)

	for(i in seq(from=1,to=length(Xpon),by=1)){

		X <- Xpon[i]

		mean <- Von_Bert(X)
		sd <- cv*Von_Bert(X)
		Von_Bert(X)
		X


		Yran <- c(mean - 4*sd, mean + 4*sd)
		Ypon <- seq(from=min(Yran),to=max(Yran),length=30)

		lines(rep(X,30),Ypon,xaxt="n",yaxt="n",xlab="",ylab="")

		Ydist <- dnorm(Ypon,mean=mean,sd=sd)

		lines(rep(X,30) + (X*10)*Ydist,Ypon,xaxt="n",yaxt="n",xlab="",ylab="")

	}

	#adds axis
	axis(1,at=pretty(Xlim),labels=pretty(Xlim),las=1)
	mtext("Age",side=1,line=2)

	axis(2,at=pretty(Ylim),labels=format(pretty(Ylim), scientific=FALSE),las=1)
	mtext("Length",side=2,line=3)

}


# #####################################################################
# First page

page1 <- function(){

layout(matrix(c(1,4,2,3), 2, 2, byrow = TRUE))

# ---------------------------
# top left

frame()

# --------------------------
# bottom left

plot(range(Xfit),range(Yfit),type='n',ylim=Ylim,xlim=Xlim,
	xaxt="n",yaxt="n",xlab="",ylab="")
lines(Xfit,Yfit,type='l',ylim=Ylim,xlim=Xlim,
	xaxt="n",yaxt="n",xlab="",ylab="")
#adds axis
axis(1,at=pretty(Xlim),labels=pretty(Xlim),las=1)
mtext("Age",side=1,line=2)
axis(2,at=pretty(Ylim),labels=format(pretty(Ylim), scientific=FALSE),las=1)
mtext("Length",side=2,line=3)

# -------------------------
# bottom right

frame()

}

# #####################################################################
# Second page

page2 <- function(){

layout(matrix(c(1,4,2,3), 2, 2, byrow = TRUE))

# ---------------------------
# top left

frame()

# --------------------------
# bottom left

curve()

# -------------------------
# bottom right

frame()

}

# #####################################################################
# Third page

page3 <- function(){

layout(matrix(c(1,4,2,3), 2, 2, byrow = TRUE))

# ---------------------------
# top left

hist(Age,breaks=30,xlim=Xlim,xlab='Age',main="")

# --------------------------
# bottom left

curve()

# -------------------------
# bottom right

frame()

}

# #####################################################################
# Fourth page

page4 <- function(){

layout(matrix(c(1,4,2,3), 2, 2, byrow = TRUE))

# ---------------------------
# top left

hist(Age,breaks=30,xlim=Xlim,xlab='Age',main='')

# --------------------------
# bottom left

curve()

# -------------------------
# bottom right

plot(dens_t$y,dens_t$x,type="l",ylim=Ylim,,xlim=range(0,dens_o$y,dens_t$y),lty="dashed",
	xlab='Density',ylab='Length')
legend("topright",c("Theoretical"),lty=c("dashed"),bty="n")

}

# #####################################################################
# Fifth page

page5 <- function(){

layout(matrix(c(1,4,2,3), 2, 2, byrow = TRUE))

# ---------------------------
# top left

hist(Age,breaks=30,xlim=Xlim)

# --------------------------
# bottom left

curve()
par( new=TRUE )
plot(range(selectivity(Yfit)),range(Yfit),type="n", xaxt="n",yaxt="n", xlab="", ylab="" )
lines(selectivity(Yfit),Yfit, xaxt="n",yaxt="n", xlab="", ylab="" )
axis(3,at=pretty(c(0,1)),labels=pretty(c(0,1)),las=1)

# -------------------------
# bottom right

plot(dens_o$y,dens_o$x,type="l",ylim=Ylim,,xlim=range(0,dens_o$y,dens_t$y),
	xlab='Density',ylab='Length')
lines(dens_t$y,dens_t$x,type="l",ylim=Ylim,lty="dashed")
legend("topright",c("Theoretical","Observed"),lty=c("dashed","solid"),bty="n")

}

# #####################################################################
# output to jpeg

jpeg("Page 1.jpeg")
page1()
dev.off()

jpeg("Page 2.jpeg")
page2()
dev.off()

jpeg("Page 3.jpeg")
page3()
dev.off()

jpeg("Page 4.jpeg")
page4()
dev.off()

jpeg("Page 5.jpeg")
page5()
dev.off()

}