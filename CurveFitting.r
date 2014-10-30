# setwd("/home/jim/Dropbox/REM/tasks/scal_lengths")

library(animation)

M <- 0.15
n_age <- data.frame(age=rep(NA,30),n_age=rep(NA,30))
n_age$age <- seq(1,30,1)
n_age[1,"n_age"] <- 80
for(i in seq(2,30,1)){
	n_age[i,"n_age"] <- ceiling(n_age[i-1,"n_age"]*exp(-M))
}

CurveFitting <- function(L_INF=205,K=0.1,T0=0.45,CV=0.1,N_AGE=n_age){

breaks <- 30

# age
Age <- vector()
for(a in n_age$age){
	Age <- c(Age,rep(n_age[a,"age"],ceiling(n_age[a,"n_age"])))
}

# define von bert function
Von_Bert <- function(t){
	return(L_INF*(1-exp(-K*(t-T0))))
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
Ylim <- range(Xfit,max(Von_Bert(Xfit) + 4*CV*Von_Bert(Xfit)))
Xlim <- range(Xfit)

# age histogram
h <- hist(Age,breaks=c(0,Xfit),plot=FALSE)

length_dens <- function(N_AGE=n_age,both=TRUE){

# calculate theoretical length distribution
TheoreticalLengthDist <- vector()
for(i in seq(1,length(Xfit),1)){
	X <- Xfit[i]
	TheoreticalLengthDist <- c(TheoreticalLengthDist,rnorm(N_AGE[i,"n_age"]*10 + 1,mean=Von_Bert(X),sd=CV*Von_Bert(X)))
}
dens_t <- density(TheoreticalLengthDist)

# calculate observed length
ObservedLengthDist <- vector()
for(i in seq(1,length(Xfit),1)){
	X <- Xfit[i]
	LenDis_i <- rnorm(N_AGE[i,"n_age"]*10 + 1,mean=Von_Bert(X),sd=CV*Von_Bert(X))
	LenSel_i <- vector()
	for(l in LenDis_i){
		LenSel_i <- c(LenSel_i,rep(l,floor(1000*selectivity(l))))
	}
	ObservedLengthDist <- c(ObservedLengthDist,LenSel_i)
}
dens_o <- density(ObservedLengthDist)

if(both == TRUE){

plot(dens_o$y,dens_o$x,type="l",ylim=Ylim,,xlim=range(0,dens_o$y,dens_t$y,0.015),
	xlab='Density',ylab='Length')
lines(dens_t$y,dens_t$x,type="l",ylim=Ylim,lty="dashed")
legend("topright",c("Theoretical","Observed"),lty=c("dashed","solid"),bty="n")

} else {

plot(dens_t$y,dens_t$x,type="l",ylim=Ylim,,xlim=range(0,dens_o$y,dens_t$y,0.015),
	xlab='Density',ylab='Length',lty="dashed")
legend("topright",c("Theoretical"),lty=c("dashed"),bty="n")

}

}

# #####################################################################################

curve <- function(){

	# age
	Age <- vector()
	for(a in n_age$age){
		Age <- c(Age,rep(n_age[a,"age"],ceiling(n_age[a,"n_age"])))
	}

	set.seed(5431)

	plot(range(Xfit),range(Yfit),type='n',ylim=Ylim,xlim=Xlim,
		xaxt="n",yaxt="n",xlab="",ylab="")
	lines(Xfit,Yfit,type='l',ylim=Ylim,xlim=Xlim,
		xaxt="n",yaxt="n",xlab="",ylab="")

	Xpon <- seq(from=5,to=25,length=5)

	for(i in seq(from=1,to=length(Xpon),by=1)){

		X <- Xpon[i]

		mean <- Von_Bert(X)
		sd <- CV*Von_Bert(X)
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

hist(Age,breaks=30,xlim=Xlim,xlab='Age',main="",ylim=range(0,300))

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

hist(Age,breaks=30,xlim=Xlim,xlab='Age',main='',ylim=range(0,300))

# --------------------------
# bottom left

curve()

# -------------------------
# bottom right

length_dens(N_AGE=n_age,both=FALSE)

}

# #####################################################################
# Fifth page

page5 <- function(){

layout(matrix(c(1,4,2,3), 2, 2, byrow = TRUE))

# ---------------------------
# top left

hist(Age,breaks=30,xlim=Xlim,main="",ylim=range(0,300))

# --------------------------
# bottom left

curve()
par( new=TRUE )
plot(range(selectivity(Yfit)),range(Yfit),type="n", xaxt="n",yaxt="n", xlab="", ylab="" )
lines(selectivity(Yfit),Yfit, xaxt="n",yaxt="n", xlab="", ylab="" )
axis(3,at=pretty(c(0,1)),labels=pretty(c(0,1)),las=1)

# -------------------------
# bottom right

length_dens(N_AGE=n_age,both=TRUE)

}

# #####################################################################
# Sixth page

page6 <- function(){

n_age_years <- list()
n_age_years[[1]] <- n_age
input_sd <- 0.3
age1 <- vector()
age1 <- c(age1,80)
for(y in seq(2,44,1)){
	M <- 0.15
	n_age_y <- data.frame(age=rep(NA,30),n_age=rep(NA,30))
	n_age_y$age <- seq(1,30,1)
	n_age_y[1,"n_age"] <- 80*exp(rnorm(1,mean=0,sd=input_sd))
	age1 <- c(age1,n_age_y[1,"n_age"])
	for(i in seq(2,30,1)){
		n_age_y[i,"n_age"] <- n_age_y[i-1,"n_age"]*exp(-M)
	}
	n_age_years[[y]] <- n_age_y
}

for(y in seq(1,44,1)){

n_age_y <- n_age_years[[y]]

Age <- vector()
for(a in n_age_y$age){
	Age <- c(Age,rep(n_age_y[a,"age"],ceiling(n_age_y[a,"n_age"])))
}

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))

# ---------------------------
# top left

hist(Age,breaks=30,xlim=Xlim,ylim=range(0,300),main="")

# ---------------------------
# top right

plot(seq(1,44,1),c(age1[seq(1,y,1)],rep(NA,44-y)),type='l',xlab="Years",ylab="Age-1 recruitment",ylim=range(age1))

# --------------------------
# bottom left

curve()
par( new=TRUE )
plot(range(selectivity(Yfit)),range(Yfit),type="n", xaxt="n",yaxt="n", xlab="", ylab="" )
lines(selectivity(Yfit),Yfit, xaxt="n",yaxt="n", xlab="", ylab="" )
axis(3,at=pretty(c(0,1)),labels=pretty(c(0,1)),las=1)

# -------------------------
# bottom right

length_dens(N_AGE=n_age_y,both=TRUE)

}

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

saveGIF(page6())

}