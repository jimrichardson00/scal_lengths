#Curve fitting example
curve <- function(a=.25,b=40,sig=.2,diff=0)
{
set.seed(5431)

#detrministic age,length values
plusGroupAge <- 30
Linf         <- 210
vonK         <- 0.12

detAge    <- seq(1,30,1)
detLength <- Linf*(1.-exp(-vonK*detAge))

# Choose 
Xran <- runif(n=10,min=5,max=max(Xdet))
e    <- rnorm(n=length(Xran),mean=0,sd=sig)
Yran<-a*Xran*exp(e-sig*sig/2)/(b+Xran)
plot(Xran,Yran,xlim=c(0,max(Xdet)),ylim=c(0,1.5*max(Ydet)),
		pch=19,xlab="Age",ylab="Length-at-age (cm)",
			cex=2,bty="l")

x1<-c(Xran[3:4],Xran[6:8])
y1<-c(Yran[3:4],Yran[6:8])
Ydet<-a*Xdet/(b+diff+Xdet)
	lines(Xdet,Ydet,lwd=2)
	lnL<-0
	#for each data point
	for(i in 1:5){
		#get point on curve
		ymean<-a*x1[i]/(b+diff+x1[i])
		#min and max of x-axis (vertical) for distribution as +- X sd's from the mean
		ymin<-ymean-2*sig*ymean;ymax<-ymean+4*sig*ymean
		#draw vertical axis
		lines(c(x1[i],x1[i]),c(ymin,ymax),lwd=2,col="blue")
		
		#x2 is just the x points
		x2<-rep(x1[i],30)
		#generate sequence of 30 increasing y-values
		y<-seq(ymin,ymax,length=30)
		#calc likelihood of each y
		py<-dnorm(log(y),mean=log(ymean),sd=sig)
		#add scaled likelihood to x values
		x2<-x2+2*py
		#plot x vs. y
		lines(x2,y,lwd=2,col="blue")
		
		lnL<-lnL+log(sig)/2-1/(2*sig*sig)*(log(ymean)-log(y1[i]))^2
	}
	print(c(b+diff,lnL))
}
