CSEGriskfigure=function(data, te=100, absolutethresh=FALSE, threshold=0.1, datalogged=FALSE, silent=FALSE, return.model=FALSE, CI.method="hessian", CI.sim=1000) {
if(!(CI.method %in% c("hessian","parametric","innovations","none"))) 
  stop("Stopped in CSEGriskfigure because allowed CI methods are none, hessian, parametric, and innovations.\n", call.=FALSE)

if(!silent){
if(datalogged) cat("Analysis assumes that data and threshold are already logged")
else cat("Data and threshold are assumed to be unlogged")
cat("\n")
if(absolutethresh) {
	cat("Using an absolute threshold for the extinction threshold"); cat("\n")
	cat("Threshold is "); cat(threshold);
	cat("\n")
	if(datalogged) thresh=exp(threshold)
}
else  {
	cat("Using an percentage threshold for the extinction threshold"); cat("\n")
	cat("Threshold is "); cat((1-threshold)*100); cat(" percent decline")
	cat("\n")
}
}
par(mfrow=c(3,2))

#set up pi
pi=get("pi",pos="package:base")

#a=read.table(data, skip=dataskip)
a=as.matrix(data)
if(dim(a)[2]!=2) stop("Stopped in CSEGriskfigure() because this function requires a data file with 2 columns; time in first column.", call.=FALSE)
if(datalogged) a[,2]=exp(a[,2]) #plot unlogged
nyr=length(a[,1])

#estimate parameters
dat=a[,2]
dat=log(dat) #because I exp-ed it above if logged
kem = MARSS(dat, model=list(U=matrix("U")),silent=TRUE) 
if(kem$par$Q!=0) kem = MARSSparamCIs(kem)	

#PANEL #1 plot the data
plot(a[,1],a[,2],type="p", bty="L", xaxp=c(a[1,1],a[nyr,1],nyr-1),
xlab="", ylab="Pop. Estimate", ylim=c(.9*min(a[,2],na.rm=TRUE),1.1*max(a[,2],na.rm=TRUE)), xlim=c(a[1,1]-1,a[nyr,1]+1))
lines(a[,1],exp(kem$states[1,]),col="red")
title(
  paste("u est = ",format(kem$par$U,digits=2)," (95% CIs ",format(kem$par.lowCI[["U"]],digits=2),", ", format(kem$par.upCI[["U"]],digits=2),")", "\n","Q est = ",format(kem$par$Q,digits=2)),cex.main=.9)

#PANEL #2 plot the cdf of extinction
tyrs = 1:te 
kal.cdf  = matrix(0,nrow=te)
N0 = exp(kem$states[1,length(dat)])
if(!absolutethresh) thresh = threshold*N0 else thresh=threshold
xd = log(N0) - log(thresh)
kal.u = kem$par$U 
kal.Q = kem$par$Q
if(kal.u<=0) p.ever = 1 else p.ever = exp(-2*kal.u*xd/kal.Q) 
if(N0 <= thresh){ kal.cdf = rep(1,te) 
  }else{
  if(is.finite(exp(2*xd*abs(kal.u)/kal.Q))){
	 sec.part = exp(2*xd*abs(kal.u)/kal.Q) * pnorm((-xd - abs(kal.u)* tyrs)/sqrt(kal.Q*tyrs))
    }else sec.part=0      
  kal.cdf = p.ever*(pnorm(( -xd + abs(kal.u)*tyrs)/ sqrt(kal.Q*tyrs)) + sec.part)
}

plot(tyrs,kal.cdf,xlab="time steps into future",ylab="probability to hit threshold",ylim=c(0,1),bty="l",type="l")
title(paste("Prob. to hit ", format(thresh,digits=1)))

if(kem$par$Q!=0){

#Add CIs to cdf of quasi-extinction plot if requested
nsim=CI.sim
plotCI=TRUE
if(CI.method=="hessian"){
  plotCI=FALSE
  kem=MARSShessian(kem)
  #The R and Q are chol transformed so = sqrt(R) and sqrt(Q) respectively
  boot.params=try(rmvnorm(n=nsim, mean=kem$parMean, sigma=kem$parSigma, method="chol"), silent=TRUE)
    if(!inherits(boot.params, "try-error")) {
    plotCI=TRUE
    boot.params=boot.params[boot.params[,which(substr(colnames(boot.params),1,1)=="Q")]>0,]
    boot.params[,which(substr(colnames(boot.params),1,1)=="Q")]=boot.params[,which(substr(colnames(boot.params),1,1)=="Q")]^2
    boot.params[,which(substr(colnames(boot.params),1,1)=="R")]=boot.params[,which(substr(colnames(boot.params),1,1)=="R")]^2
    }
  }
if(CI.method=="parametric")
  boot.params=t( MARSSboot(kem, nboot=1000, output="parameters", sim="parametric")$boot.params ) 
if(CI.method=="innovations")
  boot.params=t( MARSSboot(kem, nboot=1000, output="parameters", sim="nonparametric")$boot.params )
if(CI.method!="none" & plotCI){  
  nsim = dim(boot.params)[1]
  kal.cdf.cis  = matrix(0,nrow=nsim,ncol=te)
  for(j in 1:nsim){
	 u=boot.params[j,which(substr(colnames(boot.params),1,1)=="U")]
	 Q=boot.params[j,which(substr(colnames(boot.params),1,1)=="Q")]
	 if(u<=0){ p.ever = 1
	 }else p.ever = exp(-2*u*xd/Q)     
   if(is.finite(exp(2*xd*abs(u)/Q))){
	 sec.part = exp(2*xd*abs(u)/Q) * pnorm((-xd - abs(u)* tyrs) / sqrt(Q*tyrs))
    }else sec.part=0      
   kal.cdf.cis[j,] = p.ever*pnorm(( -xd + abs(u)*tyrs)/ sqrt(Q*tyrs)) + sec.part
  }
    CIs.ex=apply(kal.cdf.cis,2,quantile,c(0.0275,.125,.875,.9725),na.rm=TRUE)
    lines(tyrs,CIs.ex[1,],col="red")
    lines(tyrs,CIs.ex[4,],col="red")
    lines(tyrs,CIs.ex[2,],col="green",lty=2)
    lines(tyrs,CIs.ex[3,],col="green",lty=2)
    legend("topleft",c("95% CI", "75% CI", "mean"), col=c("red","green","black"), lty=c(1,2,1), bty="n")
  }

#Panel 3 Now make plot of pdf of the time to extinction
kal.pdf  = matrix(nrow=te*5)
tyrs=1:(te*5)
if(N0 <= thresh) kal.pdf[1] = 1
else {  
for (i in 1:(te*5)){      
  kal.pdf[i] = 
	xd*(1/sqrt(2*pi*kal.Q*tyrs[i]^3))*
	exp(( -(xd-abs(kal.u)*tyrs[i])^2)/(2*kal.Q*tyrs[i]))
} # end i loop
}
i=tyrs
aa = max(kal.pdf)
i.left = 1:which(kal.pdf==aa)
i.min = max(which(kal.pdf[i.left] < aa*.001))
i.right = which(kal.pdf==aa):max(tyrs)
if(min(kal.pdf[i.right])< (aa*.001))
i.max = which(kal.pdf==aa)-1+min(which(kal.pdf[i.right] < aa*.001))
else i.max = max(tyrs)
plot(tyrs[i.min:i.max],kal.pdf[i.min:i.max],xlab="time steps into future",ylab="probability to hit threshold",bty="l",type="l")
title(main="PDF of time to threshold \n given it IS reached")

#Panel 4 Now make plot of probability of hitting different thresholds  in te yrs
thshes = seq(0.5*thresh,2*thresh,(2*thresh-0.5*thresh)/100)
kal.thresh  = matrix(0,nrow=length(thshes))
for (i in 1:length(thshes)){
  xd = log(N0) - log(thshes[i]) 
  if(kal.u<=0) p.ever = 1
  else p.ever = exp(-2*kal.u*xd/kal.Q)      
  if(N0 <= thshes[i]) kal.thresh[i] = 1  
  else kal.thresh[i] = 
	p.ever*pnorm(( -xd + abs(kal.u)*te)/ sqrt(kal.Q*te)) 
	   + exp(2*xd*abs(kal.u)/kal.Q) * pnorm((-xd - abs(kal.u)* te)/sqrt(kal.Q*te))
} # end i loop
plot(thshes,kal.thresh,xlab="Number of ind. at Ne",ylab="probability to hit threshold",bty="l",ylim=c(0,1),type="l")
abline(v=thresh, col="red")
text(thresh*1.1,max(kal.thresh)*.9,"90% threshold")
title(paste("Prob. of hitting threshold", " in ", te, " time steps") )

N=matrix(0,ncol=te,nrow=10)
for(i in 1:10)
  N[i,]=N0*cumprod(exp(rnorm(te,kal.u,sqrt(kal.Q))))
matplot(t(N),type="l",bty="L",main="Sample projections",xlab="time steps into the future",ylab="N")

#Finally make TMU plot
CSEGtmufigure(N=nyr, u=kal.u, s2p=kal.Q, make.legend=FALSE)

# KW
par(mfrow=c(1, 1))
if(return.model) return(kem)
}else{
  cat("Extinction PDF plots not possible because Q estimate is 0.\n")
  if(return.model) return(kem)
}
}




