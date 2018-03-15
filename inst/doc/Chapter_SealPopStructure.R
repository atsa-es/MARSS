###################################################
### code chunk number 2: Cs00_required_libraries
###################################################
library(MARSS)


###################################################
### code chunk number 9: Cs01_set.up.data
###################################################
years = harborSeal[,1] #first col is years
#leave off Hood Canal data for now
sealData = t(harborSeal[,c(2:7,9:13)])


###################################################
### code chunk number 10: Cs02_fig1
###################################################
par(mfrow=c(4,3),mar=c(2,2,2,2))
for(i in 2:dim(harborSeal)[2]) {
    plot(years, harborSeal[,i], xlab="", ylab="", main=colnames(harborSeal)[i])
}


###################################################
### code chunk number 11: Cs03_set.up.Z.models
###################################################
#H1 stock
Z1=factor(c("wa.or","wa.or",rep("ps",4),"ca","ca","wa.or","wa.or","bc")) 
#H2 coastal+PS
Z2=factor(c(rep("coast",2),rep("ps",4),rep("coast",4),"ps")) 
#H3 N and S
Z3=factor(c(rep("N",6),"S","S","N","S","N")) 
#H4 North Coast, Inland Strait, Puget Sound, South Coast
Z4=factor(c("nc","nc","is","is","ps","ps","sc","sc","nc","sc","is"))
#H5 panmictic
Z5=factor(rep("pan",11)) 
#H6 Site
Z6=factor(1:11) #site
Z.models=list(Z1,Z2,Z3,Z4,Z5,Z6)
names(Z.models)=
     c("stock","coast+PS","N-S","NC+Strait+PS+SC","panmictic","site")


###################################################
### code chunk number 12: Cs04_Q.models
###################################################
Q.models=c("diagonal and equal", "diagonal and unequal")


###################################################
### code chunk number 13: Cs04a_other.models
###################################################
U.model="unequal"
R.model="diagonal and equal"
A.model="scaling"
B.model="identity"
x0.model="unequal"
V0.model="zero"
model.constant=list(
    U=U.model, R=R.model, A=A.model, 
    x0=x0.model, V0=V0.model, tinitx=0)


###################################################
### code chunk number 14: Cs05_run.the.models
###################################################
out.tab=NULL
fits=list()
for(i in 1:length(Z.models)){
  for(Q.model in Q.models){
     fit.model = c(list(Z=Z.models[[i]], Q=Q.model), model.constant)
     fit = MARSS(sealData, model=fit.model,
            silent=TRUE, control=list(maxit=1000))
    out=data.frame(H=names(Z.models)[i], Q=Q.model, U=U.model,
            logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
            m=length(unique(Z.models[[i]])),
            num.iter=fit$numIter, converged=!fit$convergence)
    out.tab=rbind(out.tab,out)
    fits=c(fits,list(fit))
    if(i==5) next #one m for panmictic so only run 1 Q
  }
}


###################################################
### code chunk number 15: Cs06_sort.results
###################################################
min.AICc=order(out.tab$AICc)
out.tab.1=out.tab[min.AICc,]


###################################################
### code chunk number 16: Cs07_add.delta.aicc
###################################################
out.tab.1=cbind(out.tab.1,
           delta.AICc=out.tab.1$AICc-out.tab.1$AICc[1])


###################################################
### code chunk number 17: Cs08_add.delta.aicc
###################################################
out.tab.1=cbind(out.tab.1, 
           rel.like=exp(-1*out.tab.1$delta.AICc/2))


###################################################
### code chunk number 18: Cs09_aic.weight
###################################################
out.tab.1=cbind(out.tab.1,
          AIC.weight = out.tab.1$rel.like/sum(out.tab.1$rel.like))


###################################################
### code chunk number 19: Cs10_print.table
###################################################
out.tab.1$delta.AICc = round(out.tab.1$delta.AICc, digits=2)
out.tab.1$AIC.weight = round(out.tab.1$AIC.weight, digits=3)
print(out.tab.1[,c("H","Q","delta.AICc","AIC.weight")], row.names=FALSE)


###################################################
### code chunk number 20: Cs11_fignorthsouth
###################################################
best.fit=fits[min.AICc][[1]]
matplot(years, t(best.fit$states-best.fit$states[,1]), 
        xlab="abundance index", ylab="",
        type="l",lwd=2,col="black")
legend("topleft",c("North Coastal","Inland Straits","Puget Sound","South Coastal"),lwd=2,lty=c(1:4),bty="n")


###################################################
### code chunk number 22: Cs12_new.Q.model
###################################################
for(i in 1:length(Z.models)){
    if(i==5) next #don't rerun panmictic
    for(Q.model in c("equalvarcov","unconstrained")){
      fit.model = c(list(Z=Z.models[[i]], Q=Q.model), model.constant)
      fit = MARSS(sealData, model=fit.model,
            silent=TRUE, control=list(maxit=1000))
       out=data.frame(H=names(Z.models)[i], Q=Q.model, U=U.model,
            logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
            m=length(unique(Z.models[[i]])),
            num.iter=fit$numIter, converged=!fit$convergence)
       out.tab=rbind(out.tab,out)
       fits=c(fits,list(fit))
    }
}


###################################################
### code chunk number 23: Cs13_out.tab.2
###################################################
min.AICc=order(out.tab$AICc)
out.tab.2=out.tab[min.AICc,]
fits=fits[min.AICc]
out.tab.2=cbind(out.tab.2,delta.AICc=out.tab.2$AICc-out.tab.2$AICc[1])
out.tab.2=cbind(out.tab.2,rel.like=exp(-1*out.tab.2$delta.AICc/2))
out.tab.2=cbind(out.tab.2,AIC.weight=out.tab.2$rel.like/sum(out.tab.2$rel.like))


###################################################
### code chunk number 24: Cs14_out.tab.2
###################################################
out.tab.2$AIC.weight = round(out.tab.2$AIC.weight, digits=3)
out.tab.2$delta.AICc = round(out.tab.2$delta.AICc, digits=2)
print(out.tab.2[1:10,c("H","Q","delta.AICc","AIC.weight")], row.names=FALSE)


###################################################
### code chunk number 25: Cs15_equalvarcov.weight
###################################################
c(
sum(out.tab.2$AIC.weight[out.tab.2$Q=="equalvarcov"]),
sum(out.tab.2$AIC.weight[out.tab.2$Q=="unconstrained"]),
sum(out.tab.2$AIC.weight[out.tab.2$Q=="diagonal and equal"])
)


###################################################
### code chunk number 26: Cs16_Q.mat
###################################################
Q.unc=coef(fits[[3]],type="matrix")$Q


###################################################
### code chunk number 27: Cs17_Q.diag
###################################################
diag(Q.unc)


###################################################
### code chunk number 28: Cs18_Q.corr
###################################################
h=diag(1/sqrt(diag(Q.unc)))
Q.corr=h%*%Q.unc%*%h
rownames(Q.corr)=unique(Z4)
colnames(Q.corr)=unique(Z4)

Q.corr


###################################################
### code chunk number 29: Cs19_add.hood.canal
###################################################
sealData.hc = rbind(sealData,harborSeal[,8])
rownames(sealData.hc)[12]="Hood.Canal"


###################################################
### code chunk number 30: Cs20_hood.z.models
###################################################
ZH1=factor(c("nc","nc","is","is","ps",
                     "ps","sc","sc","nc","sc","is","ps")) 
ZH2=factor(c("nc","nc","is","is","ps",
                     "ps","sc","sc","nc","sc","is","hc")) 
Z.models.hc=list(ZH1, ZH2)
names(Z.models.hc)=c("hood.in.ps","hood.separate")


###################################################
### code chunk number 31: Cs21_hood.uqr.models
###################################################
Q3=matrix(list("offdiag"),5,5)
diag(Q3)="q"
Q3[,5]=0; Q3[5,]=0; Q3[5,5]="q.hc"
Q.models=list("equalvarcov","unconstrained",Q3)
names(Q.models)=c("equalvarcov","unconstrained","hood.independent")


###################################################
### code chunk number 32: Cs22_hood-q3
###################################################
Q.models$hood.independent


###################################################
### code chunk number 33: Cs23_out.tab.hc
###################################################
out.tab.hc=NULL
fits.hc=list()
for(i in 1:length(Z.models.hc)){
  for(j in 1:length(Q.models)){
     if(i==1 & j==3) next #Q3 is only for Hood Separate model
     Q.model=Q.models[[j]]
     fit.model = c(list(Z=Z.models.hc[[i]], Q=Q.model), model.constant)
     fit = MARSS(sealData.hc, model=fit.model,
            silent=TRUE, control=list(maxit=1000))
    out=data.frame(H=names(Z.models.hc)[i], Q=names(Q.models)[j], U=U.model,
            logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
            m=length(unique(Z.models.hc[[i]])),
            num.iter=fit$numIter, converged=!fit$convergence)
    out.tab.hc=rbind(out.tab.hc, out)
    fits.hc=c(fits.hc,list(fit))
  }
}


###################################################
### code chunk number 34: Cs24_sort.aicc.hc
###################################################
min.AICc=order(out.tab.hc$AICc)
out.tab.hc=out.tab.hc[min.AICc,]
out.tab.hc=cbind(out.tab.hc, delta.AICc=out.tab.hc$AICc-out.tab.hc$AICc[1])
out.tab.hc=cbind(out.tab.hc,rel.like=exp(-1*out.tab.hc$delta.AICc/2))
out.tab.hc=cbind(out.tab.hc,AIC.weight=out.tab.hc$rel.like/sum(out.tab.hc$rel.like))


###################################################
### code chunk number 35: Cs25_out.tab.2
###################################################
out.tab.hc$AIC.weight = round(out.tab.hc$AIC.weight, digits=3)
out.tab.hc$delta.AICc = round(out.tab.hc$delta.AICc, digits=2)
print(out.tab.hc[,c("H","Q","delta.AICc","AIC.weight")], row.names=FALSE)


