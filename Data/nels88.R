cath=read.csv('nels.csv')
cath$faminc8 <- cath$faminc8-1
cath <- within(cath,{
             prop=fitted.values(glm(catholic~faminc8*math8+I(faminc8*faminc8),family=binomial))
             eotm=ifelse(catholic==1,1/prop,1/(1-prop))
})
cath <- within(cath,{
             teotm=eotm
             p=ifelse(catholic==1,eotm/sum(subset(cath,catholic==1,eotm)),eotm/sum(subset(cath,catholic==0,eotm)))
             hajek=ifelse(catholic==1,1/(nrow(cath[cath$catholic==1,])*p),1/(nrow(cath[cath$catholic==0,])*p))
})
cath$teotm[which(cath$teotm>20)]=20

## Make appendix figure 1 
par(mfrow=c(2,2))
hist(cath[which(cath$catholic==1),]$prop,breaks=seq(0,0.2,by=0.004),xlim=c(0,0.2),cex.lab=1.5,cex.axis=1.5,cex.main=1.5,main = "Histogram of propensity scores for catholic school students",xlab="Estimated Propensity Score")
hist(cath[which(cath$catholic==0),]$prop,breaks=seq(0,0.2,by=0.004),xlim=c(0,0.2),cex.lab=1.5,cex.axis=1.5,cex.main=1.5,main = "Histogram of propensity scores for public school students",xlab="Estimated Propensity Score")
hist(cath[which(cath$catholic==1),]$eotm,breaks = seq(0,50,by=1),xlim=c(0,50),cex.lab=1.5,cex.axis=1.5,cex.main=1.5,main = "Histogram of IPTW weights for catholic school students",xlab="Estimated IPTW Weights")
hist(cath[which(cath$catholic==0),]$eotm,breaks = seq(1,1.25,by=0.005),xlim=c(1,1.25),cex.lab=1.5,cex.axis=1.5,cex.main=1.5,main = "Histogram of IPTW weights for public school students",xlab="Estimated IPTW Weights")
par(mfrow=c(1,1))

summary(lm(math12~catholic+faminc8*math8,data=cath))
summary(lm(math12~catholic+faminc8*math8,weights = eotm,data=cath))
summary(lm(math12~catholic+faminc8*math8,weights = teotm,data=cath))

set.seed(12345)
library(boot)
## H-H estimator
hh=function(d,i){
  mod=lm(math12~catholic+faminc8*math8,data=d[i,],weights=d[i,]$hajek)
  return(coef(mod)[2])
}

hhest=boot(cath,hh,1000,strata=cath$catholic,weights=cath$p,parallel="snow",ncpus=8)
c(mean(hhest$t),sd(hhest$t)) ## Mean is 1.67, SD is 0.24

## Bootstrap IPTW
iptw=function(d,i){
  mod=lm(math12~catholic+faminc8*math8,data=d[i,],weights=d[i,]$eotm)
  return(coef(mod)[2])
}

biptw=boot(cath,iptw,1000,strata=cath$catholic,parallel="snow",ncpus=8)
c(mean(biptw$t),sd(biptw$t)) ## Mean is 1.61, SD is 0.25

## Bootstrap IPTW with trimmed weights
iptwt=function(d,i){
  mod=lm(math12~catholic+faminc8*math8,data=d[i,],weights=d[i,]$teotm)
  return(coef(mod)[2])
}

biptwt=boot(cath,iptwt,1000,strata=cath$catholic,parallel="snow",ncpus=8)
c(mean(biptwt$t),sd(biptwt$t)) ## Mean is 1.67, SD is 0.25

plot(density(hhest$t))
lines(density(biptw$t),col='red')
lines(density(biptwt$t),col='blue')
legend('topright',legend=c('H-H EST','BIPTW EST','BIPTW Trimmed'),lty=1,col=c('black','red','blue'))

## I create a function for bootstrapping IPTW estimator

bootiptw=function(iter,data){
  out <- numeric(iter)
  tdata <- data[data$catholic==1,]
  cdata <- data[data$catholic==0,]
  for(i in 1:iter){
    data1=tdata[sample(nrow(tdata),replace=TRUE),]
    data0=cdata[sample(nrow(cdata),replace=TRUE),]
    newdata=as.data.frame(rbind(data1,data0))
    iptw=lm(math12~catholic+math8*faminc8,weights=newdata$eotm,data=newdata)
    out[i] <- coef(iptw)[2]
  }
  return(c(mean(out),sd(out)))
}

bootiptw <- compiler::cmpfun(bootiptw)

bootiptw(1000,cath) ## 1.62 with se = 0.25

bootiptwt=function(iter,data){
  out <- numeric(iter)
  tdata <- data[data$catholic==1,]
  cdata <- data[data$catholic==0,]
  for(i in 1:iter){
    data1=tdata[sample(nrow(tdata),replace=TRUE),]
    data0=cdata[sample(nrow(cdata),replace=TRUE),]
    newdata=as.data.frame(rbind(data1,data0))
    iptw=lm(math12~catholic+math8*faminc8,weights=newdata$teotm,data=newdata)
    out[i] <- coef(iptw)[2]
  }
  return(c(mean(out),sd(out)))
}

bootiptwt <- compiler::cmpfun(bootiptwt)

bootiptwt(1000,cath) ## 1.68 with se = 0.25


## I also create a function for Hansen-Hurwitz estimator 
hhreg=function(iter,data){
  out <- numeric(iter)
  d1 <- data[data$catholic==1,]
  d0 <- data[data$catholic==0,]
  t1=rmultinom(iter,nrow(d1),prob=d1$p)
  t0=rmultinom(iter,nrow(d0),prob=d0$p)
  for(i in 1:iter){
    dt=d1[rep(1:nrow(d1),times=t1[,i]),]
    dc=d0[rep(1:nrow(d0),times=t0[,i]),]
    d=as.data.frame(rbind(dt,dc))
    mod=lm(math12~catholic+math8*faminc8,data=d,weights=d$hajek)
    out[i]=coef(mod)[2]
  }
  return(c(mean(out),sd(out)))
}

hhreg <- compiler::cmpfun(hhreg)

hhreg(1000,cath) ## 1.68 with se = 0.24


