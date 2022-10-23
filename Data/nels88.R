cath=read.csv('nels.csv')
cath$faminc8 <- cath$faminc8-1
cath <- within(cath,{
             prop=fitted.values(glm(catholic~faminc8*math8,family=binomial))
             eotm=ifelse(catholic==1,1/prop,1/(1-prop))
})
cath <- within(cath,{
             teotm=ifelse(eotm<=20,eotm,20)
             p=ifelse(catholic==1,eotm/sum(subset(cath,catholic==1,eotm)),eotm/sum(subset(cath,catholic==0,eotm)))
             hajek=ifelse(catholic==1,1/(nrow(cath[cath$catholic==1,])*p),1/(nrow(cath[cath$catholic==0,])*p))
})
cath <- within(cath,{
seotm <- ifelse(catholic==1,(sum(catholic)/length(catholic))/prop,(1-sum(catholic)/length(catholic))/(1-prop))
})

## Make appendix figure 1 
par(mfrow=c(2,2))
hist(cath[which(cath$catholic==1),]$prop,breaks=seq(0,0.2,by=0.004),xlim=c(0,0.2),cex.lab=1.5,cex.axis=1.5,cex.main=1.5,main = "Histogram of propensity scores for catholic school students",xlab="Estimated Propensity Score")
hist(cath[which(cath$catholic==0),]$prop,breaks=seq(0,0.2,by=0.004),xlim=c(0,0.2),cex.lab=1.5,cex.axis=1.5,cex.main=1.5,main = "Histogram of propensity scores for public school students",xlab="Estimated Propensity Score")
hist(cath[which(cath$catholic==1),]$eotm,breaks = seq(0,95,by=1),xlim=c(0,95),cex.lab=1.5,cex.axis=1.5,cex.main=1.5,main = "Histogram of IPTW weights for catholic school students",xlab="Estimated IPTW Weights")
hist(cath[which(cath$catholic==0),]$eotm,breaks = seq(1,1.25,by=0.005),xlim=c(1,1.25),cex.lab=1.5,cex.axis=1.5,cex.main=1.5,main = "Histogram of IPTW weights for public school students",xlab="Estimated IPTW Weights")
par(mfrow=c(1,1))

summary(lm(math12~catholic+faminc8*math8,data=cath)) ## 1.677
summary(glm(catholic~faminc8*math8,family=binomial,data=cath))
summary(lm(math12~catholic+faminc8*math8,weights = eotm,data=cath)) ## 1.60
summary(lm(math12~catholic+faminc8*math8,weights = teotm,data=cath)) ## 1.65

set.seed(12345)
library(boot)
## H-H estimator
hh=function(d,i){
  mod=lm(math12~catholic+faminc8*math8,data=d[i,],weights=d[i,]$hajek)
  return(coef(mod)[2])
}

hhest=boot(cath,hh,1000,strata=cath$catholic,weights=cath$p,parallel="multicore",ncpus=8)
c(mean(hhest$t),sd(hhest$t)) ## Mean is 1.67, SD is 0.246
quantile(hhest$t,probs = c(0.025,0.975))

## Bootstrap IPTW
iptw=function(d,i){
  dat=d[i,]
  dat=within(dat,{
    p=fitted.values(glm(catholic~math8*faminc8,family=binomial))
    eotm=ifelse(catholic==1,1/p,1/(1-p))
  })
  mod=lm(math12~catholic+faminc8*math8,data=dat,weights=dat$eotm)
  return(coef(mod)[2])
}

biptw=boot(cath,iptw,1000,strata=cath$catholic,parallel="multicore",ncpus=8)
c(mean(biptw$t),sd(biptw$t)) ## Mean is 1.59, SD is 0.248
quantile(biptw$t,probs = c(0.025,0.975))

## Bootstrap IPTW with trimmed weights
iptwt=function(d,i){
  dat=d[i,]
  dat=within(dat,{
    p=fitted.values(glm(catholic~math8*faminc8,family=binomial))
    eotm=ifelse(catholic==1,1/p,1/(1-p))
    teotm=ifelse(eotm<=20,eotm,20)
  })
  mod=lm(math12~catholic+faminc8*math8,data=dat,weights=dat$teotm)
  return(coef(mod)[2])
}

biptwt=boot(cath,iptwt,1000,strata=cath$catholic,parallel="multicore",ncpus=8)
c(mean(biptwt$t),sd(biptwt$t)) ## Mean is 1.65, SD is 0.247
quantile(biptwt$t,probs = c(0.025,0.975))

## Bootstrap IPTW with stabilized weights
iptws=function(d,i){
  dat=d[i,]
  dat=within(dat,{
    p=fitted.values(glm(catholic~math8*faminc8,family=binomial))
    ps=ifelse(catholic==1,sum(catholic)/length(catholic),1-sum(catholic)/length(catholic))
    eotm=ifelse(catholic==1,1/p,1/(1-p))
    seotm=ps*eotm
  })
  mod=lm(math12~catholic+faminc8*math8,data=dat,weights=dat$seotm)
  return(coef(mod)[2])
}

biptws=boot(cath,iptws,1000,strata=cath$catholic,parallel="multicore",ncpus=8)
c(mean(biptws$t),sd(biptws$t)) ## Mean is 1.60, SD is 0.252
quantile(biptws$t,probs = c(0.025,0.975))

library(PSweight)
pform <- catholic~math8*faminc8
oform <- math12~math8*faminc8
sand <- PSweight(ps.formula = pform,yname = "math12",zname="catholic",data = cath,augmentation = TRUE,out.formula = oform,family = 'gaussian',weight ='IPW')
mod <- summary(sand)
mod$estimates ## Estimate is 1.60, SD is 0.246

## Make a figure to compare different bootstrap estimates
plot(density(hhest$t))
lines(density(biptw$t),col='red')
lines(density(biptwt$t),col='blue')
lines(density(biptws$t),col='green')
legend('topright',legend=c('GBoot','Boot','TBoot','SBoot'),lty=1,col=c('black','red','blue','green'),cex=0.8)


## Covariate Balance
## Use the cobalt package
library(cobalt)
x <- cbind(cath$faminc8,cath$math8)
col_w_smd(x,treat = cath$catholic)
col_w_smd(x,treat = cath$catholic,weights = cath$eotm) ## Only need to report this one
col_w_smd(x,treat = cath$catholic,weights = cath$teotm)
col_w_smd(x,treat = cath$catholic,weights = cath$seotm)

