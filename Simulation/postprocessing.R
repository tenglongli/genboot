###################
## Making Tables ##
###################

###################################################
## Simulation Design for Simulation Study 1: 3*2 ##
###################################################

## s1: sample size is 200, s2: sample size is 500, s3: sample size is 1000
## m1: model with less oversize weights; m2: model with more oversize weights 
## v is 1 for m1 and 0.3 for m2


s11 <- readRDS('s1m11.rds')
s12 <- readRDS('s1m21.rds')
s13 <- readRDS('s2m11.rds')
s14 <- readRDS('s2m21.rds')
s15 <- readRDS('s3m11.rds')
s16 <- readRDS('s3m21.rds')


out1 <- list(s11,s12,s13,s14,s15,s16)


## Create a post-processing function 
postprocess <- function(d){
  ## 1-extract the bias estimates
  b <- d[,seq(1,17,by=4)]
  m <- apply(b,2,mean) ## Mean bias estimates
  ## 2-find out the approximate true standard error 
  v=apply(b,2,function(x)sum(x^2)/1000)
  ts=sqrt(v-m^2)
  ## 3-extract the mean of standard error estimates
  se=d[,seq(2,18,by=4)]
  mse=apply(se,2,mean) ## the mean of standard error estimate
  ## 4-generate the ratio between the true SE and the estimated SE
  seratio=mse/ts ## ratio between the mean SE estimate and the true SE
  ## 5-generate the coverage rate based on the SE approach
  ## Coverage rate for H-H Estimator
  hhl=b[,1]-1.96*se[,1]
  hhu=b[,1]+1.96*se[,1]
  hhcov=as.data.frame(cbind(hhl,hhu))
  ## Coverage rate for Bootstrap IPTW
  iptwl=b[,2]-1.96*se[,2]
  iptwu=b[,2]+1.96*se[,2]
  iptwcov=as.data.frame(cbind(iptwl,iptwu))
  ## Coverage rate for Bootstrap IPTW with trimmed weights
  iptwtl=b[,3]-1.96*se[,3]
  iptwtu=b[,3]+1.96*se[,3]
  iptwtcov=as.data.frame(cbind(iptwtl,iptwtu))
  ## Coverage rate for Bootstrap IPTW with stablized weights
  iptwsl=b[,4]-1.96*se[,4]
  iptwsu=b[,4]+1.96*se[,4]
  iptwscov=as.data.frame(cbind(iptwsl,iptwsu))
  ## Coverage rate for Sandwich estimator
  sandl=b[,5]-1.96*se[,5]
  sandu=b[,5]+1.96*se[,5]
  sandcov=as.data.frame(cbind(sandl,sandu))
  coverage_se <- c(nrow(hhcov[which(hhl<0&hhu>0),])/1000,nrow(iptwcov[which(iptwl<0&iptwu>0),])/1000,nrow(iptwtcov[which(iptwtl<0&iptwtu>0),])/1000,nrow(iptwscov[which(iptwsl<0&iptwsu>0),])/1000,nrow(sandcov[which(sandl<0&sandu>0),])/1000)
  ## 6-generate the coverage rate based on the percentile approach
  ## Coverage rate for H-H Estimator
  hhl=d[,3]
  hhu=d[,4]
  hhcov=as.data.frame(cbind(hhl,hhu))
  ## Coverage rate for Bootstrap IPTW
  iptwl=d[,7]
  iptwu=d[,8]
  iptwcov=as.data.frame(cbind(iptwl,iptwu))
  ## Coverage rate for Bootstrap IPTW with trimmed weights
  iptwtl=d[,11]
  iptwtu=d[,12]
  iptwtcov=as.data.frame(cbind(iptwtl,iptwtu))
  ## Coverage rate for Bootstrap IPTW with stablized weights
  iptwsl=d[,15]
  iptwsu=d[,16]
  iptwscov=as.data.frame(cbind(iptwsl,iptwsu))
  coverage_per <- c(nrow(hhcov[which(hhl<0&hhu>0),])/1000,nrow(iptwcov[which(iptwl<0&iptwu>0),])/1000,nrow(iptwtcov[which(iptwtl<0&iptwtu>0),])/1000,nrow(iptwscov[which(iptwsl<0&iptwsu>0),])/1000,nrow(sandcov[which(sandl<0&sandu>0),])/1000)
  final_result <- rbind(m,ts,mse,seratio,coverage_se,coverage_per,3.92*mse)
  rownames(final_result) <- c('Bias','TSE','Nominal SE','Ratio of SE','SE Coverage','Percentile Coverage','Width of CI')
  colnames(final_result) <- c('HH','Boot','TBoot','SBoot','Sandwich')
  return(final_result)
}


postprocess <- compiler::cmpfun(postprocess)

lapply(out1,postprocess)

###################################################
## Simulation Design for Simulation Study 2: 3*2 ##
###################################################

## s1: sample size is 200, s2: sample size is 500, s3: sample size is 1000
## m1: model with less oversize weights; m2: model with more oversize weights 
## v is 1 for m1 and 0.3 for m2

s21 <- readRDS('s1m12.rds')
s22 <- readRDS('s1m22.rds')
s23 <- readRDS('s2m12.rds')
s24 <- readRDS('s2m22.rds')
s25 <- readRDS('s3m12.rds')
s26 <- readRDS('s3m22.rds')


out2 <- list(s21,s22,s23,s24,s25,s26)
lapply(out2,postprocess)


####################
## Making Figures ##
####################


post <- function(data){
  dat <- data[,seq(1,17,by=4)]
  colnames(dat) <- c('GB','OB','TB','SB','SW')
  m <- apply(dat,2,function(x) quantile(x,probs = c(0.025,0.975)))
  m <- as.data.frame(t(m))
  colnames(m) <- c('lower','upper')
  m$estimator <- row.names(m)
  m$mb <- apply(dat,2,mean)
  return(m)
}

## Create data for visualization 
library(tidyverse)

## Appendix Figure 1
bias1 <- lapply(out1,post)
bias1 <- do.call(rbind,bias1)
bias1$scenario <- rep(c('Less Oversize Weights','More Oversize Weights'),each=5,times=3)
bias1$size <- rep(c(200,500,1000),each=10)
bias1$size <- as.factor(bias1$size)
ggplot(data = bias1)+geom_linerange(mapping = aes(x=size,ymin=lower,ymax=upper,group=estimator),position = position_dodge(width = 0.4))+
  geom_point(mapping = aes(x=size,y=mb,shape = estimator),position = position_dodge(width = 0.4),size = 2)+geom_hline(yintercept = 0,linetype="dotted")+
  labs(x = "Sample Size", y = "Bias",shape = "Mean Bias" )+facet_wrap(~scenario, nrow = 2)

## Appendix Figure 2
bias2 <- lapply(out2,post)
bias2 <- do.call(rbind,bias2)
bias2$scenario <- rep(c('Less Oversize Weights','More Oversize Weights'),each=5,times=3)
bias2$size <- rep(c(200,500,1000),each=10)
bias2$size <- as.factor(bias2$size)
ggplot(data = bias2)+geom_linerange(mapping = aes(x=size,ymin=lower,ymax=upper,group=estimator),position = position_dodge(width = 0.4))+
  geom_point(mapping = aes(x=size,y=mb,shape = estimator),position = position_dodge(width = 0.4),size = 2)+geom_hline(yintercept = 0,linetype="dotted")+
  labs(x = "Sample Size", y = "Bias",shape = "Mean Bias" )+facet_wrap(~scenario, nrow = 2)


## For comparing standard error estimates
## Still the same plot except that the dot now should represent the true standard error
## First create a data frame that contains the true standard error
post1 <- function(data){
  tse <- data[,seq(1,17,by=4)]
  tse <- apply(tse,2,function(x) sqrt(sum(x^2)/1000-(mean(x))^2))
  dat <- data[,seq(2,18,by=4)]
  colnames(dat) <- c('GB','OB','TB','SB','SW')
  m <- apply(dat,2,function(x) quantile(x,probs = c(0.025,0.975)))
  m <- as.data.frame(t(m))
  colnames(m) <- c('lower','upper')
  m$estimator <- row.names(m)
  m$tse <- tse
  return(m)
}

## Figure 1
se1 <- lapply(out1,post1)
se1 <- do.call(rbind,se1)
se1$scenario <- rep(c('Less Oversize Weights','More Oversize Weights'),each=5,times=3)
se1$size <- rep(c(200,500,1000),each=10)
se1$size <- as.factor(se1$size)
ggplot(data = se1)+geom_linerange(mapping = aes(x=size,ymin=lower,ymax=upper,group=estimator),position = position_dodge(width = 0.4))+geom_point(mapping = aes(x=size,y=tse,shape = estimator),position = position_dodge(width = 0.4),size=2)+
  labs(x = "Sample Size", y = "Standard Error",shape = "True SE" )+facet_wrap(~scenario, nrow = 2)


## Figure 2
se2 <- lapply(out2,post1)
se2 <- do.call(rbind,se2)
se2$scenario <- rep(c('Less Oversize Weights','More Oversize Weights'),each=5,times=3)
se2$size <- rep(c(200,500,1000),each=10)
se2$size <- as.factor(se2$size)
ggplot(data = se2)+geom_linerange(mapping = aes(x=size,ymin=lower,ymax=upper,group=estimator),position = position_dodge(width = 0.4))+geom_point(mapping = aes(x=size,y=tse,shape = estimator),position = position_dodge(width = 0.4),size=2)+
  labs(x = "Sample Size", y = "Standard Error",shape = "True SE" )+facet_wrap(~scenario, nrow = 2)


###############################################################
## Compare less oversized weights vs. more oversized weights ##
###############################################################

sim1=function(samplesize,rho,varv){
  library(PoisBinOrdNonNor)
  library(boot)
  library(mvtnorm)
  cath=read.csv('nels.csv')
  mat <- matrix(c(1,0.285,0.285,1),nrow = 2, byrow = T)
  pinc8 <- as.numeric(table(cath$faminc8)/sum(table(cath$faminc8)))
  m1 <- find.cor.mat.star(cor.mat = mat,no.pois = 0,no.bin = 0,no.ord = 1,no.nonn = 1,ord.list = list(pinc8),nonn.list = list(c(51.49,93.77,0.4078,-0.6815)))
  cov <- matrix(c(27.4,sqrt(27.4)*sqrt(varv)*rho,sqrt(27.4)*sqrt(varv)*rho,varv),nrow=2,byrow = T)
  s <- genPBONN(samplesize,cmat.star = m1,no.ord = 1,no.nonn = 1,ord.list = list(pinc8),nonn.list = list(c(51.49,93.77,0.4078,-0.6815)))
  r <- rmvnorm(samplesize,mean=c(0,0),sigma=cov)
  rnum <- data.frame(cbind(s,r))
  colnames(rnum)=c("z","x","u","v")
  rnum=within(rnum,{
    lc <- 0.19*z+0.047*x-0.004*x*z+0.01*(z^2)+v-4.26
    w=ifelse(lc>0,1,0)
  })
  rnum=within(rnum,{
    p=fitted.values(glm(w~x+z,family=binomial))
    y=2.15+1.677*w+0.9*x+0.946*z-0.013*x*z+u
  })
  rnum=within(rnum,{
    eotm=ifelse(w==1,1/p,1/(1-p))
    teotm=ifelse(eotm>20,20,eotm)
  })
  return(rnum$eotm)
}

sim1 <- compiler::cmpfun(sim1)

w1 <- sim1(10000,0,1)
w2 <- sim1(10000,0,0.3)

sim2=function(samplesize,varv){
  library(boot)
  library(PSweight)
  library(mvtnorm)
  cov <- matrix(c(2,1,0,0,1,1,0,0,0,0,1,0,0,0,0,varv),nrow=4,byrow = T)
  r <- rmvnorm(samplesize,mean=c(0.5,1,0,0),sigma=cov)
  rnum <- data.frame(r)
  colnames(rnum)=c("x","z","u","v")
  rnum=within(rnum,{
    lc <- 0.5+0.75*z+0.25*x+v
    w=ifelse(lc>0,1,0)
    y=1+w+x+2*z+u
  })
  rnum=within(rnum,{
    p=fitted.values(glm(w~x+z,family=binomial))
    eotm=ifelse(w==1,1/p,1/(1-p))
  })
  return(rnum$eotm)
}

sim2 <- compiler::cmpfun(sim2)

w1 <- sim2(10000,1)
w2 <- sim2(10000,0.3)
