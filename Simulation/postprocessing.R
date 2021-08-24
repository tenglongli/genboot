
###############################################################
## Jointly simulate ordinal variable and continuous variable ##
###############################################################
cath=read.csv('nels.csv')
library(PoisBinOrdNonNor)
mat <- matrix(c(1,0.285,0.285,1),nrow = 2, byrow = T)
pinc8 <- as.numeric(table(cath$faminc8)/sum(table(cath$faminc8)))
m1 <- find.cor.mat.star(cor.mat = mat,no.pois = 0,no.bin = 0,no.ord = 1,no.nonn = 1,ord.list = list(pinc8),nonn.list = list(c(51.49,93.77,0.4078,-0.6815)))


##############################
## Simulation Design: 3*2*2 ##
##############################

## s1: sample size is 1000, s2: sample size is 5000, s3: sample size is 10000
## r1: rho is 0, r2: rho is 0.1
## m1: model with less oversize weights; m2: model with more oversize weights 
## v is 1 for m1 and 0.3 for m2

## Save the output (for mis-specified propensity score model only)

d <- readRDS('s1r1m1.rds')
d1 <- readRDS('s1r1m2.rds')
d2 <- readRDS('s1r2m1.rds')
d3 <- readRDS('s1r2m2.rds')
d4 <- readRDS('s2r1m1.rds')
d5 <- readRDS('s2r1m2.rds')
d6 <- readRDS('s2r2m1.rds')
d7 <- readRDS('s2r2m2.rds')
d8 <- readRDS('s3r1m1.rds')
d9 <- readRDS('s3r1m2.rds')
d10 <- readRDS('s3r2m1.rds')
d11 <- readRDS('s3r2m2.rds')

out <- list(d,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11)
saveRDS(out,'out.rds')

## out is the file for mis-specified propensity score models
## out1 is the file for correctly specified propensity score models
## In the main text, we only show results for the mis-specified propensity score models; Results in out1 are shown in appendix.

## Create a post-processing function 
postprocess <- function(d){
  ## 1-generate the mean bias estimates
  m <- apply(d,2,mean)
  m1=m[c(1,3,5)]
  ## 2-Find out the approximate true standard error 
  v=apply(d,2,function(x)sum(x^2)/1000)
  v1=v[c(1,3,5)]
  ts=sqrt(v1-m1^2)
  ## 3-Mean of the standard error estimates
  mse=m[c(2,4,6)]
  ## 4-generate SE of the standard error estimates
  sse=apply(d[,c(2,4,6)],2,sd)
  ## Extract the standard error estimates
  hhse=d[,2]
  iptwse=d[,4]
  iptwtse=d[,6]
  ## 5-generate the proportion of underestimation for s.e. estimates
  ppust <- c(sum(hhse<ts[1])/1000,sum(iptwse<ts[2])/1000,sum(iptwtse<ts[3])/1000)
  ## 6-generate the coverage rate
  ## Coverage rate for H-H Estimator
  hh=d[,1]+1.677
  hhl=hh-1.96*hhse
  hhu=hh+1.96*hhse
  hhcov=as.data.frame(cbind(hhl,hhu))
  ## Coverage rate for Bootstrap IPTW
  iptw=d[,3]+1.677
  iptwl=iptw-1.96*iptwse
  iptwu=iptw+1.96*iptwse
  iptwcov=as.data.frame(cbind(iptwl,iptwu))
  ## Coverage rate for Bootstrap IPTW with trimmed weights
  iptwt=d[,5]+1.677
  iptwtl=iptwt-1.96*iptwtse
  iptwtu=iptwt+1.96*iptwtse
  iptwtcov=as.data.frame(cbind(iptwtl,iptwtu))
  coverage <- c(nrow(hhcov[which(hhl<1.677&hhu>1.677),])/1000,nrow(iptwcov[which(iptwl<1.677&iptwu>1.677),])/1000,nrow(iptwtcov[which(iptwtl<1.677&iptwtu>1.677),])/1000)
  final_result <- rbind(m1,ts,mse,sse,ppust,coverage)
  rownames(final_result) <- c('Bias','TSE','Mean of SE','SE of SE','Pct. of UnderEST','Coverage Rate')
  colnames(final_result) <- c('HH','IPTW','IPTW_Trimmed')
  return(final_result)
}


postprocess <- compiler::cmpfun(postprocess)


## Clean results (for mis-specified propensity score model)
out <- readRDS('out.rds')
lapply(out,postprocess)


post <- function(data){
  dat <- data[,c(1,3,5)]
  colnames(dat) <- c('GB','OB','TB')
  m <- apply(dat,2,function(x) quantile(x,probs = c(0.025,0.975)))
  m <- as.data.frame(t(m))
  colnames(m) <- c('lower','upper')
  m$estimator <- row.names(m)
  m$mb <- apply(dat,2,mean)
  return(m)
}

## Create data for visualization 
library(tidyverse)

## For mis-specified propensity score model only
bias1 <- lapply(out,post)
bias1 <- do.call(rbind,bias1)
bias1$scenario <- rep(c('Missing Confounders-No; Oversize Weights-Less','Missing Confounders-No; Oversize Weights-More',
                        'Missing Confounders-Yes; Oversize Weights-Less','Missing Confounders-Yes; Oversize Weights-More'),each=3,times=3)
bias1$size <- rep(c(1000,5000,10000),each=12)
bias1$size <- as.factor(bias1$size)

## figure 1
ggplot(data = bias1)+geom_linerange(mapping = aes(x=size,ymin=lower,ymax=upper,group=estimator),position = position_dodge(width = 0.4))+
  geom_point(mapping = aes(x=size,y=mb,shape = estimator),position = position_dodge(width = 0.4),size = 2)+geom_hline(yintercept = 0,linetype="dashed")+
  labs(x = "Sample Size", y = "Bias",shape = "Mean Bias" )+facet_wrap(~scenario, nrow = 2)

## For comparing standard error estimates
## Still the same plot except that the dot now should represent the true standard error
## First create a data frame that contains the true standard error
post1 <- function(data){
  tse <- data[,c(1,3,5)]
  tse <- apply(tse,2,function(x) sqrt(sum(x^2)/1000-(mean(x))^2))
  dat <- data[,c(2,4,6)]
  colnames(dat) <- c('GB','OB','TB')
  m <- apply(dat,2,function(x) quantile(x,probs = c(0.025,0.975)))
  m <- as.data.frame(t(m))
  colnames(m) <- c('lower','upper')
  m$estimator <- row.names(m)
  m$tse <- tse
  return(m)
}

## Again, for mis-specified propensity score model only
se1 <- lapply(out,post1)
se1 <- do.call(rbind,se1)
se1$scenario <- rep(c('Missing Confounders-No; Oversize Weights-Less','Missing Confounders-No; Oversize Weights-More',
                        'Missing Confounders-Yes; Oversize Weights-Less','Missing Confounders-Yes; Oversize Weights-More'),each=3,times=3)
se1$size <- rep(c(1000,5000,10000),each=12)
se1$size <- as.factor(se1$size)

## figure 2
ggplot(data = se1)+geom_linerange(mapping = aes(x=size,ymin=lower,ymax=upper,group=estimator),position = position_dodge(width = 0.4))+geom_point(mapping = aes(x=size,y=tse,shape = estimator),position = position_dodge(width = 0.4),size=2)+
  labs(x = "Sample Size", y = "Standard Error",shape = "True SE" )+facet_wrap(~scenario, nrow = 2)



###################################################################
## Generate figure 0: comparing the distribution of IPTW weights ##
###################################################################

hh=function(samplesize,rho,varv){
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

hh <- compiler::cmpfun(hh)

w1 <- hh(10000,0,1)
w2 <- hh(10000,0,0.3)

boxplot(w1,w2,names=c("Less Oversize Weights","More Oversize Weights"))
abline(h=20,lty=2)
