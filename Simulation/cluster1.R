## H-H estimator 
hh=function(d,i){
  mod=lm(y~w+x*z,data=d[i,],weights=d[i,]$hhw)
  return(coef(mod)[2]-1.677)
}


## Bootstrap IPTW without trimming
iptw=function(d,i){
  mod=lm(y~w+x*z,data=d[i,],weights=d[i,]$eotm)
  return(coef(mod)[2]-1.677)
}


## Bootstrap IPTW with trimming 
iptwt=function(d,i){
  mod=lm(y~w+x*z,data=d[i,],weights=d[i,]$teotm)
  return(coef(mod)[2]-1.677)
}


#####################################################################
## Simulation code: estimated propensity scores & ordinary weights ##
#####################################################################

## This version is for correctly specified propensity score model

mchh=function(mctimes,hhiter,bootiter,samplesize,rho,varv){
  library(PoisBinOrdNonNor)
  library(boot)
  library(mvtnorm)
  cath=read.csv('nels.csv')
  mat <- matrix(c(1,0.285,0.285,1),nrow = 2, byrow = T)
  pinc8 <- as.numeric(table(cath$faminc8)/sum(table(cath$faminc8)))
  m1 <- find.cor.mat.star(cor.mat = mat,no.pois = 0,no.bin = 0,no.ord = 1,no.nonn = 1,ord.list = list(pinc8),nonn.list = list(c(51.49,93.77,0.4078,-0.6815)))
  output <- mat.or.vec(mctimes,6)
  cov <- matrix(c(27.4,sqrt(27.4)*sqrt(varv)*rho,sqrt(27.4)*sqrt(varv)*rho,varv),nrow=2,byrow = T)
  for(i in 1:mctimes){
    s <- genPBONN(samplesize,cmat.star = m1,no.ord = 1,no.nonn = 1,ord.list = list(pinc8),nonn.list = list(c(51.49,93.77,0.4078,-0.6815)))
    r <- rmvnorm(samplesize,mean=c(0,0),sigma=cov)
    rnum <- data.frame(cbind(s,r))
    colnames(rnum)=c("z","x","u","v")
    rnum=within(rnum,{
      lc <- 0.19*z+0.047*x-0.004*x*z+0.01*(z^2)+v-4.26
      w=ifelse(lc>0,1,0)
    })
    rnum=within(rnum,{
      p=fitted.values(glm(w~x*z+I(z^2),family=binomial))
      y=2.15+1.677*w+0.9*x+0.946*z-0.013*x*z+u
    })
    rnum=within(rnum,{
      eotm=ifelse(w==1,1/p,1/(1-p))
      teotm=ifelse(eotm>20,20,eotm)
    })
    rnum=within(rnum,{
      hhp=ifelse(w==1,eotm/sum(subset(rnum,w==1,eotm)),eotm/sum(subset(rnum,w==0,eotm)))
      hhw=ifelse(w==1,1/(nrow(rnum[rnum$w==1,])*hhp),1/(nrow(rnum[rnum$w==0,])*hhp))
    })
    hhest=boot(rnum,hh,hhiter,strata=rnum$w,weights=rnum$hhp,parallel="snow",ncpus=24)
    hhest=c(mean(hhest$t),sd(hhest$t))
    ## Bootstrap IPTW
    biptw=boot(rnum,iptw,bootiter,strata=rnum$w,parallel="snow",ncpus=24)
    biptw=c(mean(biptw$t),sd(biptw$t))
    ## Bootstrap IPTW with trimmed weights
    tbiptw=boot(rnum,iptwt,bootiter,strata=rnum$w,parallel="snow",ncpus=24)
    tbiptw=c(mean(tbiptw$t),sd(tbiptw$t))
    output[i,]=c(hhest,biptw,tbiptw)
  }
  output <- data.frame(output)
  colnames(output)=c("HHB","HHS","BootB","BootS","TBootB","TBootS")
  return(output)
}

mchh <- compiler::cmpfun(mchh)

##############################
## Simulation Design: 3*2*2 ##
##############################

## s1: sample size is 1000, s2: sample size is 5000, s3: sample size is 10000
## r1: rho is 0, r2: rho is 0.1
## m1: model with less oversize weights (v=1); m2: model with more oversize weights (v=0.3)


s1r1m1=mchh(1000,1000,1000,1000,0,1)
saveRDS(s1r1m1,'s1r1m1.rds')
s1r2m1=mchh(1000,1000,1000,1000,0.1,1)
saveRDS(s1r2m1,'s1r2m1.rds')
s2r1m1=mchh(1000,1000,1000,5000,0,1)
saveRDS(s2r1m1,'s2r1m1.rds')
s2r2m1=mchh(1000,1000,1000,5000,0.1,1)
saveRDS(s2r2m1,'s2r2m1.rds')
s3r1m1=mchh(1000,1000,1000,10000,0,1)
saveRDS(s3r1m1,'s3r1m1.rds')
s3r2m1=mchh(1000,1000,1000,10000,0.1,1)
saveRDS(s3r2m1,'s3r2m1.rds')
s1r1m2=mchh(1000,1000,1000,1000,0,0.3)
saveRDS(s1r1m2,'s1r1m2.rds')
s1r2m2=mchh(1000,1000,1000,1000,0.1,0.3)
saveRDS(s1r2m2,'s1r2m2.rds')
s2r1m2=mchh(1000,1000,1000,5000,0,0.3)
saveRDS(s2r1m2,'s2r1m2.rds')
s2r2m2=mchh(1000,1000,1000,5000,0.1,0.3)
saveRDS(s2r2m2,'s2r2m2.rds')
s3r1m2=mchh(1000,1000,1000,10000,0,0.3)
saveRDS(s3r1m2,'s3r1m2.rds')
s3r2m2=mchh(1000,1000,1000,10000,0.1,0.3)
saveRDS(s3r2m2,'s3r2m2.rds')


