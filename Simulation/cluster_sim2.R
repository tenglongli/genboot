## H-H estimator 
hh=function(d,i){
  mod=lm(y~w+x+z,data=d[i,],weights=d[i,]$hhw)
  return(coef(mod)[2]-1)
}


## Bootstrap IPTW without trimming
iptw=function(d,i){
  dat=d[i,]
  dat=within(dat,{
    p=fitted.values(glm(w~x+z,family=binomial))
    eotm=ifelse(w==1,1/p,1/(1-p))
  })
  mod=lm(y~w+x+z,data=dat,weights=dat$eotm)
  return(coef(mod)[2]-1)
}


## Bootstrap IPTW with trimming 
iptwt=function(d,i){
  dat=d[i,]
  dat=within(dat,{
    p=fitted.values(glm(w~x+z,family=binomial))
    eotm=ifelse(w==1,1/p,1/(1-p))
    teotm=ifelse(eotm<=20,eotm,20)
  })
  mod=lm(y~w+x+z,data=dat,weights=dat$teotm)
  return(coef(mod)[2]-1)
}


## Bootstrap IPTW with stabilized weights
iptws=function(d,i){
  dat=d[i,]
  dat=within(dat,{
    p=fitted.values(glm(w~x+z,family=binomial))
    ps=ifelse(w==1,sum(w)/length(w),1-sum(w)/length(w))
    eotm=ifelse(w==1,1/p,1/(1-p))
    seotm=ps*eotm
  })
  mod=lm(y~w+x+z,data=dat,weights=dat$seotm)
  return(coef(mod)[2]-1)
}

################################################################
## Simulation Study 2: Design based on Freedman & Berk (2008) ##
################################################################

mchh=function(mctimes,hhiter,bootiter,samplesize,varv){
  library(boot)
  library(PSweight)
  library(mvtnorm)
  cov <- matrix(c(2,1,0,0,1,1,0,0,0,0,1,0,0,0,0,varv),nrow=4,byrow = T)
  output <- mat.or.vec(mctimes,18)
  for(i in 1:mctimes){
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
    rnum=within(rnum,{
      hhp=ifelse(w==1,eotm/sum(subset(rnum,w==1,eotm)),eotm/sum(subset(rnum,w==0,eotm)))
      hhw=ifelse(w==1,1/(nrow(rnum[rnum$w==1,])*hhp),1/(nrow(rnum[rnum$w==0,])*hhp))
    })
    hhest=boot(rnum,hh,hhiter,strata=rnum$w,weights=rnum$hhp,parallel="multicore",ncpus=36)
    hhest=c(mean(hhest$t,na.rm=TRUE),sd(hhest$t,na.rm=TRUE),quantile(hhest$t,probs=c(0.025,0.975),na.rm=TRUE))
    ## Bootstrap IPTW
    biptw=boot(rnum,iptw,bootiter,strata=rnum$w,parallel="multicore",ncpus=36)
    biptw=c(mean(biptw$t,na.rm=TRUE),sd(biptw$t,na.rm=TRUE),quantile(biptw$t,probs=c(0.025,0.975),na.rm=TRUE))
    ## Bootstrap IPTW with trimmed weights
    tbiptw=boot(rnum,iptwt,bootiter,strata=rnum$w,parallel="multicore",ncpus=36)
    tbiptw=c(mean(tbiptw$t,na.rm=TRUE),sd(tbiptw$t,na.rm=TRUE),quantile(tbiptw$t,probs=c(0.025,0.975),na.rm=TRUE))
    ## Bootstrap IPTW with stabilized weights
    sbiptw=boot(rnum,iptws,bootiter,strata=rnum$w,parallel="multicore",ncpus=36)
    sbiptw=c(mean(sbiptw$t,na.rm=TRUE),sd(sbiptw$t,na.rm=TRUE),quantile(sbiptw$t,probs=c(0.025,0.975),na.rm=TRUE))
    ## Robust sandwich estimates
    pform <- w~x+z
    oform <- y~x+z
    sand <- PSweight(ps.formula = pform,yname = "y",zname="w",data = rnum,augmentation = TRUE,out.formula = oform,family = 'gaussian',weight ='IPW')
    mod <- summary(sand)
    sandwich <- c(mod$estimates[1,1]-1,mod$estimates[1,2])
    output[i,]=c(hhest,biptw,tbiptw,sbiptw,sandwich)
  }
  output <- data.frame(output)
  colnames(output)=c("HHB","HHS","HHLB","HHUB","BootB","BootS","BootLB","BootUB",
                     "TBootB","TBootS","TBootLB","TBootUB","SBootB","SBootS","SBootLB","SBootUB",
                     "SandB","SandS")
  return(output)
}

mchh <- compiler::cmpfun(mchh)

############################
## Simulation Design: 3*2 ##
############################

## s1: sample size is 200, s2: sample size is 500, s3: sample size is 1000
## m1: model with less oversize weights (varv=1); m2: model with more oversize weights (varv=0.3)

s1m1=mchh(1000,1000,1000,200,1)
saveRDS(s1m1,'s1m12.rds')
s2m1=mchh(1000,1000,1000,500,1)
saveRDS(s2m1,'s2m12.rds')
s3m1=mchh(1000,1000,1000,1000,1)
saveRDS(s3m1,'s3m12.rds')
s1m2=mchh(1000,1000,1000,200,0.3)
saveRDS(s1m2,'s1m22.rds')
s2m2=mchh(1000,1000,1000,500,0.3)
saveRDS(s2m2,'s2m22.rds')
s3m2=mchh(1000,1000,1000,1000,0.3)
saveRDS(s3m2,'s3m22.rds')

