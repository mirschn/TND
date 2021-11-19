#thinking about what to do in response to reviews

#adjusting for community-level vaccine coverage
#maybe take one larger sample from each block?
popsize=5000000
seeds<-read.table("seeds.txt", header=F)$V1
dat<-datagen_int(seeds[1],popsize=popsize,return_full=T) #note: increased popsize to 5million so can adjust for block
datTND=as.data.frame(cbind(Y=dat$Infec_COVID,V=dat$V,C=dat$C,block=dat$block,f_m=dat$f_m,vaxpos=dat$vaxpos)[dat$H==1,])

mod_logreg<-function(dat,ret_CI=T){
  TNDmod<-(glm(Y~V+C,family=binomial(),data=dat)) #adjusted RR estimate
  RR<-exp(coef(TNDmod)[2])
  se<-summary(TNDmod)$coefficient[2,2]
  CI=c(0,0)
  if(ret_CI==T){
    confintmat<-confint(TNDmod)
    CI<-exp(confintmat[2,])
  }
  #Q1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C)),type="response")
  #Q0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C)),type="response")
  #return(list(mod=TNDmod, est=RR,CI=CI,Q1=Q1,Q0=Q0,se=se))
  return(list(est=(1-RR),CI=(1-CI)))
}
mod_g_control<-function(dat){
  TNDmod_g<-glm(V~C,family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y)));
  return(list(mod=TNDmod_g,g1_cont=g1_cont))
}
library(sandwich)
mod_IPW<-function(dat){
  gmodTND_cont<-mod_g_control(dat)
  gTND_cont<-gmodTND_cont$g1_cont
  w=dat$V/gTND_cont + (1-dat$V)/(1-gTND_cont)
  ipwglm=glm(Y~V,data=dat,family=binomial,weights = w)
  est=coefficients(ipwglm)[2]
  CI = c(est-1.96*sqrt(vcovHC(ipwglm)[2,2]) , est+1.96*sqrt(vcovHC(ipwglm)[2,2]) )
  est=1-exp(est)
  CI=1-exp(CI)
  return(list(est,CI))
}

mod_IPW.b<-function(dat){
  gmodTND_cont<-mod_g_control(dat)
  gTND_cont<-gmodTND_cont$g1_cont
  gTND_cont[gTND_cont<quantile(gTND_cont,0.10)]=quantile(gTND_cont,0.10)
  w=dat$V/gTND_cont + (1-dat$V)/(1-gTND_cont)
  ipwglm=glm(Y~V,data=dat,family=binomial,weights = w)
  est=coefficients(ipwglm)[2]
  CI = c(est-1.96*sqrt(vcovHC(ipwglm)[2,2]) , est+1.96*sqrt(vcovHC(ipwglm)[2,2]) )
  est=1-exp(est)
  CI=1-exp(CI)
  return(list(est,CI))
}

by(datTND$f_m,datTND$block,mean)

by(datTND$f_m,datTND$block,length)

by(datTND$f_m[datTND$Y==0],datTND$block[datTND$Y==0],length)

by(datTND,datTND$block,mod_logreg)

by(datTND,datTND$block,mod_IPW)

by(datTND,datTND$block,mod_IPW.b)

#estimation for block 2 is unstable because of tiny number of vaccinated controls
summary(mod_g_control(datTND[datTND$block==2,])$mod)
sum(datTND$block==2&datTND$Y==0&datTND$V==0) #26
sum(datTND$block==2&datTND$Y==0&datTND$V==1) #5
gTND_cont=(mod_g_control(datTND[datTND$block==2,])$g1_cont)
summary(datTND$V[datTND$block==2]/gTND_cont)

#note that 10% truncation was used in the paper -- negligeable impact except for Block 2, which had unstable estimates.

#size of each block
by(dat,dat$block,dim)


#Now overall estimates
mod_IPW_bl<-function(dat){
  gmodTND_cont<-mod_g_control_bl(dat)
  gTND_cont<-gmodTND_cont$g1_cont
  w=dat$V/gTND_cont + (1-dat$V)/(1-gTND_cont)
  ipwglm=glm(Y~V,data=dat,family=binomial,weights = w)
  est=coefficients(ipwglm)[2]
  CI = c(est-1.96*sqrt(vcovHC(ipwglm)[2,2]) , est+1.96*sqrt(vcovHC(ipwglm)[2,2]) )
  est=1-exp(est)
  CI=1-exp(CI)
  return(list(est,CI))
}


mod_logreg_bl<-function(dat,ret_CI=T){
  TNDmod<-(glm(Y~V+C+factor(block),family=binomial(),data=dat)) #adjusted RR estimate, removed vaxpos
  RR<-1-exp(coef(TNDmod)[2])
  se<-summary(TNDmod)$coefficient[2,2]
  CI=c(0,0)
  if(ret_CI==T){
    confintmat<-confint(TNDmod)
    CI<-1-exp(confintmat[2,])
  }
  Q1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C,block=dat$block)),type="response")
  Q0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C,block=dat$block)),type="response")
  return(list(est=RR,CI=CI))
}

mod_g_control<-function(dat){
  TNDmod_g<-glm(V~C+vaxpos,family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y,vaxpos=dat$vaxpos)));
  return(list(mod=TNDmod_g,g1_cont=g1_cont))
}

mod_g_control_bl<-function(dat){
  TNDmod_g<-glm(V~C+factor(block),family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y,block=dat$block)));
  return(list(mod=TNDmod_g,g1_cont=g1_cont))
}

#MVLR adjusting and not adjust for block
mod_logreg(datTND)
mod_logreg_bl(datTND)

#IPTW adjusting and not adjust for block
mod_IPW(datTND)
mod_IPW_bl(datTND)
