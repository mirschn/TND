###Run Simulation Study with partial interference

mod_logreg<-function(dat,ret_CI=F){
  TNDmod<-(glm(Y~V+C+vaxpos,family=binomial(),data=dat)) #adjusted RR estimate
  RR<-exp(coef(TNDmod)[2])
  se<-summary(TNDmod)$coefficient[2,2]
  CI=c(0,0)
  if(ret_CI==T){
    confintmat<-confint(TNDmod)
    CI<-exp(confintmat[2,])
  }
  Q1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C,vaxpos=dat$vaxpos)),type="response")
  Q0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C,vaxpos=dat$vaxpos)),type="response")
  return(list(mod=TNDmod, est=RR,CI=CI,Q1=Q1,Q0=Q0,se=se))
}


mod_g<-function(dat){
  TNDmod_g<-glm(V~C+vaxpos,family=binomial(),data=dat)
  g1=predict(TNDmod_g,type="response")
  return(list(mod=TNDmod_g,g1=g1))
}

mod_g_control<-function(dat){
  TNDmod_g<-glm(V~C+vaxpos,family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y,vaxpos=dat$vaxpos)));
  return(list(mod=TNDmod_g,g1_cont=g1_cont))
}



bootstrap<-function(dat,nbs=500,method1,method2=NA,method3=NA){
  n=dim(dat)[1]
  bsest1<-rep(NA,nbs)
  bsest2<-rep(NA,nbs)
  bsest3<-rep(NA,nbs)
  CI1=CI2=CI3=NA
  if(is.na(method2)&is.na(method3)){
    for(i in 1:nbs){
      resamps<-sample(1:n,size=n,replace=T)
      datk<-dat[resamps,]
      bsest1[i]<-method1(datk)
    }
  }
  if(!is.na(method2)&is.na(method3)){
    for(i in 1:nbs){
      resamps<-sample(1:n,size=n,replace=T)
      datk<-dat[resamps,]
      bsest1[i]<-method1(datk)
      bsest2[i]<-method2(datk)
    }
  }
  if(!is.na(method2)&!is.na(method3)){
    for(i in 1:nbs){
      resamps<-sample(1:n,size=n,replace=T)
      datk<-dat[resamps,]
      bsest1[i]<-method1(datk)
      bsest2[i]<-method2(datk)
      bsest3[i]<-method3(datk)
    }
  }
  CI1=quantile(bsest1,c(0.025,0.975))
  if(!is.na(method2)) CI2=quantile(bsest2,c(0.025,0.975))
  if(!is.na(method3)) CI3=quantile(bsest3,c(0.025,0.975))
  return(list(bs_var1=var(bsest1),CI1=CI1,bs_vect1=bsest1,
              bs_var2=var(bsest2),CI2=CI2,bs_vect2=bsest2,
              bs_var3=var(bsest3),CI3=CI3,bs_vect3=bsest3))
}

mod_IPW<-function(dat){
  gmodTND_cont<-mod_g_control(dat)
  gTND_cont<-gmodTND_cont$g1_cont
  est=mean(dat$Y*dat$V/gTND_cont)/mean(dat$Y*(1-dat$V)/(1-gTND_cont))
  return(est)
}


mod_IPW_all<-function(dat){
  gmodTND<-mod_g(dat)
  gTND<-gmodTND$g1
  est=mean(dat$Y*dat$V/gTND)/mean(dat$Y*(1-dat$V)/(1-gTND))
  return(est)
}


CI1=CI2=CI3=CI4=c(0,0)
var1=var2=var3=var4=0
seeds<-read.table("seeds.txt", header=F)$V1
for (i in c(1:1000)){ #
  CI1=CI2=CI3=CI4=c(NA,NA)
  var1=var2=var3=var4=NA
  bsobj=NA
  
  dat<-datagen_int(seeds[i])
  
  #adjusted logistic regression, ignore f_m, adjusts for vaxpos
  #This was first analysis used, in paper
  logregmod<-mod_logreg(dat,ret_CI=T)
  est1<-logregmod$est 
  var1<-(logregmod$se)^2
  CI1<-logregmod$CI
  

  
  #IPW with controls to fit g, in paper
  est2<-mod_IPW(dat)
  try( { try(
    bsobj<-bootstrap(dat,nbs=500,method1=mod_IPW,method2=mod_IPW_all),silent=TRUE);
      var2<-bsobj$bs_var1;
      CI2<-bsobj$CI1;
      var3<-bsobj$bs_var2;
      CI3<-bsobj$CI2 } ,silent=TRUE)
 
  
  #IPW with all subject used to fit g, in paper
  est3<-mod_IPW_all(dat)

  
  write(c(i,est1,var1,CI1,est2,var2,CI2,est3,var3,CI3),file="Scenint_1_results.txt",ncolumns=25,append=T)
}

res<-read.table("Scenint_1_results.txt",header=F)
1-colMeans(res)

mean(0.142849<=res$V5 & 0.142849>=res$V4) #89.2
mean(0.1519916<=res$V9 & 0.1519916>=res$V8) #92.0
mean(0.1519916<=res$V13 & 0.1519916>=res$V12) #0.00


sd(1-res[,2])
sd(1-res[,6])
sd(1-res[,10])


