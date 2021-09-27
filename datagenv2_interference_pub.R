#Common data-generating function

#popsize set at 1500000, can increase if larger ssize is needed
#create 10 clusters, each will be recruited from their own hospital. The vaccine prev in the cluster can affect individuals' probability of infection (only)
#ssize_m set at 100 from each cluster = 1000

#Defaults true marg RR = 
#OR_C<-3 #effect of vaccine on covid
#OR_WI<-1 #no effect of vaccine on W from other infection
#OR_WC<-5 #effect of vaccine on covid symptoms
#OR_H<-1.5 #effect of vaccine on hospitalization among those with symptoms
#em=0
#vaxpos is a study-level confounder so need to adjust for this


datagen_int<-function(seed=sample(1:1000000,size=1),ssize_m=100,popsize=1500000,OR_C=3,OR_WI=1,OR_WC=5,OR_H=1.5,em=0,W_spectrum=F,cfV0=F,cfV1=F,f_m_val=NA,return_full=F){
  set.seed(seed)
  
  #blocks are predefined
  block<-c(rep(1,floor(popsize*0.05)),rep(2,floor(popsize*0.05)),rep(3,floor(popsize*0.05)),rep(4,floor(popsize*0.05)),rep(5,floor(popsize*0.10)),rep(6,floor(popsize*0.10)),rep(7,floor(popsize*0.10)),rep(8,floor(popsize*0.10)),rep(9,floor(popsize*0.20)),rep(10,floor(popsize*0.20)))
  vaxpos.b<-abs(rnorm(10)) # Needs to affect V and thus f_m
  vaxpos<-c(rep(vaxpos.b[1],floor(popsize*0.05)),rep(vaxpos.b[2],floor(popsize*0.05)),rep(vaxpos.b[3],floor(popsize*0.05)),rep(vaxpos.b[4],floor(popsize*0.05)),rep(vaxpos.b[5],floor(popsize*0.10)),rep(vaxpos.b[6],floor(popsize*0.10)),rep(vaxpos.b[7],floor(popsize*0.10)),rep(vaxpos.b[8],floor(popsize*0.10)),rep(vaxpos.b[9],floor(popsize*0.20)),rep(vaxpos.b[10],floor(popsize*0.20)))
  
  
  #generate data (WITH clustering for U2)
  C<-rnorm(n=popsize) 
  U1<-rbinom(n=popsize,size=1,prob=0.5) #affects both 
  U2<-rbinom(n=popsize,size=1,prob=plogis(0.4+0.2*vaxpos)) #affects covid #shifted upwards for vaxpos blocks #it's like reckless behaviour specific to COVID
  

  V<-rbinom(prob=plogis(0.3*C-1.3*(block<=3)+(block>=6)+0.5*vaxpos),size=1,n=popsize) #prevalence is around 0.68
  #by(V,block,summary)
  
  f_m<-by(V,block,mean)
  f_m<-c(rep(f_m[1],floor(popsize*0.05)),rep(f_m[2],floor(popsize*0.05)),rep(f_m[3],floor(popsize*0.05)),rep(f_m[4],floor(popsize*0.05)),rep(f_m[5],floor(popsize*0.10)),rep(f_m[6],floor(popsize*0.10)),rep(f_m[7],floor(popsize*0.10)),rep(f_m[8],floor(popsize*0.10)),rep(f_m[9],floor(popsize*0.20)),rep(f_m[10],floor(popsize*0.20)))

  
  #Infection (with something) has some common risk factors U1 and C
  Infec<-rbinom(prob=plogis(0.2+0.5*C-5+0.5*U1),size=1,n=popsize) #current prevalence around 0.007
  
  #Infected with COVID #vaccine more effective with high vaccination rates
  Infec_COVID<- rbinom(prob=plogis( -3.5-vaxpos-1.2*log(OR_C)*V*(f_m)^2  + C+em*V*C+log(2)*U2-2*U1), size=1,n=popsize) #0.018
  #mean(Infec_COVID)
  
  #symptoms based on infection
  
  W=W1=W2=rep(0,popsize)

    W1[Infec==1]<-rbinom(prob=plogis(-4+0.5*C[Infec==1]-log(OR_WI)*V[Infec==1]-0.5*U1[Infec==1]),size=1, n=sum(Infec==1))
    W2[Infec_COVID==1]<-rbinom(prob=plogis(-1.5+1*C[Infec_COVID==1]-log(OR_WC)*V[Infec_COVID==1]-1*U1[Infec_COVID==1]+0.5*U2[Infec_COVID==1]*(1-V[Infec_COVID==1])),size=1, n=sum(Infec_COVID))
    W=(W1|W2)+0
    #mean(W[Infec==1|Infec_COVID==1]) #13%
    #mean(W[Infec_COVID==1]) #25%
    #mean(W[Infec==1]) #3%
    #mean(W[Infec_COVID==1&V==1]) #14%
    
    #hospitalization
    H=rep(0,popsize)
    H[W==1]<-rbinom(prob=plogis(0.5*C[W==1]-log(OR_H)*V[W==1]-0.5*U1[W==1]+f_m[W==1]),size=1,n=sum(W==1))
    #mean(H[W==1]) #77% with severe symptoms go to hospital
    #get idea of true values of cond RR
    #qpglm<-glm(H*Infec_COVID~C+V,family=quasipoisson) 1-exp(coef)=77%
  

  
  #selection on outcome for testing (does not condition on infectious status, just being in the hospital)
  R<-c(sample((1:popsize)[block==1&H==1],min(ssize_m,sum(block==1&H==1))), #sample randomly from each hospital to control the study size
       sample((1:popsize)[block==2&H==1],min(ssize_m,sum(block==2&H==1))), 
       sample((1:popsize)[block==3&H==1],min(ssize_m,sum(block==3&H==1))), 
       sample((1:popsize)[block==4&H==1],min(ssize_m,sum(block==4&H==1))), 
       sample((1:popsize)[block==5&H==1],min(ssize_m,sum(block==5&H==1))), 
       sample((1:popsize)[block==6&H==1],min(ssize_m,sum(block==6&H==1))), 
       sample((1:popsize)[block==7&H==1],min(ssize_m,sum(block==7&H==1))), 
       sample((1:popsize)[block==8&H==1],min(ssize_m,sum(block==8&H==1))),
       sample((1:popsize)[block==9&H==1],min(ssize_m,sum(block==9&H==1))), 
       sample((1:popsize)[block==10&H==1],min(ssize_m,sum(block==10&H==1))) )
       
  
  if(return_full==F){
    dat<-as.data.frame(cbind(Y=Infec_COVID,V=V,C=C,block=block,f_m=f_m,vaxpos=vaxpos)[R,])
  } else{dat<-as.data.frame(cbind(Infec_COVID=Infec_COVID,Infec=Infec,H=H,W=W,V=V,C=C,block=block,f_m=f_m,vaxpos=vaxpos))}
  return(dat)
}

datagen_int_cf<-function(seed=sample(1:1000000,size=1),ssize_m=1000000,OR_C=3,OR_WI=1,OR_WC=5,OR_H=1.5,em=0,W_spectrum=F,cfV0=F,cfV1=F,f_m_val=NA){
  set.seed(seed)
  popsize=ssize_m
  #assign person 1 to 0 or 1. Assign the rest of their cluster according to Bern(f_m_val). But instead of repeating, just do it all at once since others don't have impact beyond f_m
  #blocks are predefined
  block<-c(rep(1,floor(popsize*0.05)),rep(2,floor(popsize*0.05)),rep(3,floor(popsize*0.05)),rep(4,floor(popsize*0.05)),rep(5,floor(popsize*0.10)),rep(6,floor(popsize*0.10)),rep(7,floor(popsize*0.10)),rep(8,floor(popsize*0.10)),rep(9,floor(popsize*0.20)),rep(10,floor(popsize*0.20)))
  vaxpos.b<-abs(rnorm(10)) # Needs to affect V and thus f_m
  vaxpos<-c(rep(vaxpos.b[1],floor(popsize*0.05)),rep(vaxpos.b[2],floor(popsize*0.05)),rep(vaxpos.b[3],floor(popsize*0.05)),rep(vaxpos.b[4],floor(popsize*0.05)),rep(vaxpos.b[5],floor(popsize*0.10)),rep(vaxpos.b[6],floor(popsize*0.10)),rep(vaxpos.b[7],floor(popsize*0.10)),rep(vaxpos.b[8],floor(popsize*0.10)),rep(vaxpos.b[9],floor(popsize*0.20)),rep(vaxpos.b[10],floor(popsize*0.20)))
  
  
  #generate data (WITH clustering for U2)
  C<-rnorm(n=popsize) #shifted upwards for vaxpos blocks
  U1<-rbinom(n=popsize,size=1,prob=0.5) #affects both 
  U2<-rbinom(n=popsize,size=1,prob=plogis(0.4+0.2*vaxpos)) #affects covid #it's like reckless behaviour specific to COVID
  
  
  V=rep(NA,ssize_m)
  if(cfV0==T){V=rep(0,ssize_m)}; if(cfV1==T){V=rep(1,ssize_m)}
  f_m=f_m_val #between 0 and 1
  
  #Infection (with something) has some common risk factors U1 and C
  Infec<-rbinom(prob=plogis(0.2+0.5*C-5+0.5*U1),size=1,n=popsize) #current prevalence around 0.007
  
  #Infected with COVID
  Infec_COVID<- rbinom(prob=plogis( -3.5-vaxpos-1.2*log(OR_C)*V*(f_m)^2  + C+em*V*C+log(2)*U2-2*U1), size=1,n=popsize) #0.013
  #mean(Infec_COVID)
  
  #symptoms based on infection
  
  W=W1=W2=rep(0,popsize)
  
  W1[Infec==1]<-rbinom(prob=plogis(-4+0.5*C[Infec==1]-log(OR_WI)*V[Infec==1]-0.5*U1[Infec==1]),size=1, n=sum(Infec==1))
  W2[Infec_COVID==1]<-rbinom(prob=plogis(-1.5+1*C[Infec_COVID==1]-log(OR_WC)*V[Infec_COVID==1]-1*U1[Infec_COVID==1]+0.5*U2[Infec_COVID==1]*(1-V[Infec_COVID==1])),size=1, n=sum(Infec_COVID))
  W=(W1|W2)+0
  #mean(W[Infec==1|Infec_COVID==1]) #13%
  #mean(W[Infec_COVID==1]) #25%
  #mean(W[Infec==1]) #3%
  #mean(W[Infec_COVID==1&V==1]) #14%
  
  #hospitalization
  H=rep(0,popsize)
  H[W==1]<-rbinom(prob=plogis(0.5*C[W==1]-log(OR_H)*V[W==1]-0.5*U1[W==1]+f_m),size=1,n=sum(W==1))
  
  dat<-as.data.frame(cbind(Infec_COVID=Infec_COVID,Infec=Infec,H=H,W=W,V=V,C=C,vaxpos=vaxpos))
  return(dat)
}


