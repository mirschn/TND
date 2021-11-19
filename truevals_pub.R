##Scenario 1 

#true conditional RR for hospitalization due to COVID
orvect<-rep(NA,50)
for (j in 1:50){
  datfull<-datagen(popsize=1.5*10^6,return_full=T)
  Qmod<-glm(I(Infec_COVID*W*H)~V+C,data=datfull,family=binomial)
  orvect[j]<-exp(Qmod$coef[2])
}
mean(orvect,na.rm=T) #0.04
median(orvect) #
hist(orvect,nclass=20)

#true marginal RR
orvect2<-rep(NA,50)
for (j in 1:50){
  datfull1<-datagen(cfV1=T,return_full=T)
  datfull0<-datagen(cfV0=T,return_full=T)
  orvect2[j]<-mean(datfull1$H*datfull1$Infec_COVID)/mean(datfull0$H*datfull0$Infec_COVID)
}
mean(orvect2) #0.04
hist(orvect2)

#true marginal RR for Infec_COVID
orvect2<-rep(NA,50)
for (j in 1:50){
  datfull1<-datagen(cfV1=T,return_full=T)
  datfull0<-datagen(cfV0=T,return_full=T)
  orvect2[j]<-mean(datfull1$Infec_COVID)/mean(datfull0$Infec_COVID)
}
mean(orvect2) #0.16
hist(orvect2)

#true marginal RR for Infec_COVID * W (severe COVID disease)
orvect2<-rep(NA,50)
for (j in 1:50){
  datfull1<-datagen(cfV1=T,return_full=T)
  datfull0<-datagen(cfV0=T,return_full=T)
  orvect2[j]<-mean(datfull1$W*datfull1$Infec_COVID)/mean(datfull0$W*datfull0$Infec_COVID)
}
mean(orvect2) #0.04
hist(orvect2)



##Scenario 2 

#true conditional RR
orvect<-rep(NA,50)
for (j in 1:50){
  datfull<-datagen(popsize=1.5*10^6,return_full=T,em=1)
  Qmod<-glm(I(Infec_COVID*W*H)~V+C,data=datfull,family=binomial)
  orvect[j]<-exp(Qmod$coef[2]) 
}
mean(orvect,na.rm=T) #0.23
median(orvect) #
hist(orvect,nclass=20)

#true marginal RR
orvect2<-rep(NA,50)
for (j in 1:50){
  datfull1<-datagen(cfV1=T,return_full=T,em=1)
  datfull0<-datagen(cfV0=T,return_full=T,em=1)
  orvect2[j]<-mean(datfull1$H*datfull1$Infec_COVID)/mean(datfull0$H*datfull0$Infec_COVID)
}
mean(orvect2) #0.25
hist(orvect2)

#true marginal RR for Infec_COVID
orvect2<-rep(NA,50)
for (j in 1:50){
  datfull1<-datagen(cfV1=T,return_full=T,em=1)
  datfull0<-datagen(cfV0=T,return_full=T,em=1)
  orvect2[j]<-mean(datfull1$Infec_COVID)/mean(datfull0$Infec_COVID)
}
mean(orvect2) #0.56
hist(orvect2)

#true marginal RR for Infec_COVID * W
orvect2<-rep(NA,50)
for (j in 1:50){
  datfull1<-datagen(cfV1=T,return_full=T,em=1)
  datfull0<-datagen(cfV0=T,return_full=T,em=1)
  orvect2[j]<-mean(datfull1$W*datfull1$Infec_COVID)/mean(datfull0$W*datfull0$Infec_COVID)
}
mean(orvect2) #0.23
hist(orvect2)

## Scenatio 2 again -- subgroup effects

#true marginal RR in subset where C in -inf, -0.674 (Q1)
orvect1<-rep(NA,50)
orvect2<-rep(NA,50)
orvect3<-rep(NA,50)
orvect4<-rep(NA,50)
for (j in 1:50){
  datfull1<-datagen(cfV1=T,return_full=T,em=1)
  datfull0<-datagen(cfV0=T,return_full=T,em=1)
  orvect1[j]<-mean(datfull1$H[datfull1$C<(-0.6745)]*datfull1$Infec_COVID[datfull1$C<(-0.6745)])/mean(datfull0$H[datfull0$C<(-0.6745)]*datfull0$Infec_COVID[datfull0$C<(-0.6745)])
  orvect2[j]<-mean(datfull1$H[datfull1$C>=(-0.6745) & datfull1$C<(0)]*datfull1$Infec_COVID[datfull1$C>=(-0.6745) & datfull1$C<(0)])/mean(datfull0$H[datfull0$C>=(-0.6745) & datfull0$C<(0)]*datfull0$Infec_COVID[datfull0$C>=(-0.6745) & datfull0$C<(0)])
  orvect3[j]<-mean(datfull1$H[datfull1$C>=(0) & datfull1$C<(0.6745)]*datfull1$Infec_COVID[datfull1$C>=(0) & datfull1$C<(0.6745)])/mean(datfull0$H[datfull0$C>=(0) & datfull0$C<(0.6745)]*datfull0$Infec_COVID[datfull0$C>=(0) & datfull0$C<(0.6745)])
  orvect4[j]<-mean(datfull1$H[datfull1$C>=(0.6745)]*datfull1$Infec_COVID[datfull1$C>=(0.6745)])/mean(datfull0$H[datfull0$C>=(0.6745)]*datfull0$Infec_COVID[datfull0$C>=(0.6745)])
  
  }
1-mean(orvect1)
1-median(orvect1)# 
hist(1-orvect1)

#95th percentile of C is about 1.65
1-mean(datfull1$H[datfull1$C>=(1.648105)]*datfull1$Infec_COVID[datfull1$C>=(1.648105)])/mean(datfull0$H[datfull0$C>=(1.648105)]*datfull0$Infec_COVID[datfull0$C>=(1.648105)])


#true marginal RR in subset where C in -0.674, 0 (Q2)
1-mean(orvect2)
1-median(orvect2)# 
hist(1-orvect2)

#true marginal RR in subset where C in 0, 0.674 (Q3)
1-mean(orvect3)
1-median(orvect3)# 
hist(1-orvect3)

#true marginal RR in subset where C in 0.674+ (Q4)
1-mean(orvect4)
1-median(orvect4)# 
hist(1-orvect4)


##Scenario 3

#direct effect of vaccination, holding f_m fixed

mrrvect<-rep(NA,50)
for (j in 1:50){
  datfull1<-datagen_int_cf(cfV1=T,f_m_val=0.77,ssize_m=2000000)
  datfull0<-datagen_int_cf(cfV0=T,f_m_val=0.77,ssize_m=2000000)
  mrrvect[j]=1-mean(datfull1$H*datfull1$Infec_COVID)/mean(datfull0$H*datfull0$Infec_COVID) 
}
hist(mrrvect)
mean(mrrvect)
median(mrrvect)
#0.75 median=mean=0.86
#0.5 median=mean=0.8
#0.25 mean=median=0.76
#0.37 mean=median=0.77
#0.24 mean=median=0.756
#0.32 mean=median=0.77
#0.69 mean=median = 0.84
#0.58 mean=median = 0.82
#0.79 mean=median= 0.87
#0.83 mean=median= 0.88
#0.78 mean=median= 0.87
#0.77 mean=median=0.87


#true conditional RR 
#true mRR for infection
rrvect<-rep(NA,50)
mrrvect<-rep(NA,50)
for (j in 1:50){
  datfull<-datagen_int(popsize=1.5*10^6,return_full=T)
  Qmod<-glm(I(Infec_COVID)~V+C+vaxpos,data=datfull,family=binomial(link="logit"))#rare outcome, so same as rr
  rrvect[j]<-exp(Qmod$coef[2]) #0.52 #1-0.52=0.48 
  Q1<-predict(Qmod,newdata=as.data.frame(cbind(C=C,V=1)),type="response")
  Q0<-predict(Qmod,newdata=as.data.frame(cbind(C=C,V=0)),type="response")
  mrrvect[j]<-mean(Q1)/mean(Q0) #0.15# 1-0.15=0.85
}
mean(rrvect,na.rm=T) #0.142849
median(rrvect) #0.1416554
mean(mrrvect,na.rm=T) #0.1532558
median(mrrvect) #0.1519916
hist(mrrvect,nclass=10) #looks fine for both