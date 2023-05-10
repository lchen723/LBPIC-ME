install.packages("bayesSurv")
install.packages("mice")
install.packages("DMwR2")
library(bayesSurv)
library(mice);library(DMwR2)
data(tandmob2, package="bayesSurv")
#?tandmob2
AA<-tandmob2[-c(which(is.na(tandmob2[,31]))),]
A<-AA[-c(which(is.na(AA[,32]))),]
data<-A[,c(-1,-2,-4,-(21:30),-(33:132))]
#interval<-which(is.na(data[,19]))
C<-mice(data,m=3,method = "pmm",maxit = 10,seed=163)
data1<-complete(C)
data1<-data1[-c(2001:(length(A[,1]))),]
Y<-cbind(data1[,18],data1[,19])

data2<-data1[,c(-(18:19))]
dataF<-centralImputation(cbind(Y,data2))
# for(i in 1:length(dataF[,1])){
#   for(j in 1:43){
#     if(dataF[i,j]==0) dataF[i,j]=-1
#     
#   }
# }
#dataF<-dataF[c(1:10),]
TtutaSet<-NULL;dataF1=dataF;
for(kfold in 1:5){
  n=length(dataF1[,1]);p=41;betazero=t(rep(0,p));namita=0.03
  dataF2=dataF1[-c((((kfold-1)*(n/5))+1):(((kfold)*(n/5)))),]
  #step2 想法，把x也分为两部分,分别计算出确切观测到的Y，以及censorY，用L、R、u来估计。
  #X1<-X
  LR<-cbind(dataF2[,1],dataF2[,2])
  #以下是SIMEXboost部分
  betarbind<-NULL;eta<-c(0,0.25,0.5,0.75,1);B=5
  for(cc in 1:length(eta)){
    betaihat<-NULL
    for(i in 1:B){
      X1<-t(dataF2[,c(3:43)])
      Xchange1<-X1[c(1,(5:41)),]
      orthogonal<-matrix(c(sqrt(0.5),sqrt(0.5),sqrt(0.5),-sqrt(0.5)),nrow=2,ncol=2,byrow=TRUE)
      orthogonalT<-t(orthogonal)
      piepsilon<-orthogonalT%*%(diag(c(1^(eta[cc]),0.9^(eta[cc]))))%*%orthogonal
      Xchange1pSet<-NULL
      for(s in 1:length(Xchange1[,1])){
        exchange1<-runif(length(X1[1,]),0,1)
        change1=which(exchange1<piepsilon[2])
        Xchange1p<-Xchange1[s,]
        #写一个1：p的回圈，把38个变数包在里面，每一个都改
        if(length(change1)!=0) for(o in 1:length(change1)){
          if(Xchange1[s,][change1[o]]==(0)) Xchange1p[change1[o]]=1
          else Xchange1p[change1[o]]=(0)
        }
        Xchange1pSet<-rbind(Xchange1pSet,Xchange1p)
      }
      # Xchange2<-X1[c(2:4),]
      # Xchange2pSet<-NULL
      # for(s in 1:length(Xchange2[,1])){
      #   exchange2<-runif(n,0,1)
      #   change2=which(exchange2<piepsilon[2])
      #   Xchange2p<-Xchange2[s,]
      #   #写一个1：p的回圈，把38个变数包在里面，每一个都改
      #   if(length(change2)!=0) for(o in 1:length(change2)){
      #     if(Xchange2[s,][change2[o]]==(0)) Xchange2p[change2[o]]=sample(1:2,1)
      #     else Xchange2p[change2[o]]=(0)
      #   }
      #   Xchange2pSet<-rbind(Xchange2pSet,Xchange2p)
      # }
      
      Wetab<-rbind(Xchange1pSet[1,],X1[c(2:4),],Xchange1pSet[c(2:38),])
      #Wetab1<-Wetab[,-l1out]
      X1censor<-Wetab
      TT1censor<-log(LR[,1]);TT1<-TT1censor
      functionx<-function(x){
        b=1/(exp(1)*(exp(1)-1))
        a=(log(x)/x)*(b/((pi)*(((log(x)/x)^2)+b^2)))
        return(a)
      }
      function2<-function(x){
        b=1/(exp(1)*(exp(1)-1))
        a=(b/((pi)*(((log(x)/x)^2)+b^2)))
        return(a)
      }
      Iter=50;time=0
      for(j in 1:Iter){ 
        betazero1=betazero
        Lbeta1<-NULL
        for(i in 1:length(TT1censor)){
          Lbeta<-(LR[,1][i]*exp(-sum(X1censor[,i]*betazero1)))
          Lbeta1<-cbind(Lbeta1,Lbeta)
        }
        Rbeta1<-NULL
        for(i in 1:length(TT1censor)){
          Rbeta<-(LR[,2][i]*exp(-sum(X1censor[,i]*betazero1)))
          Rbeta1<-cbind(Rbeta1,Rbeta)
        }
        dFtgeneration<-NULL
        for(i in 1:length(TT1censor)){
          a<-runif(1000,Lbeta1[i],Rbeta1[i])
          a1<-functionx(a)
          dFt=sum((Rbeta1[i]-Lbeta1[i])*a1)/length(a1)
          b<-runif(1000,Lbeta1[i],Rbeta1[i])
          b1<-function2(b)
          dFt1<-sum((Rbeta1[i]-Lbeta1[i])*b1)/length(b1)
          censorYi<-dFt/dFt1
          dFtgeneration<-rbind(dFtgeneration,censorYi)
        }
        for(i in 1:length(dFtgeneration)){
          if(dFtgeneration[i]=="Inf") dFtgeneration[i]=0
          if(dFtgeneration[i]=="-Inf") dFtgeneration[i]=0
          if(dFtgeneration[i]=="NaN") dFtgeneration[i]=0
        }
        ubetastep2<-NULL
        for(i in 1:length(dFtgeneration)){
          a<-X1censor[,i]*dFtgeneration[i]-namita*betazero1
          ubetastep2<-rbind(ubetastep2,a)
        }
        #compute u
        ubetastep<-ubetastep2
        u<-NULL
        for(i in 1:p){
          u<-c(u,sum(ubetastep[,i])/length(TT1))
        }
        for(i in 1:p){
          Delta=u
          #if(abs(Delta[i])>0.6*max(abs(Delta))) betazero[i]=betazero[i]+0.01*sign(Delta[i])#Delta[i]=Delta[i]+0.05*sign(Delta[i])
          if(abs(Delta[i])>0*max(abs(Delta))) betazero[i]=betazero[i]+0.02*sign(Delta[i])
          #if(abs(Delta[i])==max(abs(Delta))) betazero[i]=betazero[i]+0.005*Delta[i]
          #if(time==Iter/2) betazero[which(abs(betazero)<0.05)]=0
          if(time==Iter-1) betazero[which(abs(betazero)<0.02)]=0
        }
        time=time+1
      }
      betaihat<-rbind(betaihat,betazero)
    }
    betarbind<-rbind(betarbind,betaihat)
  }
  meanbeta0hat<-NULL;meanbeta0.25hat<-NULL;meanbeta0.5hat<-NULL;meanbeta0.75hat<-NULL;meanbeta1hat<-NULL
  beta0hat<-betarbind[1:B,];beta0.25hat<-betarbind[(B+1):(2*B),];
  beta0.5hat<-betarbind[((2*B)+1):(3*B),];beta0.75hat<-betarbind[((3*B)+1):(4*B),];
  beta1hat<-betarbind[((4*B)+1):(5*B),]
  for (i in 1:p){
    a<-mean(beta0hat[,i])
    b<-mean(beta0.25hat[,i])
    c<-mean(beta0.5hat[,i])
    d<-mean(beta0.75hat[,i])
    e<-mean(beta1hat[,i])
    meanbeta0hat<-cbind(meanbeta0hat,a)
    meanbeta0.25hat<-cbind(meanbeta0.25hat,b)
    meanbeta0.5hat<-cbind(meanbeta0.5hat,c)
    meanbeta0.75hat<-cbind(meanbeta0.75hat,d)
    meanbeta1hat<-cbind(meanbeta1hat,e)
  }
  # for(i in 1:p){
  #   if(abs(meanbeta0hat[i])<0.05) meanbeta0hat[i]=0
  #   if(abs(meanbeta0.25hat[i])<0.05) meanbeta0.25hat[i]=0
  #   if(abs(meanbeta0.5hat[i])<0.05) meanbeta0.5hat[i]=0
  #   if(abs(meanbeta0.75hat[i])<0.05) meanbeta0.75hat[i]=0
  #   if(abs(meanbeta1hat[i])<0.05) meanbeta1hat[i]=0
  # }
  gamma1hat<-NULL;gamma0hat<-NULL;#gamma2hat<-NULL;#gamma3hat<-NULL
  for(i in 1:p){
    x<-c(0,0.25,0.5,0.75,1)
    y<-c(meanbeta0hat[i],meanbeta0.25hat[i],
         meanbeta0.5hat[i],meanbeta0.75hat[i],meanbeta1hat[i])
    c<-lm(y~x)
    #c<-lm(y~x+I(x*x))
    gamma0<-as.numeric(c$coefficients[1])
    gamma1<-as.numeric(c$coefficients[2])
    #gamma2<-as.numeric(c$coefficients[3])
    #gamma3<-as.numeric(c$coefficients[4])
    
    gamma0hat<-c(gamma0hat,gamma0)
    gamma1hat<-c(gamma1hat,gamma1)
    #gamma2hat<-c(gamma2hat,gamma2)
  }
  betacorrect<-gamma0hat-gamma1hat
  #betacorrect<-gamma0hat-gamma1hat+gamma2hat
  for(i in 1:p){
    if(abs(betacorrect[i])<0.005) betacorrect[i]=0
  }
  kfoldprec<-dataF1[c((((kfold-1)*(n/5))+1):((kfold*(n/5)))),][3:41]
  for(kk in 1:length(kfoldprec[,1])){
    Ttuta<-sum(kfoldprec[kk,]*betacorrect)
    TtutaSet<-c(TtutaSet,Ttuta)
  }
  
}
TT<-data.frame(exp(TtutaSet))
#View(cbind(Y[1:435,],TT,TtutaSet))
Set<-NULL
for(i in 1:length(TT[,1])){
  if(Y[i,1]<TT[i,]&Y[i,2]>TT[i,]) Set=c(Set,TT[i,])
  
}
interval<-which(Y[,1]<TT&Y[,2]>TT)
# sum((dataF[1,][3:41])*betazero);betazero
pin<-length(which(Y[,1]<TT&Y[,2]>TT))/length(Y[,1])
TtutaSet1<-TtutaSet[-interval]
Y1<-Y[-interval,]
doutSet<-NULL
for(i in 1:length(TtutaSet1)){
  dout<-min(abs(Y1[i,1]-exp(TtutaSet1[i]))/(Y1[i,1]),abs(Y1[i,2]-exp(TtutaSet1[i]))/Y1[i,2])
  doutSet<-c(doutSet,dout)
}
mean(doutSet);1-pin



