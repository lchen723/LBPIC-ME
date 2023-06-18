library(mice);library(DMwR2);library(survival);library(survminer)

## Note: From the datasets, responses for tooth label 11 are columns 21 and 22. Following an order in the main text, labels are 21, 31, 41, and so on, and the corresponding columns are (23,24),(25,26),(27,28), and so on.

## keep tooth label 11 from columns 21 and 22 and remove NA to get interval-censored responses
AA<-tandmob2[-c(which(is.na(tandmob2[,21]))),] 
A<-AA[-c(which(is.na(AA[,22]))),]

## remove the first three non-informative columns and other teeth labels
data<-A[,c(-1,-2,-4,-(23:32),-(33:132))]  
C<-mice(data,m=3,method = "pmm",maxit = 10,seed=163)
data1<-complete(C)
data1<-data1[-c(1741:(length(A[,1]))),]
Y<-cbind(data1[,18],data1[,19])       
data2<-data1[,c(-(18:19))]


dataF<-centralImputation(cbind(Y,data2))

TtutaSet<-NULL;dataF1=dataF;

## 5-fold cross validation
for(kfold in 1:5){
  n=length(dataF1[,1]);p=41;betazero=t(rep(0,p));namita=0.03
  dataF2=dataF1[-c((((kfold-1)*(n/5))+1):(((kfold)*(n/5)))),]
  #step2 Ie¡P¡L¢G?¢XNxO2¡POIaA?2?¡PO,¡PO¡Óe?AEa3oE¡PCD1U2a£g?£gAY¢G?OO?¢XcensorY¢G?OAL!¢FR!¢FuA¡¦1A?A!¢G
  #X1<-X
  LR<-cbind(dataF2[,1],dataF2[,2])
  #OOIAECSIMEXboost2?¡PO
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
        #D¡¦O?¡Mo1¢Gop£gA?OE|¢G?¢XN38¡Mo¡ÓaEy¢XuOUAiAa¢G?A?O?¡Mo??¡MA
        if(length(change1)!=0) for(o in 1:length(change1)){
          if(Xchange1[s,][change1[o]]==(0)) Xchange1p[change1[o]]=1
          else Xchange1p[change1[o]]=(0)
        }
        Xchange1pSet<-rbind(Xchange1pSet,Xchange1p)
      }

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

  gamma1hat<-NULL;gamma0hat<-NULL;#gamma2hat<-NULL;#gamma3hat<-NULL
  for(i in 1:p){
    x<-c(0,0.25,0.5,0.75,1)
    y<-c(meanbeta0hat[i],meanbeta0.25hat[i],
         meanbeta0.5hat[i],meanbeta0.75hat[i],meanbeta1hat[i])
    c<-lm(y~x)
    #c<-lm(y~x+I(x*x))
    gamma0<-as.numeric(c$coefficients[1])
    gamma1<-as.numeric(c$coefficients[2])

    
    gamma0hat<-c(gamma0hat,gamma0)
    gamma1hat<-c(gamma1hat,gamma1)

  }
  betacorrect<-gamma0hat-gamma1hat

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

Set<-NULL
for(i in 1:length(TT[,1])){
  if(Y[i,1]<TT[i,]&Y[i,2]>TT[i,]) Set=c(Set,TT[i,])
  
}
interval<-which(Y[,1]<TT&Y[,2]>TT)

pin<-length(which(Y[,1]<TT&Y[,2]>TT))/length(Y[,1])
TtutaSet1<-TtutaSet[-interval]
Y1<-Y[-interval,]
doutSet<-NULL
for(i in 1:length(TtutaSet1)){
  dout<-min(abs(Y1[i,1]-exp(TtutaSet1[i]))/(Y1[i,1]),abs(Y1[i,2]-exp(TtutaSet1[i]))/Y1[i,2])
  doutSet<-c(doutSet,dout)
}
mean(doutSet);1-pin



