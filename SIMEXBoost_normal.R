library(SIMEXBoost)
library(MASS)

betahat = NULL
for(K in 1:100) {

X1 = matrix(rnorm((100)*400),nrow=400,ncol=100,byrow=TRUE)
data=ME_Data(X1,beta=c(1,1,1,rep(0,dim(X1)[2]-3)),pr0=0.5,
type="AFT-normal",
sigmae=diag(0.1,dim(X1)[2]))
Y = data$response
Xstar = data$ME_covariate

B1 = SIMEXBoost(Y,Xstar,B=50,zeta=c(0,0.25,0.5,0.75,1),
type="AFT-normal",Iter=20,sigmae=diag(0.1,dim(X1)[2]))$BetaHatCorrect


betahat = rbind(betahat,B1)

}

colMeans(betahat)

SIMEXBoostBeta = colMeans(betahat)
SIMEXBoostBeta[which(abs(SIMEXBoostBeta)<0.02)]=0




