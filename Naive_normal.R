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

B1 = Boost_VSE(Y,Xstar,type="AFT-normal",Iter=20)$BetaHat


betahat = rbind(betahat,B1)

}

colMeans(betahat)

NaiveBeta = colMeans(betahat)




