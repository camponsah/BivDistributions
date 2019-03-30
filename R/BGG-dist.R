

##### Simulate bivariate gamma mixture and geometric random variables

# n is number of observations
rbgammageo<- function(n,alpha,beta,p){
  N<-rgeom(n,p)+1
  X<-rgamma(n,shape =alpha*N,rate = beta )
  return(data.frame(X,N))
}

# data is bivariate vector  (X,N) vector representing observations from bivariate gamma mixture and geometric model
dbgammageo<- function(data,alpha,beta,p){
  N<-data[,2]
  X<-data[,1]
  M<-((1-p)^(N-1))*(p*beta^(alpha*N))*(X^(alpha*N-1))*exp(-beta*X)/gamma(alpha*N)
  return(M)
}

# q is bivariate vector  (X,N) vector quantiles from bivariate gamma mixture and geometric model
pbgammageo<- function(q,alpha,beta,p, lower.tail=TRUE,log.p=FALSE){
  N<-q[,2]
  X<-q[,1]
  M<-NULL
  t=1
  for (i in 1:length(N)) {
    k=seq(1,N[i])
    S0<-0
    for (j in k) {
      S0<-S0 + p*((1-p)^(j-1))*pracma::gammainc(beta*X[i],j*alpha)[3] 
    }
    M[t]<- S0
    t=t+1
  }
  if (lower.tail==FALSE & log.p==TRUE){
    M<-log(1-M)
    return(M)
  }
  else if (lower.tail==TRUE & log.p==TRUE){
    M<-log(M)
    return(M) 
  }
  else if (lower.tail==TRUE & log.p==FALSE){
    return(M) 
  }
  else {
    return(1-M)
  }
}

Data.df<-rbgammageo(10, alpha = 2, beta = 100, p=0.1)
pbgammageo(Data.df,alpha = 2, beta = 100, p=0.1)
