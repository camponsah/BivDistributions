##### Simulate bivariate exponential and geometric random variables

# n is number of observations
rbexpgeo<- function(n,beta,p){
  N<-rgeom(n,p)+1
  X<-rgamma(n,shape =N,rate = beta )
  return(data.frame(X,N))
}


# data is bivariate vector  (X,N) vector representing observations from BEG model
dbexpgeo<- function(data,beta,p){
  N<-data[,2]
  X<-data[,1]
  M<-(p*beta^N)*((X*(1-p))^(N-1))*exp(-beta*X)/gamma(N)
  return(M)
}


# q is bivariate vector  (X,N) vector quantiles from BEG model
pbexpgeo<- function(q,beta,p, lower.tail=TRUE,log.p=FALSE){
  N<-q[,2]
  X<-q[,1]
  S1<-ppois(N-1,lambda = ((1-p)*beta*X),lower.tail = TRUE)
  S2<-ppois(N-1,lambda = (beta*X),lower.tail = FALSE)
  M<-1-exp(-p*beta*X)*S1 -((1-p)^N)*S2
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


Data.df<-rbexpgeo(10, beta = 100, p=0.1)
pbexpgeo(Data.df, beta = 100, p=0.1)
