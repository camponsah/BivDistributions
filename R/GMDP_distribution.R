
##### Simulate bivariate gamma mixture and discrete Pareto random variables
# n is number of observations
rgammamixdispareto<- function(n,alpha,beta,delta,p){
  u<-runif(n)
  sigma=-1/(delta*log(1-p))
  N<-ceiling(sigma*((1-u)^(-delta) -1))
  X<-rgamma(n,shape =alpha*N,rate = beta )
  return(data.frame(X,N))
}


# data is bivariate vector  (X,N) vector representing observations from BEG model
dgammamixdispareto<- function(data,alpha,beta,delta,p){
  N<-data[,2]
  X<-data[,1]
  p1<-(beta^(alpha*N))*(X^(alpha*N-1))*exp(-beta*X)/gamma(alpha*N)
  p21<-(1-delta*(N-1)*log(1-p))^(-1/delta)
  p22<-(1-delta*N*log(1-p))^(-1/delta)
  M<-p1*(p21-p22)
  return(M)
}

# q is bivariate vector  (X,N) vector quantiles from BEG model
pgammamixdispareto<- function(q,alpha,beta,delta,p, lower.tail=TRUE,log.p=FALSE){
  N<-q[,2]
  X<-q[,1]
  M<-NULL
  t=1
  for (i in 1:length(N)) {
    k=seq(1,N[i])
    S0<-0
    for (j in k) {
      p21<-(1-delta*(N[i]-1)*log(1-p))^(-1/delta)
      p22<-(1-delta*N[i]*log(1-p))^(-1/delta)
      S0<-S0 + (p21-p22)*pracma::gammainc(beta*X[i],j*alpha)[3] 
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

## Examples
Data.df<-rgammamixdispareto(n=10, alpha = 1,beta = 2, delta = 0.1,p=0.4)
pgammamixdispareto(Data.df,alpha = 1,beta = 2, delta = 0.1,p=0.4)
