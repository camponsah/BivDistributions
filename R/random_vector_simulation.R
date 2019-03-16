
#'
##### Simulate discrete pareto random variable
rdpareto<-function(n,delta,p){
  u<-runif(n)
  sigma=-1/(delta*log(1-p))
  return(ceiling(sigma*((1-u)^(-delta) -1)))
}

##### Simulate bivariate exponential and geometric random variables
rbeg<- function(n,beta,p){
  N<-rgeom(n,p)+1
  X<-rgamma(n,shape =N,rate = beta )
  return(data.frame(X,N))
}


##### Simulate bivariate gamma mixture and geometric random variables
rbgg<- function(n,alpha,beta,p){
  N<-rgeom(n,p)+1
  X<-rgamma(n,shape =alpha*N,rate = beta )
  return(data.frame(X,N))
}


##### Simulate bivariate gamma mixture and discrete Pareto random variables
rgmdp<- function(n,alpha,beta,delta,p){
  u<-runif(n)
  sigma=-1/(delta*log(1-p))
  N<-ceiling(sigma*((1-u)^(-delta) -1))
  X<-rgamma(n,shape =alpha*N,rate = beta )
  return(data.frame(X,N))
}


##### Simulate bivariate Lomax and geometric random variables
rblg<- function(n,alpha,beta,p){
  N<-rgeom(n,p)+1
  X<-rgamma(n,shape =N,rate = beta )/rgamma(n,shape =1/alpha,rate = 1/alpha )
  return(data.frame(X,N))
}


