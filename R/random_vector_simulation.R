
##### Simulate bivariate Lomax and geometric random variables
rblomaxgeo<- function(n,alpha,beta,p){
  N<-rgeom(n,p)+1
  X<-rgamma(n,shape =N,rate = beta )
  X<-X/rgamma(1,shape =1/alpha,rate = 1/alpha )
  return(data.frame(X,N))
}



