

## n is number of observation
rdpareto<-function(n,delta,p){
  u<- runif(n)
  sigma<- - 1/(delta*log(1 - p))
  return(ceiling(sigma*((1 - u)^(- delta) - 1)))
}



# N is vector of observations from discrete pareto distribution
ddpareto<-function(N,delta,p,log.p=FALSE){
 W1<- (1 - delta*(N - 1)*log(1 - p))^( - 1/delta)
 W2<- (1 - delta*N*log(1 - p))^( - 1/delta)
 M<- W1 - W2
 if (log.p == FALSE){
   return(M)
   }
 else{
   M<- log(M)
   return(M)
   }
}

# q is vector of quantiles values from discrete pareto distribution
pdpareto<- function(q,delta,p,lower.tail=TRUE,log.p=FALSE){
  sigma = -1/(delta*log(1 - p))
  M<- 1- ((1 + q/sigma)^(- 1/delta))
if (lower.tail == TRUE & log.p == FALSE){
  return(M)
  }
else if (lower.tail == TRUE & log.p == TRUE){
  M<- log(M)
  return(M)
  }
else if (lower.tail == FALSE & log.p == TRUE){
  M<- log(1 - M)
  return(M)
  }
else {
  return(1 - M)
  }
}


# prob is vector of probabilities
qdpareto<-function(prob,delta,p){
  sigma <- - 1/(delta*log(1 - p))
  return(ceiling(sigma*((1-  prob)^(- delta) - 1)))
}

#Examples
rdpareto(n=100,delta = 0.1,p=0.5)
qdpareto(seq(0.1,0.9,0.1), delta = 2,p=0.1)
