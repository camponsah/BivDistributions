#' library(distr)



DP.EM <- function(N,gamma1,delta, K)
{
  # initialization
  n <- length(N)
  Z<-rgamma(n,shape = 1/delta, rate = gamma1)
  
  like <- (1-exp(-Z))*exp(-Z*N)
  deviance <- -2*sum(log(like))
  res <- matrix(NA,K+1,4)
  res[1,] <- c(0, gamma1,delta, deviance)
  
  for (k in 1:K) {
    # S step
    a<- 1/(1+((N-1)/gamma1))
    b<- 1/(1+(N/gamma1))
    P<-a^(1/delta)-b^(1/delta)
    C1<- (1/(delta*P))*(gamma1^(1/delta))* ((gamma1+N-1)^(-(1/delta +1))-(gamma1+N)^(-(1/delta +1)))
    t1<- digamma(1/delta)-log(gamma1+N-1)
    t2<- digamma(1/delta)-log(gamma1+N)
    C2<- (1/P)*(gamma1^(1/delta))* (t1*((gamma1+N-1)^(-1/delta))-t2*((gamma1+N)^(-1/delta)))
    
    # M step
    oldgamma1<- gamma1
    gamma1<- 1/(delta*mean(C1))
    delta<- 1/distr:: igamma(log(gamma1)+mean(C2))
    
    # -2 x LL
    like <- (1-exp(-Z))*exp(-Z*N)
    deviance <- -2*sum(log(like))
    
    # add results to output
    res[k+1,] <- c(k, gamma1,delta, deviance)
    Z<-rgamma(n,shape = 1/delta, rate = gamma1)
  }
  res <- data.frame(res)
  names(res) <- c("iteration","gamma","delta","deviance")
  out <- list(parameters=c(q,beta, p), deviance=deviance, res=res)
  return(out)
}

p=0.5; delta=0.5; gamma1=-1/(delta*log(1-p)) ; K=100
n<-100
u<-runif(n)
N<- ceiling(gamma1*((1-u)^(-delta) -1))
fit <- DP.EM(N,0.5,0.1, K)
fit$res[K,]



