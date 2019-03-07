#' library(distr)



DP.EM <- function(N,gamma1,delta, K)
{
  # initialization
  n <- length(N)
  Z<-rgamma(n,shape = 1/delta, rate = gamma1)
  Deviance <- -2*sum(log((1-exp(-Z))*exp(-Z*N)))
  Output <- matrix(NA,K+1,4)
  Output[1,] <- c(0, gamma1,delta, Deviance)
  
  for (k in 1:K) {
    ### E step
    a<- 1/(1+((N-1)/gamma1))
    b<- 1/(1+(N/gamma1))
    P<-a^(1/delta)-b^(1/delta)
    C1<- (1/(delta*P))*(gamma1^(1/delta))* ((gamma1+N-1)^(-(1/delta +1))-(gamma1+N)^(-(1/delta +1)))
    t1<- digamma(1/delta)-log(gamma1+N-1)
    t2<- digamma(1/delta)-log(gamma1+N)
    C2<- (1/P)*(gamma1^(1/delta))* (t1*((gamma1+N-1)^(-1/delta))-t2*((gamma1+N)^(-1/delta)))
    
    #### M step
    oldgamma1<- gamma1
    gamma1<- 1/(delta*mean(C1))
    delta<- 1/distr:: igamma(log(gamma1)+mean(C2))
    
    Deviance <- -2*sum(log((1-exp(-Z))*exp(-Z*N)))
    
    # Output
    Output[k+1,] <- c(k, gamma1,delta, Deviance)
    Z<-rgamma(n,shape = 1/delta, rate = gamma1)
  }
  Output <- data.frame(Output)
  names(Output) <- c("iteration","gamma","delta","Deviance")
  result <- list(par=c(gamma1,delta), Deviance=Deviance, data=Output)
  return(result)
}

p=0.5; delta=0.5; gamma1=-1/(delta*log(1-p)) ; K=100
n<-100
u<-runif(n)
N<- ceiling(gamma1*((1-u)^(-delta) -1))
fit <- DP.EM(N,0.5,0.1, K)
fit$par



