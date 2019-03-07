#' library(distr)



DP.EM <- function(N,delta, p, K)
{
  # initialization
  n <- length(N)
  gamma1=-1/(delta*log(1-p))
  Z<-rgamma(n,shape = 1/delta, rate = gamma1)
  Deviance <- -2*sum(log((1-exp(-Z))*exp(-Z*N)))
  Output <- matrix(NA,K+1,4)
  Output[1,] <- c(0,delta, p, Deviance)
  
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
    p<-1- exp(-1/(gamma1*delta))
    # Output
    Output[k+1,] <- c(k, delta, p, Deviance)
    Z<-rgamma(n,shape = 1/delta, rate = gamma1)
  }
  Output <- data.frame(Output)
  names(Output) <- c("iteration","delta","p","Deviance")
  result <- list(par=c(delta,p), Deviance=Deviance, data=Output)
  return(result)
}

p=0.1; delta=0.2; K=100
n<-1000
u<-runif(n)
gamma1=-1/(delta*log(1-p))
N<- ceiling(gamma1*((1-u)^(-delta) -1))
fit <- DP.EM(N,0.5,0.5, K)
fit$par



