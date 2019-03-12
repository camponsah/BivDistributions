rm(list = ls())
#' 
##### Simulate discrete pareto random variable
rdpareto<-function(n,delta,p){
   u<-runif(n)
  sigma=-1/(delta*log(1-p))
  return(ceiling(sigma*((1-u)^(-delta) -1)))
}


#' Algorithmn
EM.dpareto <- function(data,delta, p, maxiter=2000, tol=1e-16)
{
  # initialization
  N<-data
  n <- length(N)
  sigma=-1/(delta*log(1-p))
  ##Likelihood function of Discrete Pareto
  like<-function(N,sigma, delta){
    LL<- (1/(1+((N-1)/sigma)))^(1/delta)-(1/(1+((N)/sigma)))^(1/delta)
    return(LL)
  }
  ### function to optimize for delta
  func_delta<-function(delta){
    sig<- -1/(delta*mean(C1)) 
    ll<- (n/delta)*log(sig) - n*lgamma(1/delta) - sum((sig-1+N)*C1) +((1/delta)-1)*sum(C2)
    return(-sum(log(ll)))
  }
  Devianceold<-0
  Deviancenew <- sum(log(like(N,sigma = sigma,delta = delta)))
  Outi<-NULL;outd<-NULL;outp<-NULL;outD<-NULL; k=1 
  Outi[1]<-0; outd[1]<-delta; outp[1]<-p; outD<-Deviancenew
  ### log-like function
  Derivative<-function(delta){
    dd<-log(1/delta)-digamma(1/delta)-log(mean(C1))+mean(C2)
    dd
  }
  while(abs(Deviancenew-Devianceold)>tol && k <= maxiter){ 
    ### E step
    a<- 1/(1+((N-1)/sigma))
    b<- 1/(1+(N/sigma))
    P<-a^(1/delta)-b^(1/delta)
    C1<- (1/(delta*P))*(sigma^(1/delta))* ((sigma+N-1)^(-(1/delta +1))-(sigma+N)^(-(1/delta +1)))
    t1<- digamma(1/delta)-log(sigma+N-1)
    t2<- digamma(1/delta)-log(sigma+N)
    C2<- (1/P)*(sigma^(1/delta))* (t1*((sigma+N-1)^(-1/delta))-t2*((sigma+N)^(-1/delta)))
    
    #### M step
    delta<- newtonRaphson(Derivative, delta)$root #nlm(func_delta, p=delta)$estimate#  
    sigma<- 1/(delta*mean(C1))
    p<-1- exp(-1/(sigma*delta))
    Devianceold<-Deviancenew
    Deviancenew <- sum(log(like(N,sigma = sigma,delta = delta)))
    # Output
    k=k+1
    Outi[k]<-k; outd[k]<-delta; outp[1]<-p; outD<-Deviancenew
  }
  Output <- data.frame(Outi,outd,outp,outD)
  names(Output) <- c("iteration","delta","p","log-lik values")
  result <- list(par=c(delta,p), Deviance=Deviancenew, data=Output)
  return(result)
}


N<-rdpareto(100,delta = 2,p=0.6)

fit <- EM.dpareto(N,delta = 1, p=0.5)
fit$par
tail(fit$data)




