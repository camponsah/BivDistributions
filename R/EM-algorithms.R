#' 
#' Algorithmn
EM.dpareto <- function(data,delta, p, maxiter=400, tol=1e-12)
{
  # initialization
  N<-data
  n <- length(N)
  gamma1=-1/(delta*log(1-p))
  ##Likelihood function of Discrete Pareto
  like<-function(N,gamma1, delta){
    LL<- (1/(1+((N-1)/gamma1)))^(1/delta)-(1/(1+((N)/gamma1)))^(1/delta)
    return(sum(log(LL)))
  }
  ### function to optimize for delta
  func_delta<-function(omega){ # omega=1/delta
    ll<- n*omega*log(omega/mean(C1)) -n*lgamma(omega)-sum((omega/mean(C1)+N-1)*C1)+(omega-1)*sum(C2)
    return(-ll)
  }
  Devianceold<-0
  Deviancenew <- sum(log(like(N,gamma1 = gamma1,delta = delta)+tol)) 
  Outi<-NULL;outd<-NULL;outp<-NULL;outD<-NULL; k=1 
  Outi[1]<-0; outd[1]<-delta; outp[1]<-p; outD<-Deviancenew
  ### log-like function
  Derivative<-function(delta){
    dd<-log(1/delta)-digamma(1/delta)-log(mean(C1))+mean(C2)
    dd
  }
  while(abs(Deviancenew-Devianceold)>tol && k <= maxiter){ 
    ### E step
    a<- 1/(1+((N-1)/gamma1))
    b<- 1/(1+(N/gamma1))
    P<-a^(1/delta)-b^(1/delta)
    C1<- (1/(delta*P))*(gamma1^(1/delta))* ((gamma1+N-1)^(-(1/delta +1))-(gamma1+N)^(-(1/delta +1)))
    t1<- digamma(1/delta)-log(gamma1+N-1)
    t2<- digamma(1/delta)-log(gamma1+N)
    C2<- (1/P)*(gamma1^(1/delta))* (t1*((gamma1+N-1)^(-1/delta))-t2*((gamma1+N)^(-1/delta)))
    
    #### M step
    delta<- 1/nlm(func_delta, p=1/delta)$estimate# newtonRaphson(Derivative, delta)$root # 
    gamma1<- 1/(delta*mean(C1))
    p<-1- exp(-1/(gamma1*delta))
    Devianceold<-Deviancenew
    Deviancenew <- sum(log(like(N,gamma1 = gamma1,delta = delta)+tol))
    # Output
    k=k+1
    Outi[k]<-k; outd[k]<-delta; outp[1]<-p; outD<-Deviancenew
  }
  Output <- data.frame(Outi,outd,outp,outD)
  names(Output) <- c("iteration","delta","p","log-lik values")
  result <- list(par=c(delta,p), Deviance=Deviancenew, data=Output)
  return(result)
}
N<- rdpareto(50,delta = 1,p=0.5)
fit <- EM.dpareto(N,delta = 0.5, p=0.5)
fit$par

