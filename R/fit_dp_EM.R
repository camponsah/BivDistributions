#' 
#' Algorithmn
dpareto_em <- function(data,delta, p, maxiter=100, tol=1e-12)
{
  require(pracma)
  # initialization
  N<-data
  n <- length(N)
  gamma1<- -1/(delta*log(1-p))
  eta<- 1/delta
  ##Likelihood function of Discrete Pareto
  log_like<-function(N,gamma1, eta){
    LL<- (1/(1+((N-1)/gamma1)))^(eta)-(1/(1+(N/gamma1)))^(eta)
    return(sum(log(LL)))
  }
  ### function to solve for eta
  func_eta<-function(eta){
    ll<- eta*log(eta+tol)-lgamma(eta)-eta -eta*log(mean(a)+tol) +eta*mean(c)
      #log(eta+tol)-digamma(eta)-log(a+tol)+c
    return(-ll) 
  }
  Devianceold<-0
  Deviancenew <- log_like(N,gamma1,eta)
  Outi<-NULL;outd<-NULL;outp<-NULL;outD<-NULL; k=1 
  Outi[1]<-0; outd[1]<-delta; outp[1]<-p; outD[1]<-Deviancenew
  ##
  while((abs(Deviancenew-Devianceold)>tol) & (k <= maxiter)){ 
    ### E step
    const<- 1/(gamma(eta)*((gamma1+N-1)^(-eta) - (gamma1+N)^(-eta)))
    a<- const * gamma(eta+1) * ((gamma1+N-1)^(-(eta+1)) - (gamma1+N)^(-(eta+1)))
    t1<- digamma(eta)-log(gamma1+N-1)
    t2<- digamma(eta)-log(gamma1+N)
    c<-  const*gamma(eta) * (t1*(gamma1+N-1)^(-eta) -t2*(gamma1+N)^(-eta))
    #### M step
    eta<- stats:: nlm(f=func_eta,p=eta)$estimate
    delta<- 1/eta
    p<-exp(-mean(a))
    gamma1<- - 1/(delta*log(1-p))
    Devianceold<-Deviancenew
    Deviancenew <- log_like(N,gamma1,eta)
    # Output
    k=k+1
    Outi[k]<-k; outd[k]<-delta; outp[k]<-p; outD[k]<-Deviancenew
  }
  Output <- data.frame(Outi,outd,outp,outD)
  names(Output) <- c("iteration","delta","p","log-lik values")
  result <- list(par=c(delta,p), Deviance=Deviancenew, data=Output)
  return(result)
}
N<- rdpareto(100,delta = 1,p=0.5)
fit <- dpareto_em(N,delta = 2, p=0.5, maxiter = 1000)
fit$par

