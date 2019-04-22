

dpareto_em <- function(data, delta = 0.1, p = 0.5, maxiter = 500, tol = 1e-16){
  N<- data
  n <- length(N)
  gamma_1<- - 1/(delta*log(1 - p))
  eta<- 1/delta
  # log-likelihood for discrete Pareto
  log_like<- function(N, delta, p){
   W1<- (1 - delta*(N - 1)*log(1 - p))^( - 1/delta)
   W2<- (1 - delta*N*log(1 - p))^( - 1/delta)
    return(sum(log(W1 - W2)+tol))
  }
  # Function to optimize to eta
  func_eta<-function(eta){
    ll<- eta*log(eta + tol) - lgamma(eta)- eta - eta*log(mean(a) + tol) + eta*mean(c)
    return( -ll)
  }
  # log-liklihood calcultion
  Devianceold<- 0
  Deviancenew <- log_like(N, delta, p)
  # Intitialize vectors for storing data
  Outi<- NULL; outd<- NULL; outp<-NULL; outD<- NULL; k = 1
  Outi[1]<- 0; outd[1]<- delta; outp[1]<- p; outD[1]<- Deviancenew
  while((abs(Deviancenew - Devianceold) > tol) & (k <= maxiter)){
    ### E step
    const<- 1/(gamma(eta)*((gamma_1 + N - 1)^(-eta) - (gamma_1 + N)^(-eta)))
    a<- const * gamma(eta + 1) * ((gamma_1 + N - 1)^(-(eta + 1)) - (gamma_1 + N)^(-(eta + 1)))
    t1<- digamma(eta) - log(gamma_1 + N - 1)
    t2<- digamma(eta) - log(gamma_1 + N)
    c<-  const*gamma(eta) * (t1*(gamma_1 + N - 1)^(-eta) - t2*(gamma_1 + N)^(-eta))
    #### M step
    constant<-mean(c)-log(mean(a))
    if(constant<0){
    eta <-stats:: nlm(f=func_eta,p=eta, ndigit = 12)$estimate
    }
    else{
      eta <- Inf
      }
    delta<- 1/eta
    p<- 1 - exp( - mean(a))
    gamma_1<- - 1/(delta*log(1 - p))
    Devianceold<-Deviancenew
    Deviancenew <- log_like(N, delta, p)
    # Output
    k<- k + 1
    Outi[k]<- k; outd[k]<- delta; outp[k]<- p; outD[k]<- Deviancenew
  }
  #Output data
  Output <- data.frame(Outi,outd,outp,outD)
  names(Output) <- c("iteration","delta","p","log-lik values")
  result <- list(par=c(delta,p), Deviance=Deviancenew, data=Output)
  return(result)
}

#Example
N<- rdpareto(10,delta =0.9 ,p=0.5)

fit <- dpareto_em(N, maxiter = 1000)
fit$par

