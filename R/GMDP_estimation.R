
#### Estimation of BGG parameters
gmdp_fit <- function(data,delta=1, p=0.5,level=0.95, method="EM") ## data has to be a vector (X,N)
{
  N<-data[,2]
  X<-data[,1]
  n<-nrow(data)
  qt<-(1-level)/2
  z<-abs(qnorm(qt))
  ### Log likelihood function
  log.lik <- function(par) { #par[1]=alpha, par[2]=beta
    #b=par*mean(N)/mean(X)
    B<-par*mean(N)/mean(X)
    ll<- par*mean(N)*log(B+1e-16)-B*mean(X)+mean(par*N*log(X)-lgamma(par*N))
    return(-ll)
  }
  log.lik.DP <- function(par) { #par[1]=p and par[2]=delta
    ll1<- (1/(1-par[1]*(N-1)*log(1-par[2])))^(1/par[1])
    ll2<- (1/(1-par[1]*N*log(1-par[2])))^(1/par[1])
    D<- sum(log(ll1-ll2))
    -D
  }
  ## Estimate alpha
  alpha<- mean(N*((mean(X))^2)/var(X))
  constant<- log(mean(N)/Mean(X)) + mean(N*log(X/N))/mean(N)
  if (constant<0){
    a<- nlm(log.lik,p=alpha)$estimate
  }
  else{
    a<- Inf
    message("MLE of alpha does not exist!")
  }
  b <- a*mean(N)/mean(X) ## estimate beta
  if( method == "EM"){
    source("fit_dp_EM.R")
    fit <- dpareto_em(N)
    delta <- fit$par$delta
    p <- fit$par$p
  }else{
    fit <-optim(par = c(delta, p), fn = log.lik.DP)
    delta <- fit$par[1]
    p <- fit$par[2]
  }

  Output<-data.frame(t(matrix(c(a,b,delta,p))))
  colnames(Output)<-c("alpha"," beta", "delta", "p")
  return(Output)
}

#Example
Data.df<-rgammamixdpareto(n=1000,alpha=1,beta=2,delta=0.45, p=0.8)
fit<-gmdp_fit(data = Data.df)
fit
