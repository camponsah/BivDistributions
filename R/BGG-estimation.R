
#### Estimation of BGG parameters
fit.BGG <- function(data,level=0.95) ## data has to be a vector (X,N)
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
    ll<- par*mean(N)*log(B)-B*mean(X)+mean(par*N*log(X)-lgamma(par*N))
    return(-ll)
  }
  ## Estimate alpha
  alpha<-mean(X)/mean(N)
 a<- nlm(log.lik,p=alpha)$estimate #optim(par =alpha, fn=log.lik)
  #a<-para$par #estimate for alpha
  b<-a*mean(N)/mean(X) ## estimate beta
  p<-1/mean(N) # estimate p
  ### Sum to infinity function
  sumToInfinity<- function(p,a){
    j=1
    error=1
    S1=0
    while(error>0.00001){ 
      S=S1+ p*j*j *((1-p)^(j-1))* psigamma(j*a, deriv = 1)
      j=j+1
      error=abs(S-S1)
      S1=S
    }
    return(S1)
  }
  sigmaaa<- sumToInfinity(p, a)
  sigmabb<- a/(p*b*b)
  sigmaab<- -1/(p*b)
  sigmapp<-1/(p*p*(1-p))
  J= solve(matrix(c(sigmaaa,sigmaab,0,sigmaab,sigmabb,0,0,0,sigmapp),byrow = 3, ncol = 3))
  J<-data.frame(J)
  colnames(J)<- c("alpha","beta","p")
  row.names(J)<- c("alpha","beta","p")
  lowera<-a-z*sqrt(J[1,1]/n)
  lowerb<-b-z*sqrt(J[2,2]/n)
  lowerp<-p-z*sqrt(J[3,3]/n)
  uppera<-a+z*sqrt(J[1,1]/n)
  upperb<-b+z*sqrt(J[2,2]/n)
  upperp<-p+z*sqrt(J[3,3]/n)
  log.like<-N*a*log(b)-lgamma(a*N)+(a*N-1)*log(X)-b*X+log(p)+(N-1)*log(1-p)
  Deviance<- -2*sum(log.like)
  Output<-data.frame(matrix(c(a,b,p)),matrix(c(lowera,lowerb,lowerp)),matrix(c(uppera,upperb,upperp)))
  colnames(Output)<-c("estimate",paste(level*100,"%", " lower bound", sep=""),
                      paste(level*100,"%", " upper bound", sep=""))
  row.names(Output)<- c("alpha","beta","p")
  result <- list(Estimates=Output,Deviance=Deviance, Inverse.Fisher.Matrix=J)
  return(result)
}  

#Example
Data.df<-rbgg(1000,10,2,0.45)
fit<-fit.BGG(data = Data.df)
fit
