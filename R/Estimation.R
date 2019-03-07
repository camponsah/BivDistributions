

#### Estimation of BEG parameters
fit.BEG <- function(data,level=0.95) ## data has to be a vector (X,N)
{
  N<-data[,2]
  X<-data[,1]
  qt<-(1-level)/2
  z<-abs(qnorm(qt))
  b<-mean(X)/mean(N)
  p<-1/mean(N)
  lowerb<-b-z*sqrt(p*b^2/length(data[,1]))
  lowerp<-p-z*sqrt(p*p*(1-p)/length(data[,2]))
  upperb<-b+z*sqrt(p*b^2/length(data[,1]))
  upperp<-p+z*sqrt(p*p*(1-p)/length(data[,2]))
  log.like<-((p*b^N)/gamma(N))*((X*(1-p))^(N-1))*exp(-b*X)
  Deviance<- -2*sum(log(log.like))
  Output<-data.frame(matrix(c(b,p)),matrix(c(lowerb,lowerp)),matrix(c(upperb,upperp)))
  colnames(Output)<-c("estimate",paste(level*100,"%", " lower bound", sep=""),
                      paste(level*100,"%", " upper bound", sep=""))
  row.names(Output)<- c("beta","p")
  result <- list(Estimates=Output,Deviance=Deviance)
  return(result)
}  


#### Estimation of BEG parameters

fit.BGG <- function(data,alpha,level=0.95) ## data has to be a vector (X,N)
{
  N<-data[,2]
  X<-data[,1]
  qt<-(1-level)/2
  z<-abs(qnorm(qt))
  ### Log likelihood function
    log.lik <- function(par) { #par[1]=p and par[2]=delta
    #N<-data[,2]
    #X<-data[,1]
    b=par*mean(N)/mean(X)
    p=1/mean(N)
    ll1<- b^(N*par)*(1/gamma(N*par))*X^(N*par -1)
    ll2<- exp(-b*X)*p*(1-p)^(N-1)
    D<- sum(log(ll1*ll2))
    -D
  }
  ## Estimate alpha
  a<-as.numeric(nlm(log.lik, p=alpha)$estimate)
  b<-a*mean(X)/mean(N) ## estimate beta
  p<-1/mean(N) # estimate p
  ### Sum to infinity function
  sumToInfinity<- function(p,a){
    j=1
    error=1
    S1=0
    while(error>0.00001){ 
      S=S1+ p*j^2 *((1-p)^(j-1))* psigamma(j*a, deriv = 1)
      j=j+1
      error=S-S1
      S1=S
    }
    return(S1)
  }
  sigmaaa<- sumToInfinity(p, a)
    sigmabb<- a/(p*b^2)
  sigmaab<--1/(p*b)
  sigmapp<-1/(p*p*(1-p))
 J= solve(matrix(c(sigmabb,sigmaab,0,sigmaab,sigmaaa,0,0,0,sigmapp),byrow = 3, ncol = 3))
 J<-data.frame(J)
 colnames(J)<- c("alpha","beta","p")
 row.names(J)<- c("alpha","beta","p")
 lowera<-a-z*sqrt(J[2,2]/length(X))
  lowerb<-b-z*sqrt(J[1,1]/length(X))
  lowerp<-p-z*sqrt(J[3,3]/length(N))
  uppera<-a+z*sqrt(J[2,2]/length(X))
  upperb<-b+z*sqrt(J[1,1]/length(X))
  upperp<-p-z*sqrt(J[3,3]/length(N))
  log.like<-b^(N*a)*(1/gamma(N*a))*(X^(N*a -1))*exp(-b*X)*p*(1-p)^(N-1)
  Deviance<- -2*sum(log(log.like))
  Output<-data.frame(matrix(c(a,b,p)),matrix(c(lowera,lowerb,lowerp)),matrix(c(uppera,upperb,upperp)))
  colnames(Output)<-c("estimate",paste(level*100,"%", " lower bound", sep=""),
                      paste(level*100,"%", " upper bound", sep=""))
  row.names(Output)<- c("alpha","beta","p")
  result <- list(Estimates=Output,Deviance=Deviance, Inverse.Fisher.Matrix=J)
  return(result)
}  


# Example
#fit<-fit.BEG(Data.df,level = 0.9)
#fit$Estimates
#fit<-fit.BGG(data = Data.df,alpha = 1, level = 0.99)
#fit$Estimates
