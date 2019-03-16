

#### Estimation of BEG parameters
fit.BEG <- function(data,level=0.95) ## data has to be a vector (X,N)
{
  N<-data[,2]
  X<-data[,1]
  n<-nrow(data)
  qt<-(1-level)/2
  z<-abs(qnorm(qt))
  b<-mean(N)/mean(X)
  p<-1/mean(N)
  sigmabb<-1/(b*b*p)
  sigmapp<-1/((1-p)*p*p)
  J= solve(matrix(c(sigmabb,0,0,sigmapp),byrow = 2, ncol = 2))
  lowerb<-b-z*sqrt(J[1,1]/n)
  lowerp<-p-z*sqrt(J[2,2]/n)
  upperb<-b+z*sqrt(J[1,1]/n)
  upperp<-p+z*sqrt(J[2,2]/n)
  log.like<-log(p)+ N*log(b)-lgamma(N)+(N-1)*log((1-p)*X)-b*X
  Deviance<- -2*sum(log.like)
  Output<-data.frame(matrix(c(b,p)),matrix(c(lowerb,lowerp)),matrix(c(upperb,upperp)))
  colnames(Output)<-c("estimate",paste(level*100,"%", " lower bound", sep=""),
                      paste(level*100,"%", " upper bound", sep=""))
  row.names(Output)<- c("beta","p")
  result <- list(Estimates=Output,Deviance=Deviance,Inverse.Fisher.Matrix=J)
  return(result)
}  


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
    p=1/mean(N)
    ll<- N*par[1]*log(par[2]+1e-8)-lgamma(par[1]*N)+(par[1]*N-1)*log(X)-par[2]*X
    return(-sum(ll))
  }
  ## Estimate alpha
  beta<-mean(N)/mean(X)
  alpha<- 1
   para<-optim(par = c(alpha, beta), fn=log.lik)
   a<-para$par[1] #estimate for alpha
   b<-para$par[2] ## estimate beta
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


#### Estimation of BLG parameters

fit.BLG <- function(data,alpha=NULL,beta=NULL,level=0.95, tol=1e-8) ## data has to be a vector (X,N)
{
  if (is.null(alpha)| is.null(beta)){
    message("Initial values for alpha or beta are rquired for optimization")
  }else {
  N<-data[,2]
  X<-data[,1]
  qt<-(1-level)/2
  z<-abs(qnorm(qt))
  ### Log likelihood function to optimze for  alpha and beta
  log.lik <- function(par) { #par[1]=alpha and par[2]=beta
    N1<-NULL
    j=1
    for (i in N){
      t<-seq(0,(i-1),1)
      N1[j]<-sum(log(1+par[1]*t))
      j=j+1
    }
    ll<- mean(N1) -mean((N + 1/par[1])*log(1+par[1]*par[2]*X)) + mean(N)*log(par[2]+tol)
    return(ll)
  }
  ## Estimate alpha
  para<-optim(par = c(alpha,beta),fn= log.lik)
  a<-para$par[1]
  b<-para$par[2]
  p<-1/mean(N) # estimate p
  ### Covariance matrix
  W1<- (1-p)*mean(N*N/((1+a*N)^2))/p +(2/(a^2))*(1+(1-p)*mean(1/(1+a*N))/p)
  W2<- (2/(a^2))*mean(N/(1+a*N))+(1/a)*mean(N*(N+1)/(1+a*(N+1)))
  sigmaaa<- W1-W2
  sigmabb<- (1/(b*b))*mean(N/(1+a+a*N))
  sigmaab<- -(1/b)*mean(N/((1+a*N)*(1+a+a*N)))
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
  upperp<-p-z*sqrt(J[3,3]/n)
  log_like<- lgamma(N +1/a)-lgamma(1/a)-lgamma(N)+log(p)+(N-1)*log((1-p)*X)+N*log(a*b)-(1/a +N)*log(1+a*b*X)
  Deviance<- -2*sum(log_like)
  Output<-data.frame(matrix(c(a,b,p)),matrix(c(lowera,lowerb,lowerp)),matrix(c(uppera,upperb,upperp)))
  colnames(Output)<-c("estimate",paste(level*100,"%", " lower bound", sep=""),
                      paste(level*100,"%", " upper bound", sep=""))
  row.names(Output)<- c("alpha","beta","p")
  result <- list(Estimates=Output,Deviance=Deviance, Inverse.Fisher.Matrix=J)
  return(result)
  }
}  

# Example
Data.df<-rbeg(2,100,0.45) 
fit1<-fit.BEG(Data.df,level = 0.95)
fit1$Estimates

Data.df<-rblg(100,alpha = 1,beta = 100,p=0.4)
fit2<-fit.BLG(data = Data.df,alpha = 2,beta = 10, level = 0.95)
fit2$Estimates

Data.df<-rbgg(10,100,2,0.45)
fit3<-fit.BGG(data = Data.df)
fit3$Estimates
