#' method=c(Gnedenk, Bryson, Kozubowski)
#### Estimation of BEG parameters
ep.test <- function(data,sigma1=1,method="Gnedenk",level=0.95,tol=1e-12) ## data has to be a vector (X,N)
{
  X<-data
  n=length(X)
  qt<-(1-level)
  ## 1
  if (method=="Gnedenk"){
    r<- n-2
    S<-(n-seq(1,n)+1)*c(0,diff(sort(X)))
    test<-mean(S[1:r])/mean(S[(r+1):n])
    p_value<-pf(test, df1=2*r, df2=2*(n-r), lower.tail = FALSE) 
    
  }
  else if (method=="Bryson")
  {
    a<- mean(X)*max(X)
    b<-(n-1)*(exp(mean(log(X+(max(X)/(n-1))))))^2
    test<-a/b
    p_value<-pf(test, df1=2*(length(X)-1), df2=2, lower.tail = FALSE) 
  }
  ##2
  else if (method=="Kozubowski"){
   log.like<-function(par){ #par[1]=omega , par[2]=s
    ll<- 1+log(par+ tol ) + log(mean(log(1+X/par))+tol) + mean(log(1+X/par))
    return(-ll)
   }
   sig<- nlm(p=mean(X),f=log.like)$estimate
  ######################################
  ######################################
  Sn<- sig*mean(log(1+X/sig))
  Wn<-mean(log(1+X/sig))
  L0<- -n*(1+log(mean(X)))
  L1<- -n*(log(Sn+tol) +(1+(1/Wn))*mean(log(1+Wn*X/Sn)))
  test<- -2*(L0+L1)
  p_value<-0.5*pchisq(test,df=1, lower.tail = FALSE)
  }
  ###Output
  Output<-data.frame(test,p_value)
  colnames(Output)<-c("Test statistic","p-value")
  result <- Output
  return(result)
}

library(EnvStats)
#X<-rpareto(100, location=1, shape = 20)
X<-rexp(300, rate = 1)
ep.test(X, method = "Gnedenk")
ep.test(X, method = "Kozubowski")
