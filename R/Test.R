#' method=c(Gnedenk, Bryson, Kozubowski)
#### Estimation of BEG parameters
ep.test <- function(data,sigma1=1,method="Gnedenk",level=0.95) ## data has to be a vector (X,N)
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
   Q<-function(sigma1){
    ll<- 1+log(sigma1)+log(mean(log(1+X/sigma1))) +mean(log(1+X/sigma1))
    -ll
  }
  sig<-as.numeric(nlm(Q, p=sigma1)$estimate)
  if((n==1)|| (sig==Inf)){
   Sn<-mean(X) 
   Wn<-0
  }else {
  Sn<- sig*mean(log(1+X/sig))
  Wn<-Sn/sig
  }
  log.lik<- - n*(log(Sn)+(1+1/Wn)*mean(log(1+Wn*X/Sn)))
  #Deviance<- -2*sum(log.lik)
  L0<- (exp(1)*mean(X))^(-n)
  L1<-exp(log.lik)
  test<- -2*log(L0/L1)
  p_value<-pchisq(test,df=2, lower.tail = FALSE)
  }
  ###Output
  Output<-data.frame(test,p_value)
  colnames(Output)<-c("F statistic","p-value")
  result <- Output
  return(result)
}

library(EnvStats)
#X<-rpareto(100, location=1, shape = 20)
X<-rexp(3, rate = 10)
ep.test(X, method = "Gnedenk")
ep.test(X, method = "Kozubowski")
