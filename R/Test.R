
#### Estimation of BEG parameters
ep.test <- function(data,sigma1=1,level=0.95) ## data has to be a vector (X,N)
{
  X<-data
  qt<-(1-level)
   Q<-function(sigma1){
    ll<- 1+log(sigma1)+log(mean(log(1+X/sigma1))) +mean(log(1+X/sigma1))
    -ll
  }
  sig<-as.numeric(nlm(Q, p=sigma1)$estimate)
  if(sig==Inf){
   Sn<-mean(X) 
   Wn<-0
  }else {
  Sn<- sig*mean(log(1+X/sig))
  Wn<-Sn/sig
  }
  log.lik<- - length(X)*(log(Sn)+(1+1/Wn)*mean(log(1+Wn*X/Sn)))
  Deviance<- -2*sum(log.lik)
  L0<- - length(X)*(1+log(mean(X)))
  L1<-log.lik
  test<- -2*(L0-L1)
  C<-qchisq(1-2*qt,df=1)
  p_value<-pchisq(test,df=2, lower.tail = FALSE)
  Output<-data.frame(test,p_value)
  colnames(Output)<-c("chi-square statistic","p-value")
  result <- list(estimates=Output,deviance=Deviance)
  return(result)
}

ep.test(Data.df$X)
