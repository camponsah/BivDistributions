

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


# Example
Data.df<-rbeg(n=1000,beta = 10,p=0.45) 
fit<-fit.BEG(Data.df,level = 0.95)
fit
