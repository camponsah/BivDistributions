

#### Estimation of BEG parameters
fit.BEG <- function(data,level=0.95) ## data has to be a vector (X,N)
{
  N<-data[,2]
  X<-data[,1]
  alpha<-(1-level)/2
  z<-abs(qnorm(alpha))
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

# Example
fit<-fit.BEG(Data.df,level = 0.9)
