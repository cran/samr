\name{samr.missrate}
\alias{samr.missrate}
\title{Estimate the  miss rate table for a SAM analysis}
\description{
Estimates the  miss rate table, showing the local false negative rate, for a SAM
analysis
}
\usage{
samr.missrate(samr.obj, del, delta.table, quant=NULL)
}
\arguments{
\item{samr.obj}{ Object returned from call to samr}
\item{del}{Value of delta to define cutoff rule}
\item{delta.table}{ Object returned from call to samr.compute.delta.table}
\item{quant}{Optional vector of quantiles to used in the miss rate calculation}
}
\references{Tusher, V.,  Tibshirani, R.  and Chu, G. (2001):
Significance analysis of microarrays applied to the ionizing radiation response"  PNAS 2001 98: 5116-5121, (Apr 24).
http://www-stat.stanford.edu/~tibs/sam

Taylor, J., Tibshirani, R. and Efron. B. (2005). 
The ``Miss rate'' for the analysis of gene expression data.  Biostatistics 2005 6(1):111-117. 
}

\author{Jun Li and Balasubrimanian Narasimhan and Robert Tibshirani}


\examples{
#generate some example data
set.seed(100)
x<-matrix(rnorm(1000*20),ncol=20)
dd<-sample(1:1000,size=100)

u<-matrix(2*rnorm(100),ncol=10,nrow=100)
x[dd,11:20]<-x[dd,11:20]+u

y<-c(rep(1,10),rep(2,10))

data=list(x=x,y=y, geneid=as.character(1:nrow(x)),
genenames=paste("g",as.character(1:nrow(x)),sep=""), logged2=TRUE)


samr.obj<-samr(data,  resp.type="Two class unpaired", nperms=100)

delta.table<-samr.compute.delta.table(samr.obj)
del<- 0.3
siggenes.table<- samr.compute.siggenes.table(samr.obj, del, data, delta.table)


samr.missrate(samr.obj, del, delta.table)
}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{survival}
\keyword{ts}
\keyword{nonparametric}

