\name{samr.compute.delta.table}
\alias{samr.compute.delta.table}
\title{Compute delta table for SAM analysis}
\description{
Computes tables of thresholds, cutpoints and corresponding
False Discovery rates for SAM (Significance analysis of microarrays)
analysis
}
\usage{

samr.compute.delta.table(samr.obj, min.foldchange=0, dels=NULL, nvals=50)
}
\arguments{
  \item{samr.obj}{ Object returned from call to samr}
 \item{min.foldchange}{The minimum fold change desired; should be >1;
default is zero, meaning no fold change criterion is applied}
\item{dels}{vector of delta values used. 
Delta is the vertical distance from the 45 degree line to the upper and lower
parallel lines that define the SAM threshold rule.
By default, for array data, 50 values
are chosen in the relevant operating change for delta.
For sequencing data, the maximum number of effective delta values
are chosen automatically according to the data.
}

\item{nvals}{Number of delta values used. 
For array data, the default value is 50. For sequencing data,
the value will be chosen automatically.}
}
\details{Returns a table of the FDR and upper and lower cutpoints for various values of
delta, for a SAM analysis.
}

\references{Tusher, V.,  Tibshirani, R.  and Chu, G. (2001): 
Significance analysis of microarrays applied to the ionizing radiation response"  PNAS 2001 98: 5116-5121, (Apr 24).
http://www-stat.stanford.edu/~tibs/sam}
\author{ Balasubrimanian Narasimhan and Robert Tibshirani}


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


samr.obj<-samr(data,  resp.type="Two class unpaired", nperms=50)

delta.table<- samr.compute.delta.table(samr.obj)

}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{survival}
\keyword{ts}
\keyword{nonparametric}
