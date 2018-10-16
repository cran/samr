\name{samr.plot}
\alias{samr.plot}

\title{Make Q-Q plot for SAM analysis}
\description{
Makes the Q-Q plot for a SAM analysis}

\usage{
samr.plot(samr.obj, del, min.foldchange=0) 
}
\arguments{
  \item{samr.obj}{ Object returned from call to samr}
\item{del}{Value of delta to use. 
Delta is the  vertical  distance from the 45 degree line to the upper and lower
parallel lines that define the SAM threshold rule.}
 \item{min.foldchange}{The minimum fold change desired; should be >1;
default is zero, meaning no fold change criterion is applied}
}

\details{Creates the Q-Q plot fro a SAm analysis, marking features (genes)
that are significant, ie. lie outside a slab around teh 45 degreee line of
width delta. A gene must also pass the  min.foldchange criterion 
to be called significant, if this criterion is specified in the call.
}

\references{Tusher, V.,  Tibshirani, R.  and Chu, G. (2001): 
Significance analysis of microarrays applied to the ionizing radiation response"  PNAS 2001 98: 5116-5121, (Apr 24).
http://www-stat.stanford.edu/~tibs/sam}
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


samr.obj<-samr(data,  resp.type="Two class unpaired", nperms=50)


samr.plot(samr.obj, del=.3)
}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{survival}
\keyword{ts}
\keyword{nonparametric}
