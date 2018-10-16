\name{samr.pvalues.from.perms}
\alias{samr.pvalues.from.perms}
\title{Report estimated p-values for each gene, from a  SAM analysis}
\description{Report estimated p-values for each gene, from set of
permutations in a SAM analysis
}
\usage{
samr.pvalues.from.perms(tt, ttstar)
}
\arguments{
\item{tt}{Vector of gene scores, returned by samr in component tt}
\item{ttstar}{Matrix of gene scores (p by nperm) from nperm permutations.
Returned by samr in component ttstar}
}
\references{Taylor, J. and Tibshirani, R.   (2005):
A tail strength measure for assessing the overall significance in a dataset.
Submitted.
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

pv=samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar)
}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{survival}
\keyword{nonparametric}

