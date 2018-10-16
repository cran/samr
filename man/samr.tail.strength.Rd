\name{samr.tail.strength}
\alias{samr.tail.strength}
\title{Estimate tail strength for a dataset, from a SAM analysis}
\description{Estimate tail strength for a dataset, from a SAM analysis
}
\usage{
samr.tail.strength(samr.obj)
}
\arguments{
\item{samr.obj}{Object returned by samr}
}
\value{
A list with components
\item{ts}{Estimated tail strength. A number less than or equal to 1.
Zero means all genes are null; 1 means all genes are differentially expressed.},
\item{se.ts}{Estimated standard error of tail strength.}
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

samr.tail.strength(samr.obj)

}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{survival}
\keyword{nonparametric}

