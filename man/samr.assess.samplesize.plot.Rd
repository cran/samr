\name{samr.assess.samplesize.plot}
\alias{samr.assess.samplesize.plot}
\title{Make a plot of the results from samr.assess.samplesize}
\description{
Plots of the results from samr.assess.samplesize
}
\usage{
samr.assess.samplesize.plot(samr.assess.samplesize.obj, logx=TRUE, 
call.win.metafile=FALSE)
}
\arguments{
\item{samr.assess.samplesize.obj}{ Object returned from call to samr.assess.samplesize}
\item{logx}{Should logs be used on the horizontal (\# of genes) axis? Default TRUE} 
\item{call.win.metafile}{Used by Excel interface}
}
\references{Tusher, V.,  Tibshirani, R.  and Chu, G. (2001):
Significance analysis of microarrays applied to the ionizing radiation response"  PNAS 2001 98: 5116-5121, (Apr 24).
http://www-stat.stanford.edu/~tibs/sam
}


\details{
Plots results: FDR (or 1-power) and FNR (or 1-type 1 error) from samr.assess.samplesize}



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

log2=function(x){log(x)/log(2)}

# run SAM first
samr.obj<-samr(data,  resp.type="Two class unpaired", nperms=100)

# assess current sample size (20), assuming 1.5fold difference on the log base 2 scale

samr.assess.samplesize.obj<- samr.assess.samplesize(samr.obj, data, log2(1.5))

samr.assess.samplesize.plot(samr.assess.samplesize.obj)
}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{survival}
\keyword{ts}
\keyword{nonparametric}

