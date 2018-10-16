\name{samr.compute.siggenes.table}
\alias{samr.compute.siggenes.table}
\title{Compute significant genes table}
\description{
 Computes significant genes table, starting with samr object "samr.obj" and  delta.table "delta.table"
}
\usage{
samr.compute.siggenes.table(samr.obj, del, data, delta.table, 
min.foldchange=0, all.genes=FALSE, compute.localfdr=FALSE)
}
\arguments{
\item{samr.obj}{ Object returned from call to samr}
\item{del}{Value of delta to define cutoff rule}
\item{data}{Data object, same as that used in call to samr}
\item{delta.table}{ Object returned from call to
  samr.compute.delta.table}
\item{min.foldchange}{The minimum fold change desired; should be >1;
default is zero, meaning no fold change criterion is applied}
\item{all.genes}{Should all genes be listed? Default FALSE}
\item{compute.localfdr}{Should the local fdrs be computed (this can take some time)? Default FALSE}
}
\value{
 return(list(genes.up=res.up,
genes.lo=res.lo, 
color.ind.for.multi=color.ind.for.multi, 
ngenes.up=ngenes.up, 
ngenes.lo=ngenes.lo))

\item{genes.up}{Matrix of significant genes having posative correlation with the outcome.
For survival data, genes.up are those genes having positive correlation with risk-
 that is, increased expression corresponds to higher risk (shorter survival).}
\item{genes.lo}{Matrix of significant genes having negative correlation with the outcome.
For survival data,genes. lo are those whose increased expression corresponds to lower risk (longer survival).}
\item{color.ind.for.multi}{For multiclass response: a matrix with entries +1
if the class mean is larger  than the overall mean at the 95% significance
levels, -1 if less, and zero otehrwise. This is useful in 
determining which class or classes causes a feature to be significant}
\item{ngenes.up}{Number of significant genes with  positive correlation}
\item{ngenes.lo}{Number of significant genes with  negative correlation}
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


samr.obj<-samr(data,  resp.type="Two class unpaired", nperms=100)

delta.table<-samr.compute.delta.table(samr.obj)
del<- 0.3
siggenes.table<- samr.compute.siggenes.table(samr.obj, del, data, delta.table)


}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{survival}
\keyword{ts}
\keyword{nonparametric}

