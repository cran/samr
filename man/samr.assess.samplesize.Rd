\name{samr.assess.samplesize}
\alias{samr.assess.samplesize}
\title{Assess the sample size for a SAM analysis}
\description{
Estimate the false discovery rate, 
false negative rate, power and type I error for a SAM
analysis. Currently
implemented only for two class (unpaired or paired), one-sample and survival problems). 
}
\usage{
samr.assess.samplesize(samr.obj, data, dif, samplesize.factors=c(1,2,3,5), 
min.genes = 10, max.genes = nrow(data$x)/2)
}
\arguments{
\item{samr.obj}{ Object returned from call to samr}
\item{data}{Data list, same as that passed to samr.train} 
\item{dif}{Change in  gene expression between groups 1 and 2, for genes that
are differentially expressed. For log base 2 data, a value of 1 means
a 2-fold change. For One-sample problems, dif is the number of units away from
zero for differentially expressed genes. For survival data, dif is the numerator of the Cox score statistic (this info is provided in the output of samr).}
\item{samplesize.factors}{Integer vector of length 4, indicating the
sample sizes to be examined.
The values are factors that multiply the original sample size.
So the value 1 means a sample size of ncol(data$x), 2 means a sample size of ncol(data$x), etc.}
\item{min.genes}{Minimum number of genes that are assumed to truly changed  in the population}
\item{max.genes}{Maximum number of genes that are assumed to truly changed  in the population}
}
\references{Tusher, V.,  Tibshirani, R.  and Chu, G. (2001):
Significance analysis of microarrays applied to the ionizing radiation response"  PNAS 2001 98: 5116-5121, (Apr 24).
http://www-stat.stanford.edu/~tibs/sam

Taylor, J., Tibshirani, R. and Efron. B. (2005). 
The ``Miss rate'' for the analysis of gene expression data.  Biostatistics 2005 6(1):111-117. 

 A more complete  description
is given in the SAM manual
at http://www-stat.stanford.edu/~tibs/SAM
}

\details{Estimates  false discovery rate,
false negative rate, power and type I error for a SAM
analysis. The argument samplesize.factor allows the use to assess the
effect of  varying the 
 sample size (total number of samples).  A  detailed  description
of this calculation
is given in the SAM manual
at http://www-stat.stanford.edu/~tibs/SAM}


\value{
A list with components
\item{Results}{A matrix with columns: number of genes- both  the number differentially expressed genes in the population and number called significant;
cutpoint- the threshold used for the absolute  SAM score d; FDR, 1-power- the median
false discovery rate, also equal to the power for each gene; FDR-90perc-
 the
upper 90th percentile of the FDR; FNR, Type 1 error- the false negative rate, also equal to 
the type I error for each gene; FNR-90perc- the upper 90th percentile of the FNR }
\item{dif.call}{Change in  gene expression between groups 1 and 2, that was provided in the
call to samr.assess.samplesize}
\item{difm}{The average difference in SAM score d for the genes
differentially expressed vs unexpressed}
\item{samplesize.factor}{The  samplesize.factor that was passed to  samr.assess.samplesiz}
\item{n}{Number of samples in input data (i.e. ncol of x component in data)}
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

log2=function(x){log(x)/log(2)}

# run SAM first
samr.obj<-samr(data,  resp.type="Two class unpaired", nperms=100)

# assess current sample size (20), assuming 1.5fold difference on log base 2 scale

samr.assess.samplesize.obj<- samr.assess.samplesize(samr.obj, data, log2(1.5))

# assess the effect of doubling the sample size

samr.assess.samplesize.obj2<- samr.assess.samplesize(samr.obj, data, log2(1.5))

}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{survival}
\keyword{ts}
\keyword{nonparametric}

