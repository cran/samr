\name{SAM}
\alias{SAM}
\title{Significance analysis of microarrays - simple user interface}
\description{
Correlates a large number of features (eg genes) with an outcome
variable, such as a group indicator, quantitative variable or survival time.
This is a simple user interface for the samr function applied to array data.
For sequencing data applications, see the function SAMseq.
}
\usage{
SAM(x,y=NULL,censoring.status=NULL,
resp.type=c("Quantitative","Two class unpaired","Survival","Multiclass",
"One class", "Two class paired","Two class unpaired timecourse",
"One class timecourse","Two class paired timecourse", "Pattern discovery"),
geneid = NULL,
genenames = NULL,
s0=NULL, 
s0.perc=NULL, 
nperms=100, 
center.arrays=FALSE, 
testStatistic=c("standard","wilcoxon"), 
time.summary.type=c("slope","signed.area"), 
regression.method=c("standard","ranks"), 
return.x=TRUE, 
knn.neighbors=10, 
random.seed=NULL,
logged2 = FALSE,
fdr.output = 0.20,
eigengene.number = 1)
}
\arguments{
  \item{x}{Feature matrix: p (number of features) by n (number of samples),
one observation per column (missing values allowed)}
\item{y}{n-vector of outcome measurements}
\item{censoring.status}{n-vector of censoring censoring.status (1= died or event occurred, 0=survived, or event was censored), needed for a censored survival outcome}
\item{resp.type}{Problem type:
"Quantitative" for a continuous parameter;
"Two class unpaired"; 
"Survival" for censored survival outcome; 
"Multiclass": more than 2 groups;  
"One class" for a single group; 
"Two class paired" for two classes with paired observations;
"Two class unpaired timecourse", 
"One class time course",
"Two class.paired timecourse",
"Pattern discovery"}
\item{geneid}{Optional character vector of geneids for output.}
\item{genenames}{Optional character vector of genenames for output.}
\item{s0}{Exchangeability factor  for denominator of test statistic; Default
is automatic choice. Only used for array data.}
\item{s0.perc}{Percentile of standard  deviation values to use for s0; default is
automatic choice; -1 means s0=0 (different from s0.perc=0, meaning
s0=zeroeth percentile of  standard  deviation values= min of sd values.
Only used for array data.}
\item{nperms}{Number of permutations used to estimate false discovery rates}
\item{center.arrays}{Should the data for each sample (array) be median centered
at the outset? Default =FALSE. Only used for array data.},
\item{testStatistic}{Test statistic to use in two class unpaired case.Either
"standard" (t-statistic) or ,"wilcoxon" (Two-sample wilcoxon or Mann-Whitney
test). Only used for array data.} 
\item{time.summary.type}{Summary measure for each time course: "slope", or "signed.area"). Only used for array data.} 
\item{regression.method}{Regression method for quantitative case:
"standard", (linear least squares) or "ranks" (linear least squares 
on ranked data). Only used for array data.}
\item{return.x}{Should the matrix of feature values be returned?
Only useful for time course data, where x contains summaries of the features
over time. Otherwise x is the same as the input data data\$x}
\item{knn.neighbors}{Number of nearest neighbors to use for imputation
of missing features values. Only used for array data.}
\item{random.seed}{Optional initial seed for random number generator (integer)}
\item{logged2}{Has the data been transformed by log (base 2)? This information is used only 
for computing fold changes }
\item{fdr.output}{(Approximate) False Discovery Rate cutoff for output in significant genes table}
\item{eigengene.number}{Eigengene to be used  (just for resp.type="Pattern discovery")}
}

\details{This is a simple, user-friendly interface to the samr package used on array data.
It calls samr, samr.compute.delta.table and samr.compute.siggenes.table.
samr detects differential expression for  microarray data,
and sequencing data,
and other data with a large number of features. samr is the R package
that is called by the "official" SAM Excel Addin.
The format of the response vector y and the calling sequence
is illustrated in the examples below. A more complete  description
is given in the SAM manual
at http://www-stat.stanford.edu/~tibs/SAM}



\value{
A list with components
\item{samr.obj}{Output of samr. See documentation for samr for details.}
\item{siggenes.table}{Table of significant genes, output of  samr.compute.siggenes.table.
This has components: genes.up--- matrix of significant genes having positive correlation with the outcome and
genes.lo---matrix of significant genes having negative correlation with the outcome.
For survival data, genes.up are those genes having positive correlation with risk-
 that is, increased expression corresponds to higher risk (shorter survival)
genes.lo are those whose increased expression corresponds to lower risk (longer survival).}
\item{delta.table}{Output of  samr.compute.delta.table.}
\item{del}{Value of delta (distance from 45 degree line in SAM plot) for 
 used for creating delta.table and siggenes.table. Changing the input value fdr.output
will change the resulting del.}
\item{call}{The calling sequence}
}

 \references{Tusher, V.,  Tibshirani, R.  and Chu, G. (2001): 
Significance analysis of microarrays applied to the ionizing radiation response  PNAS 2001 98: 5116-5121, (Apr 24). 
http://www-stat.stanford.edu/~tibs/SAM

Li, Jun and Tibshirani, R. (2011). Finding consistent patterns: a nonparametric
approach for identifying differential expression in
RNA-Seq data. To appear, Statistical Methods in Medical Research.
}
\author{Jun Li and Balasubrimanian Narasimhan and Robert Tibshirani}


\examples{

######### two class unpaired comparison
# y must take values 1,2

set.seed(100)
x<-matrix(rnorm(1000*20),ncol=20)
dd<-sample(1:1000,size=100)

u<-matrix(2*rnorm(100),ncol=10,nrow=100)
x[dd,11:20]<-x[dd,11:20]+u
y<-c(rep(1,10),rep(2,10))

samfit<-SAM(x,y,resp.type="Two class unpaired")

# examine significant gene list

print(samfit)

# plot results
plot(samfit)

########### two class paired

# y must take values  -1, 1, -2,2 etc, with (-k,k) being a pair

set.seed(100)
x<-matrix(rnorm(1000*20),ncol=20)
dd<-sample(1:1000,size=100)

u<-matrix(2*rnorm(100),ncol=10,nrow=100)
x[dd,11:20]<-x[dd,11:20]+u

y=c(-(1:10),1:10)



samfit<-SAM(x,y, resp.type="Two class paired",fdr.output=.25)




#############quantitative response


set.seed(30)
p=1000
x<-matrix(rnorm(p*20),ncol=20)
y<-rnorm(20)
x[1:20,y>0]=x[1:20,y>0]+4
a<-SAM(x,y,resp.type="Quantitative",nperms=50,fdr.output=.5)






###########survival data
# y is numeric; censoring.status=1 for failures, and 0 for censored

set.seed(84048)
x=matrix(rnorm(1000*50),ncol=50)
x[1:50,26:50]= x[1:50,26:50]+2
x[51:100,26:50]= x[51:100,26:50]-2

y=abs(rnorm(50))
y[26:50]=y[26:50]+2
censoring.status <- sample(c(0,1),size=50,replace=TRUE)

a<-SAM(x,y,censoring.status=censoring.status,resp.type="Survival",
nperms=20)



################multi-class example
# y takes values 1,2,3,...k where k= number of classes

set.seed(84048)
x=matrix(rnorm(1000*10),ncol=10)

y=c(rep(1,3),rep(2,3),rep(3,4))
x[1:50,y==3]=x[1:50,y==3]+5

a <- SAM(x,y,resp.type="Multiclass",nperms=50)






##################### pattern discovery
# here there is no outcome y; the desired eigengene is indicated by 
# the argument eigengene.numbe in the data object

set.seed(32)
x=matrix(rnorm(1000*9),ncol=9)
mu=c(3,2,1,0,0,0,1,2,3)
b=3*runif(100)+.5
x[1:100,]=x[1:100,]+ b%*%t(mu)



d=list(x=x,eigengene.number=1,
geneid=as.character(1:nrow(x)),genenames=paste("gene", as.character(1:nrow(x))))


a <- SAM(x, resp.type="Pattern discovery", nperms=50)


#################### timecourse data

# elements of y are of the form  kTimet  where k is the class label and t
# is the time; in addition, the   suffixes Start or End indicate the first
# and last observation in a given time course
# the class label can be that for a two class unpaired, one class or
# two class paired problem

set.seed(8332)
y=paste(c(rep(1,15),rep(2,15)),"Time",rep(c(1,2,3,4,5,1.1,2.5, 3.7, 4.1,5.5),3),
sep="")
start=c(1,6,11,16,21,26)
for(i in start){
y[i]=paste(y[i],"Start",sep="")
}
for(i in  start+4){
y[i]=paste(y[i],"End",sep="")
}
x=matrix(rnorm(1000*30),ncol=30)
x[1:50,16:20]=x[1:50,16:20]+matrix(3*c(0,1,2,3,4),ncol=5,nrow=50,byrow=TRUE)
x[1:50,21:25]=x[1:50,21:25]+matrix(3*c(0,1,2,3,4),ncol=5,nrow=50,byrow=TRUE)
x[1:50,26:30]=x[1:50,26:30]+matrix(3*c(0,1,2,3,4),ncol=5,nrow=50,byrow=TRUE)

x[51:100,16:20]=x[51:100,16:20]-matrix(3*c(0,1,2,3,4),ncol=5,nrow=50,byrow=TRUE)
x[51:100,21:25]=x[51:100,21:25]-matrix(3*c(0,1,2,3,4),ncol=5,nrow=50,byrow=TRUE)
x[51:100,26:30]=x[51:100,26:30]-matrix(3*c(0,1,2,3,4),ncol=5,nrow=50,byrow=TRUE)

a<- SAM(x,y,  resp.type="Two class unpaired timecourse",
 nperms=100, time.summary.type="slope")


}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{survival}
\keyword{ts}
\keyword{nonparametric}

