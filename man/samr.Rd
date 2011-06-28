\name{samr}
\alias{samr}
\title{Significance analysis of microarrays}
\description{
Correlates a large number of features (eg genes) with an outcome
variable, such as a group indicator, quantitative variable or  survival time.
NOTE: for most users, the interface function SAM--- which calls samr-- will be more 
convenient for array data, and the interface function SAMseq-- which also calls samr--
will be more convenient for sequencing data.
}
\usage{
samr(data, resp.type=c("Quantitative","Two class unpaired",
"Survival","Multiclass", "One class", "Two class paired",
"Two class unpaired timecourse", "One class timecourse",
"Two class paired timecourse", "Pattern discovery"),
assay.type=c("array","seq"), s0=NULL, s0.perc=NULL, nperms=100, 
center.arrays=FALSE, testStatistic=c("standard","wilcoxon"), 
time.summary.type=c("slope","signed.area"), 
regression.method=c("standard","ranks"), return.x=FALSE, 
knn.neighbors=10, random.seed=NULL, nresamp=20,nresamp.perm=NULL, 
xl.mode=c("regular","firsttime","next20","lasttime"), 
xl.time=NULL,  xl.prevfit=NULL)
}
\arguments{
  \item{data}{Data object with components x- p by n matrix of features,
one observation per column (missing values allowed); y- n-vector of outcome measurements;
censoring.status- n-vector of censoring censoring.status (1= died or event occurred, 0=survived, or event was censored), needed
for a censored survival outcome}
  \item{resp.type}{Problem type:
"Quantitative" for a continuous parameter (Available for both array and sequencing data);
"Two class unpaired" (for both array and sequencing data); 
"Survival" for censored survival outcome (for both array and sequencing data); 
"Multiclass": more than 2 groups (for both array and sequencing data);  
"One class" for a single group (only for array data); 
"Two class paired" for two classes with paired observations (for both array and sequencing data);
"Two class unpaired timecourse" (only for array data), 
"One class time course" (only for array data),
"Two class.paired timecourse" (only for array data),
or "Pattern discovery" (only for array data)}
\item{assay.type}{Assay type: "array" for microarray data, "seq" for counts from sequencing}
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
\item{nresamp}{For assay.type="seq", number of resamples used to construct test statistic. Default 20. Only used for sequencing data.}
\item{nresamp.perm}{For assay.type="seq", number of resamples used to construct test statistic for permutations. Default is equal to nresamp and it must be at most nresamp. Only used for sequencing data.}
\item{xl.mode}{Used by Excel interface}
\item{xl.time}{Used by Excel interface}
\item{xl.prevfit}{Used by Excel interface}
}

\details{Carries out a SAM analysis. Applicable to microarray data,
sequencing data,
and other data with a large number of features. This is the R package
that is called by the "official" SAM Excel package v2.0.
The format of the response vector y and the calling sequence
is illustrated in the examples below. A more complete  description
is given in the SAM manual
at http://www-stat.stanford.edu/~tibs/SAM}



\value{
A list with components
\item{n}{Number of observations}
\item{x}{Data matrix p by n (p=\# genes or features). Equal to the  matrix data\$x  in the original call to samr except
for (1) time course analysis, where is contains the summarized data
or (2) quantitative outcome with rank regression, where it contains
the data transformed to ranks. Hence it is null except for in time course analysis. }
\item{y}{Vector of n outcome values. equal the values data\$y in the original call to samr, except
for (1) time course analysis, where is contains the summarized y 
or (2) quantitative outcome with rank regression, where it contains
the y values transformed to ranks}
\item{argy}{The  values data\$y in the original call to samr} 
\item{censoring.status}{Censoring status indicators if applicable}
\item{testStatistic}{Test Statistic used},
\item{nperms}{Number of permutations requested}
\item{nperms.act}{Number of permutations actually used. Will be <
nperms when \# of possible permutations  <= nperms (in which case
all permutations are done)}
\item{tt}{tt=numer/sd, the  vector of p test statistics for original data}
\item{numer}{Numerators for tt}
\item{sd}{Denominators for tt. Equal to standard deviation for feature plus s0}
\item{s0}{Computed  exchangeability factor}
\item{s0.perc}{Computed percentile of standard deviation values.
s0= s0.perc percentile of the gene standard deviations}
\item{eva}{p-vector of expected values for tt under permutation sampling}
\item{perms}{nperms.act by n matrix  of permutations used. Each row is 
a permutation of 1,2...n}
\item{permsy}{nperms.act by n matrix  of permutations used. Each row is
a permutation of y1,y2,...yn. Only one of perms or permys is non-Null, depending on resp.type} 
\item{all.perms.flag}{Were all possible permutations used?}
\item{ttstar}{p by nperms.aca matrix t of test statistics from permuted data. Each column if sorted in descending order}
\item{ttstar0}{p by nperms.act matrix  of test statistics from permuted data. Columns are in order of data}
\item{eigengene.number}{The number of the eigengene (eg 1,2,..) that was requested for Pattern discovery}
\item{eigengene}{Computed eigengene}
\item{pi0}{Estimated proportion of non-null features (genes)}
\item{foldchange}{p-vector of foldchanges for original data}
\item{foldchange.star}{p by nperms.act matrix estimated foldchanges from permuted data}
\item{sdstar.keep}{n by nperms.act matrix of standard deviations
from each permutation}
\item{censoring.status.star.keep}{n by nperms.act matrix of 
censoring.status indicators from each permutation}
\item{resp.type}{The response type used. Same as resp.type.arg, except for
time course data, where time data is summarized and then treated as non-time course. Eg if resp.type.arg="oneclass.timecourse" then  resp.type="oneclass"}
\item{resp.type.arg}{The response type requested in the call to samr}
\item{stand.contrasts}{For multiclass data, p by nclass matrix of  standardized differences
between the class mean and the overall mean}
\item{stand.contrasts.star}{For multiclass data, p by nclass by nperms.act array of standardized contrasts for permuted datasets}
\item{stand.contrasts.95}{For multiclass data, 2.5% and 97.5% percentiles
of standardized contrasts. Useful for determining which class contrast
for significant genes, are large}
\item{depth}{For array.type="seq", estimated sequencing depth for each sample.}
\item{call}{calling sequence}
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

data=list(x=x,y=y, geneid=as.character(1:nrow(x)),
genenames=paste("g",as.character(1:nrow(x)),sep=""), logged2=TRUE)


samr.obj<-samr(data,  resp.type="Two class unpaired", nperms=100)

delta=.4
samr.plot(samr.obj,delta)

delta.table <- samr.compute.delta.table(samr.obj)

siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, data, delta.table)

# sequence data

set.seed(3)
x<-abs(100*matrix(rnorm(1000*20),ncol=20))
x=trunc(x)
y<- c(rep(1,10),rep(2,10))
x[1:50,y==2]=x[1:50,y==2]+50
data=list(x=x,y=y, geneid=as.character(1:nrow(x)),
genenames=paste("g",as.character(1:nrow(x)),sep=""))

samr.obj<-samr(data,  resp.type="Two class unpaired",assay.type="seq",  nperms=100)

delta=5
samr.plot(samr.obj,delta)

delta.table <- samr.compute.delta.table(samr.obj)

siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, data, delta.table)





########### two class paired

# y must take values  -1, 1, -2,2 etc, with (-k,k) being a pair

set.seed(100)
x<-matrix(rnorm(1000*20),ncol=20)
dd<-sample(1:1000,size=100)

u<-matrix(2*rnorm(100),ncol=10,nrow=100)
x[dd,11:20]<-x[dd,11:20]+u
y=c(-(1:10),1:10)


d=list(x=x,y=y, geneid=as.character(1:nrow(x)),
genenames=paste("g",as.character(1:nrow(x)),sep=""), logged2=TRUE)


samr.obj<-samr(d,  resp.type="Two class paired", nperms=100)




#############quantitative response

# y must take numeric values

set.seed(84048)
x=matrix(rnorm(1000*9),ncol=9)

mu=c(3,2,1,0,0,0,1,2,3)
b=runif(100)+.5
x[1:100,]=x[1:100,]+ b%*%t(mu)

y=mu

d=list(x=x,y=y, 
geneid=as.character(1:nrow(x)),genenames=paste("gene", as.character(1:nrow(x))))

samr.obj =samr(d,  resp.type="Quantitative", nperms=50)



########### oneclass
# y is a vector of ones

set.seed(100)
x<-matrix(rnorm(1000*20),ncol=20)
dd<-sample(1:1000,size=100)

u<-matrix(2*rnorm(100),ncol=10,nrow=100)
x[dd,11:20]<-x[dd,11:20]+u

y<-c(rep(1,20))

data=list(x=x,y=y, geneid=as.character(1:nrow(x)),
genenames=paste("g",as.character(1:nrow(x)),sep=""), logged2=TRUE)

samr.obj<-samr(data,  resp.type="One class", nperms=100)



###########survival data
# y is numeric; censoring.status=1 for failures, and 0 for censored

set.seed(84048)
x=matrix(rnorm(1000*50),ncol=50)
x[1:50,26:50]= x[1:50,26:50]+2
x[51:100,26:50]= x[51:100,26:50]-2

y=abs(rnorm(50))
y[26:50]=y[26:50]+2
censoring.status=sample(c(0,1),size=50,replace=TRUE)
d=list(x=x,y=y,censoring.status=censoring.status,
geneid=as.character(1:1000),genenames=paste("gene", as.character(1:1000)))

samr.obj=samr(d,  resp.type="Survival", nperms=20)


################multi-class example
# y takes values 1,2,3,...k where k= number of classes

set.seed(84048)
x=matrix(rnorm(1000*10),ncol=10)
x[1:50,6:10]= x[1:50,6:10]+2
x[51:100,6:10]= x[51:100,6:10]-2

y=c(rep(1,3),rep(2,3),rep(3,4))
d=list(x=x,y=y,geneid=as.character(1:1000),
genenames=paste("gene", as.character(1:1000))) 

samr.obj <- samr(d,  resp.type="Multiclass")



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

data=list(x=x,y=y, geneid=as.character(1:nrow(x)),
genenames=paste("g",as.character(1:nrow(x)),sep=""), logged2=TRUE)

samr.obj<- samr(data,  resp.type="Two class unpaired timecourse",
 nperms=100, time.summary.type="slope")


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


samr.obj=samr(d,  resp.type="Pattern discovery", nperms=50)


}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{survival}
\keyword{ts}
\keyword{nonparametric}

