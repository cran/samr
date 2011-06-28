\name{SAMseq}
\alias{SAMseq}
\title{Significance analysis of sequencing data - simple user interface}
\description{
Correlates a large number of features (eg. genes) with an outcome
variable, such as a group indicator, quantitative variable or survival time.
This is a simple user interface for the samr function applied to sequencing data.
For array data applications, see the function SAM.

}
\usage{
SAMseq(x, y, censoring.status = NULL, 
resp.type = c("Quantitative", "Two class unpaired", 
"Survival", "Multiclass", "Two class paired"), 
geneid = NULL, genenames = NULL, nperms = 100, 
random.seed = NULL, nresamp = 20, fdr.output = 0.20)
}
\arguments{
\item{x}{Feature matrix: p (number of features) by n (number of samples),
one observation per column (missing values allowed)}
\item{y}{n-vector of outcome measurements}
\item{censoring.status}{n-vector of censoring censoring.status (1=died or event occurred, 0=survived, or event was censored), needed for a censored survival outcome}
\item{resp.type}{Problem type:
"Quantitative" for a continuous parameter;
"Two class unpaired" for two classes with unpaired observations; 
"Survival" for censored survival outcome; 
"Multiclass": more than 2 groups;  
"Two class paired" for two classes with paired observations.}
\item{geneid}{Optional character vector of geneids for output.}
\item{genenames}{Optional character vector of genenames for output.}
\item{nperms}{Number of permutations used to estimate false discovery rates}
\item{random.seed}{Optional initial seed for random number generator (integer)}
\item{nresamp}{Number of resamples used to construct test statistic. Default 20.}
\item{fdr.output}{(Approximate) False Discovery Rate cutoff for output in significant genes table}

}

\details{This is a simple, user-friendly interface to the samr package used on sequencing data.
It automatically disables arguments/features that do not apply to sequencing data.
It calls samr, samr.compute.delta.table and samr.compute.siggenes.table.
samr detects differential expression for microarray data,
and sequencing data,
and other data with a large number of features. samr is the R package
that is called by the "official" SAM Excel Addin.
The format of the response vector y and the calling sequence
is illustrated in the examples below. A more complete description
is given in the SAM manual
at http://www-stat.stanford.edu/~tibs/SAM}



\value{
A list with components
\item{samr.obj}{Output of samr. See documentation for samr for details }
\item{siggenes.table}{Table of significant genes, output of  samr.compute.siggenes.table.
This has components: genes.up---matrix of significant genes having positive correlation with the outcome and
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
set.seed(100)
mu <- matrix(100, 1000, 20)
mu[1:100, 11:20] <- 200
mu <- scale(mu, center=FALSE, scale=runif(20, 0.5, 1.5))
x <- matrix(rpois(length(mu), mu), 1000, 20)
y <- c(rep(1, 10), rep(2, 10))

samfit <- SAMseq(x, y, resp.type = "Two class unpaired") 

# examine significant gene list
print(samfit)

# plot results
plot(samfit)

######### two class paired comparison
set.seed(100)
mu <- matrix(100, 1000, 20)
mu[1:100, 11:20] <- 200
mu <- scale(mu, center=FALSE, scale=runif(20, 0.5, 1.5))
x <- matrix(rpois(length(mu), mu), 1000, 20)
y <- c(-(1:10), 1:10)

samfit <- SAMseq(x, y, resp.type = "Two class paired") 

# examine significant gene list
print(samfit)

# plot results
plot(samfit)

######### Multiclass comparison
set.seed(100)
mu <- matrix(100, 1000, 20)
mu[1:20, 1:5] <- 120
mu[21:50, 6:10] <- 80
mu[51:70, 11:15] <- 150
mu[71:100, 16:20] <- 60
mu <- scale(mu, center=FALSE, scale=runif(20, 0.5, 1.5))
x <- matrix(rpois(length(mu), mu), 1000, 20)
y <- c(rep(1:4, rep(5, 4)))

samfit <- SAMseq(x, y, resp.type = "Multiclass") 

# examine significant gene list
print(samfit)

# plot results
plot(samfit)

######### Quantitative comparison
set.seed(100)
mu <- matrix(100, 1000, 20)
y <- runif(20, 1, 3)
mu[1 : 100, ] <- matrix(rep(100 * y, 100), ncol=20, byrow=TRUE)
mu <- scale(mu, center=FALSE, scale=runif(20, 0.5, 1.5))
x <- matrix(rpois(length(mu), mu), 1000, 20)
samfit <- SAMseq(x, y, resp.type = "Quantitative")

# examine significant gene list
print(samfit)

# plot results
plot(samfit)

######### Survival comparison
set.seed(100)
mu <- matrix(100, 1000, 20)
y <- runif(20, 1, 3)
mu[1 : 100, ] <- matrix(rep(100 * y, 100), ncol=20, byrow=TRUE)
mu <- scale(mu, center=FALSE, scale=runif(20, 0.5, 1.5))
x <- matrix(rpois(length(mu), mu), 1000, 20)
y <- y + runif(20, -0.5, 0.5)
censoring.status <- as.numeric(y < 2.3)
y[y >= 2.3] <- 2.3
samfit <- SAMseq(x, y, censoring.status = censoring.status, 
resp.type = "Survival")

# examine significant gene list
print(samfit)

# plot results
plot(samfit)
}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{survival}
\keyword{ts}
\keyword{nonparametric}

