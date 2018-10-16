\name{samr.estimate.depth}
\alias{samr.estimate.depth}
\title{estimate the sequencing depth}
\description{
Estimate the sequencing depth of each experiment for sequencing data.
}
\usage{
samr.estimate.depth(x)}
\arguments{
  \item{x}{the original count matrix. p by n matrix of features, one observation per column.}
}

\details{normalize the data matrix so that each number looks roughly like
Gaussian distributed and each experiment has the same sequencing depth.
To do this, we first use Anscombe transformation to stablize the variance
and makes each number look like Gaussian,
and then divide each experiment by the square root of the sequencing depth.}



\value{
\item{depth}{sequencing depth of each experiment. a vector of length n.}
}

 \references{Tusher, V.,  Tibshirani, R.  and Chu, G. (2001): 
Significance analysis of microarrays applied to the ionizing radiation response  PNAS 2001 98: 5116-5121, (Apr 24). 
http://www-stat.stanford.edu/~tibs/SAM}

\author{Jun Li and Balasubrimanian Narasimhan and Robert Tibshirani}


\examples{
set.seed(100)
mu <- matrix(100, 1000, 20)
mu[1:100, 11:20] <- 200
mu <- scale(mu, center=FALSE, scale=runif(20, 0.5, 1.5))
x <- matrix(rpois(length(mu), mu), 1000, 20)
y <- c(rep(1, 10), rep(2, 10))
data=list(x=x,y=y, geneid=as.character(1:nrow(x)),
genenames=paste("g",as.character(1:nrow(x)),sep=""))
depth <- samr.estimate.depth(data$x)
}

\keyword{univar}% at least one, from doc/KEYWORDS

