\name{aggregatePvals}
\alias{aggregatePvals}
\title{
Aggregate p-values from gene set over-representation tests.
}
\description{
This function takes as input a matrix of p-values for example obtained 
from a GSEA on multiple phenotypes, with a row for each gene set and a 
column for each phenotype and aggregates the p-values by row (i.e. one 
aggregated p-value for each gene set) according to Fisher or Stouffer's 
methods.
}
\usage{
aggregatePvals(pvalMatrix, method="fishers", pAdjustMethod="BH", 
order=TRUE)
}
\arguments{
  \item{pvalMatrix}{
a numeric matrix of p-values, with rows named according to the gene set 
(rows correspond to gene sets, and columns to multiple p-values to be 
aggregated for that gene set)
}
  \item{method}{
a single character value of "stouffers" or "fishers"
}
  \item{pAdjustMethod}{
a single character value specifying the p-value adjustment method to be 
used (see 'p.adjust' for details)
}
  \item{order}{
a single logical value: TRUE or FALSE. If it is TRUE, the results table 
will be ordered according to the aggregated p-values.
}
}
\value{
a matrix with a row for each gene set and two columns:
"Aggregated.p.value" and "Adjusted.aggregated.p.value"
}
\details{
The Fisher method combines the p-values into an aggregated chi-squared 
statistic equal to -2*sum(log(Pk)) were we have k=1,..,K  pvalues 
independently distributed as uniform on the unit interval under the null 
hypothesis.  The resulting p-values are calculated by comparing this 
chi-squared statistic to a chi-squared distribution with 2K degrees of 
freedom.
The Stouffer method computes a z-statistics assuming that the sum of the 
quantiles (from a standard normal distribution) corresponding to the 
p-values are distributed as N(0,K).
}
\author{
Jack Rose and Camille Terfve
}
\examples{
p1 <- runif(100, min=0, max=1)
p2 <- runif(100, min=0, max=1)
names(p1) <- as.character(sample(x=seq(from=1, to=100, by=1), size=100, 
replace=FALSE))
pmatrix <- cbind(p1, p2)
rownames(pmatrix) <- names(p1)
aggP <- aggregatePvals(pvalMatrix=pmatrix, method="stouffers")
}
