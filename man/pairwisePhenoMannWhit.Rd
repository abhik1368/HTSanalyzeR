\name{pairwisePhenoMannWhit}
\alias{pairwisePhenoMannWhit}
\title{
Mann-Whitney U test for shift in location of genes from gene sets on a pair
of phenotypes
}
\description{
This function performs a Mann-Whitney U test for shift in location of genes
from gene sets, on a pair of phenotypes. It looks for gene sets that are
represented towards the two different ends of two ranked lists of genes,
i.e. whose phenotype distribution is located around two different ends
of the two phenotype vectors, rather than spread on the whole list in both
lists.
}
\usage{
pairwisePhenoMannWhit(gl1, gl2, gsc, minGeneSetSize=15, pAdjustMethod="BH")
}
\arguments{
  \item{gl1}{
a named numeric or integer vector where names are gene identifiers of the
same type as the ones in the gene set collection, and values are the
measurements on phenotype one corresponding to those genes. This vector MUST
be ordered (decreasing or increasing)
}
  \item{gl2}{
a named numeric or integer vector where names are gene identifiers of the
same type as the ones in the gene set collection, and values are the
measurements on phenotype two corresponding to those genes. This vector MUST
be ordered
}
  \item{gsc}{
a list of gene sets, each of which in the list is a character vector of
gene identifiers. 
}
  \item{minGeneSetSize}{
a single numeric or integer value specifying the minimum size required
for a gene set to be considered.
}
  \item{pAdjustMethod}{
a single character value specifying the p-value adjustment method to be
used (see 'p.adjust' for details)
}
}
\value{
a table with a row for each gene set, containing the p-value for the
Mann-Whitney U test, and the adjusted p-value. The table is ordered by the
p-value column.
}
\author{
Camille Terfve, Xin Wang
}
\seealso{
\code{\link[HTSanalyzeR:pairwiseGseaPlot]{pairwiseGseaPlot}}, \code{\link[HTSanalyzeR:pairwiseGsea]{pairwiseGsea}}
}
\examples{
gl1 <- runif(100, min=-5, max=5)
gl2 <- runif(100, min=-5, max=5)
names(gl1) <- as.character(sample(x=seq(from=1, to=100, by=1), size=100,
replace=FALSE))
names(gl2) <- names(gl1)
gs1 <- sample(names(gl1), size=20, replace=FALSE)
gs2 <- sample(names(gl1), size=20, replace=FALSE)
gsc <- list(subset1=gs1, subset2=gs2)
pwGSC <- pairwisePhenoMannWhit(gl1=gl1, gl2=gl2, gsc=gsc, minGeneSetSize=2)
}
