\name{pairwiseGsea}
\alias{pairwiseGsea}
\title{
GSEA on a pair of phenotypes
}
\description{
This function performs pairwise GSEA: it looks for gene sets that are
specifically overrepresented towards the two different ends of two ranked
lists of genes, in a gene set collection.
}
\usage{
pairwiseGsea(gl1, gl2, gsc, exponent=1, nPermutations=1000,
minGeneSetSize=15, pAdjustMethod="BH")
}
\arguments{
  \item{gl1}{
a named numeric or integer vector where names are gene identifiers of
the same type as the ones in the gene set collection, and values are
the measurement on phenotype one corresponding to those genes this vector
MUST be ordered (decreasing or increasing)
}
  \item{gl2}{
a named numeric or integer vector where names are gene identifiers of
the same type as the ones in the gene set collection, and values are the
measurement on phenotype two corresponding to those genes this vector MUST
be ordered
}
  \item{gsc}{
a list of gene sets, each of which in the list is a character vector. 
}
  \item{exponent}{
a single numeric or integer value (set as 1 by default) specifying the
exponent of the GSEA method. 
}
  \item{nPermutations}{
a single numeric or integer value specifying the number of permutation
tests for each gene set
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
a table with a row for each gene set, containing the p-values for the
GSEA, and the observed scores for each of the phenotypes independently.
The table is ordered by the p-value column.
}
\details{
phenotypes one and two must be measured on the same genes, i.e. the two
vectors gl1 and gl2 must have the same length and their names must match,
but the two vectors must be ordered separately, i.e.one phenotype vector
is ordered based on the values of that phenotype only
}
\references{
Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L.,
Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander,
E. S. & Mesirov, J. P. (2005) \emph{ Gene set enrichment analysis: A
knowledge-based approach for interpreting genome-wide expression profiles.}
Proc. Natl. Acad. Sci. USA 102, 15545-15550.
}
\author{
Camille Terfve, Xin Wang
}
\seealso{
\code{\link[HTSanalyzeR:pairwiseGseaPlot]{pairwiseGseaPlot}}, \code{\link[HTSanalyzeR:pairwisePhenoMannWith]{pairwisePhenoMannWith}}
}
\examples{
gl1 <- runif(100, min=-5, max=5)
gl1 <- gl1[order(gl1, decreasing=TRUE)]
gl2 <- runif(100, min=-5, max=5)
gl2 <- gl2[order(gl2, decreasing=TRUE)]
names(gl1) <- as.character(sample(x=seq(from=1, to=100, by=1), size=100,
replace=FALSE))
names(gl2) <- names(gl1)
gs1 <- sample(names(gl1), size=20, replace=FALSE)
gs2 <- sample(names(gl1), size=20, replace=FALSE)
gsc <- list(subset1=gs1, subset2=gs2)
pwGSEAscore <- pairwiseGsea(gl1=gl1, gl2=gl2, gsc=gsc)
}
