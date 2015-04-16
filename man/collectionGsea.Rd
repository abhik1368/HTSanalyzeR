\name{collectionGsea}
\alias{collectionGsea}
\title{
Compute observed and permutation-based enrichment scores for a collection
(list) of gene sets
}
\description{
This function computes observed and permutation-based scores associated
with a gene set enrichment analysis for a collection of gene sets.
}
\usage{
collectionGsea(collectionOfGeneSets, geneList, exponent=1, nPermutations=
1000, minGeneSetSize=15, verbose=TRUE)
}
\arguments{
  \item{collectionOfGeneSets}{
a list of gene sets. Each gene set in the list is a character vector of
gene identifiers. 
}
  \item{geneList}{
a numeric or integer vector which has been named and ordered. It cannot
contain any duplicates nor NAs.
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
  \item{verbose}{
a single logical value specifying to display detailed messages (when
verbose=TRUE) or not (when verbose=FALSE)
}
}
\value{
\item{Observed.scores}{The observed scores for the given gene sets (a 
named vector)}
\item{Permutation.scores}{The scores for the permutation tests (one
column for each permutation and a row for each gene set)}
}
\references{
Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L.,
Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E.
S. & Mesirov, J. P. (2005) \emph{ Gene set enrichment analysis:
A knowledge-based approach for interpreting genome-wide expression
profiles.} Proc. Natl. Acad. Sci. USA 102, 15545-15550.
}
\author{
Camille Terfve, Xin Wang
}
\seealso{
\code{\link[HTSanalyzeR:FDRcollectionGsea]{FDRcollectionGsea}}
}
\examples{
##example 1
gl <- runif(100, min=0, max=5)
gl <- gl[order(gl, decreasing=TRUE)]
names(gl) <- as.character(sample(x=seq(from=1, to=100, by=1), size=100,
replace=FALSE))
gs1 <- sample(names(gl), size=20, replace=FALSE)
gs2 <- sample(names(gl), size=20, replace=FALSE)
gsc <- list(subset1=gs1, subset2=gs2)
GSCscores <- collectionGsea(collectionOfGeneSets=gsc, geneList=gl,
exponent=1, nPermutations=1000, minGeneSetSize=5)
GSCpvalues <- permutationPvalueCollectionGsea(permScores=
GSCscores$Permutation.scores, dataScores=GSCscores$Observed.scores)
##example 2
\dontrun{
library(org.Dm.eg.db)
library(KEGG.db)
##load phenotype vector (see the vignette for details about the
##preprocessing of this data set)
data("KcViab_Data4Enrich")
DM_KEGG <- KeggGeneSets(species="Dm")
GSCscores <- collectionGsea(collectionOfGeneSets=DM_KEGG, geneList=
KcViab_Data4Enrich, exponent=1, nPermutations=1000, minGeneSetSize=100)
GSCpvalues <- permutationPvalueCollectionGsea(permScores=
GSCscores$Permutation.scores, dataScores=GSCscores$Observed.scores)
}
}
