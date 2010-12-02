\name{permutationPvalueCollectionGsea}
\alias{permutationPvalueCollectionGsea}
\title{
Compute the GSEA p-values for a collection (list) of gene sets
}
\description{
Compute the nominal p-value associated with a GSEA for a collection of
gene sets, from the outputs of collectionGsea
}
\usage{
permutationPvalueCollectionGsea(permScores, dataScores)
}
\arguments{
  \item{permScores}{
a numeric matrix of permutation-based scores (rows=gene sets, named)
output from function 'collectionGsea' 
}
  \item{dataScores}{
a named numeric vector of observed scores output from function
'collectionGsea' 
}
}
\value{
a vector of p-values (named), one for each gene set
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
collectionGsea
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
GSCpvals <- permutationPvalueCollectionGsea(permScores=
GSCscores$Permutation.scores, dataScores=GSCscores$Observed.scores)
##example 2
\dontrun{
library(KEGG.db)
library(org.Dm.eg.db)
##load phenotype vector (see the vignette for details about the
##preprocessing of this data set)
data("data4enrich")
DM_KEGG <- KeggGeneSets(species="Dm")
GSCscores <- collectionGsea(collectionOfGeneSets=DM_KEGG, geneList=
data4enrich, exponent=1, nPermutations=1000, minGeneSetSize=100)
GSCpvals <- permutationPvalueCollectionGsea(permScores=
GSCscores$Permutation.scores, dataScores=GSCscores$Observed.scores)
}
}
