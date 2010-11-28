\name{FDRcollectionGsea}
\alias{FDRcollectionGsea}
\title{
Compute the FDR associated with a GSEA for a collection (list) of gene sets
}
\description{
This function computes the FDR associated with a permutation-based p-value from the GSEA on a list of gene sets
}
\usage{
FDRcollectionGsea(permScores, dataScores)
}
\arguments{
  \item{permScores}{
a numeric matrix of permutation-based scores resulting from the output of collectionGsea
}
  \item{dataScores}{
a named numeric vector of observed scores resulting from the output of collectionGsea
}
}

\value{
a named numeric vector of FDR, one for each gene set
}
\references{
Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S. & Mesirov, J. P. (2005)  \emph{Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles.} Proc. Natl. Acad. Sci. USA 102, 15545-15550.
}
\author{
Camille Terfve, Xin Wang
}
\seealso{
collectionGsea, permutationPvalueCollectionGsea
}
\examples{
##example 1
gl<-runif(100, min=0, max=5)
gl<-gl[order(gl, decreasing=TRUE)]
names(gl)<-as.character(sample(x=seq(from=1, to=100, by=1), size=100, replace=FALSE))
gs1<-sample(names(gl), size=20, replace=FALSE)
gs2<-sample(names(gl), size=20, replace=FALSE)
gscs<-list(gs1=gs1, gs2=gs2)
scores.collection<-collectionGsea(collectionOfGeneSets=gscs, geneList=gl, exponent=1, nPermutation=1000, minGeneSetSize=5)
FDR<-FDRcollectionGsea(permScores=scores.collection$Permutation.scores,dataScores=scores.collection$Observed.scores)
##example 2 (see the vignette for details about the preprocessing of this data set)
\dontrun{
library(org.Dm.eg.db)
library(KEGG.db)
data("data4enrich")
kegg.droso<-KeggGeneSets(species="Dm")
scores.collection<-collectionGsea(collectionOfGeneSets=kegg.droso, geneList=data4enrich,exponent=1,nPermutations=1000,minGeneSetSize=150)
FDRcollectionkegg<-FDRcollectionGsea(permScores=scores.collection$Permutation.scores,dataScores=scores.collection$Observed.scores)
}
}
