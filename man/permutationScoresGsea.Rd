\name{permutationScoresGsea}
\alias{permutationScoresGsea}
\title{
This function computes permutation-based scores for GSEA
}
\description{
This function computes permutation-based scores for GSEA, it is a modification of the permutation.pvalue function, that does not compute a pvalue but returns the scores of the permutationsMatrix. 
}
\usage{
permutationScoresGsea(permutationsMatrix, geneList, geneSet, exponent)
}
\arguments{
  \item{permutationsMatrix}{
A matrix of permutationsMatrix from permutationsMatrixgenerator
}
  \item{geneList}{
An ordered and named data vector (no duplicates nor NAs)
}
  \item{geneSet}{
A gene set in geneSet format
}
  \item{exponent}{
An exponent for gseascores
}
}
\value{
A vector of permutation-based scores
}
\references{
Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S. & Mesirov, J. P. (2005) \emph{ Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles.} Proc. Natl. Acad. Sci. USA 102, 15545-15550.
}
\author{
Camille Terfve
}
\seealso{
permutationPvalueGsea, permutationPvalueCollectionGsea, gseaScores
}
\examples{
gl<-runif(100,min=0,max=5)
gl<-order(gl)
names(gl)<-as.character(sample(x=seq(from=1,to=100,by=1),size=100,replace=FALSE))
gs1<-GeneSet(sample(names(gl),size=20,replace=FALSE), setName="subset1")
gs2<-GeneSet(sample(names(gl),size=20,replace=FALSE), setName="subset2")
gsc<-GeneSetCollection(c(gs1,gs2))
perm<-permutationsGeneratorGsea(nPermutations=1000,geneList=gl)
pscores<-permutationScoresGsea(permutationsMatrix=perm, geneList=gl, 
geneSet=geneIds(gsc[[2]]), exponent=1)
}
