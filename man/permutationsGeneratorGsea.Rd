\name{permutationsGeneratorGsea}
\alias{permutationsGeneratorGsea}
\title{
Generate random permutations of a gene list
}
\description{
This function generates random permutations of the identifiers in a gene list.  It is used in the GSEA process.
}
\usage{
permutationsGeneratorGsea(nPermutations, geneList)
}
\arguments{
  \item{nPermutations}{
Number of permutations to be made
}
  \item{geneList}{
Ordered gene list as a named vector, with names=gene IDs 
}
}
\value{
Matrix of dimension nPermutations*geneListIDs:each row is a random permutation of the gene identifiers in the gene list
}
\author{
Camille Terfve
}
\seealso{
gseaScores, permutationPvalueGsea,permutationPvalueCollectionGsea
}
\examples{
gl<-runif(100,min=0,max=5)
gl<-order(gl)
names(gl)<-as.character(sample(x=seq(from=1,to=100,by=1),size=100,replace=FALSE))
perm<-permutationsGeneratorGsea(nPermutations=1000,geneList=gl)
}
