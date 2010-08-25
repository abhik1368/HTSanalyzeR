\name{overlapComputing}
\alias{overlapComputing}
\title{
Overlap between gene list and gene set/collection
}
\description{
This function computes the overlap between a gene list and a gene set or a collection of gene sets.  This might be useful e.g. to check that the identifier mapping was done correctly.
}
\usage{
overlapComputing(geneList, geneSet)
}
\arguments{
  \item{geneList}{
A gene list as a single named vector, with duplicates removed.
}
  \item{geneSet}{
A gene set as a vector of identifiers or as a GeneSet object (GSEABase)
}
}
\value{
None, just prints on the screen the number of identifiers matched, their proportion etc.
}
\author{
Camille Terfve
}
\examples{
gl<-runif(100,min=0,max=5)
gl<-order(gl)
names(gl)<-as.character(sample(x=seq(from=1,to=100,by=1),size=100,replace=FALSE))
gs1<-GeneSet(sample(names(gl),size=20,replace=FALSE), setName="subset1")
gs2<-GeneSet(sample(names(gl),size=20,replace=FALSE), setName="subset2")
gsc<-GeneSetCollection(c(gs1,gs2))
overlapComputing(geneList=gl,geneSet=gsc[["subset1"]])
overlapComputing(geneList=gl,geneSet=gsc)
}
