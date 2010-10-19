\name{overlapComputing}
\alias{overlapComputing}
\title{
Overlap between gene list and gene set
}
\description{
This function computes the overlap between a gene list and a gene set.  This might be useful e.g. to check that the identifier mapping was done correctly.
}
\usage{
overlapComputing(geneList, geneSet)
}
\arguments{
  \item{geneList}{
A gene list as a single named vector, with duplicates removed.
}
  \item{geneSet}{
A gene set as a vector of identifiers
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
gs1<-sample(names(gl),size=20,replace=FALSE)
gs2<-sample(names(gl),size=20,replace=FALSE)
gsc<-list(subset1=gs1,subset2=gs2)
overlapComputing(geneList=gl,geneSet=gsc[["subset1"]])
}
