\name{celAnnotationConvertor}
\alias{celAnnotationConvertor}
\title{
Conversion between identifier categories.
}
\description{
This function converts an initial named data vector to the same vector but with a different identifier category. This function also works on a matrix, with rows=gene id's. This function removes the genes for which no mapping was found.
}
\usage{
celAnnotationConvertor(geneList, initialIDs = c("Ensembl.transcript",
 "Ensembl.prot", "Ensembl.gene", "Entrez.gene", "RefSeq", "Symbol", "GenBank"), 
finalIDs = c("Ensembl.transcript", "Ensembl.prot", "Ensembl.gene", 
"Entrez.gene", "RefSeq", "Symbol", "GenBank", "wormbase"), 
keepMultipleMappings = TRUE)
}

\arguments{
  \item{geneList}{
Named data vector or data matrix with rows named by gene identifiers.
}
  \item{initialIDs}{
The type of initial identifiers wanted. One of: "Ensembl.transcript", "Ensembl.prot", "Ensembl.gene", "Entrez.gene", "RefSeq", "Symbol", "GenBank"
}
  \item{finalIDs}{
The type of final identifiers wanted. One of "Ensembl.transcript", "Ensembl.prot", "Ensembl.gene", "Entrez.gene", "RefSeq", "Symbol", "GenBank", "wormbase"
}
  \item{keepMultipleMappings}{
If TRUE, the function keeps the entries with multiple mappings (first mapping is kept). If FALSE, the entries with multiple mappings will be discarded.
}
}
\details{
This function relies on the org.Ce.eg.db package and therefore only maps 
\itemize{
\item from any identifier to an Entrez gene id
\item from an Entrez gene ID to any identifier
}	
The entries that could not be mapped to any identifiers are removed from the resulting data vector/matrix.
}
\value{
The same vector-matrix but with names-rownames corresponding to a different type of identifiers.
}
\author{
Camille Terfve
}
\seealso{
drosoAnnotationConvertor, mammalAnnotationConvertor
}
\examples{
library(org.Ce.eg.db)
x<-runif(10)
names(x)<-names(as.list(org.Ce.egSYMBOL2EG))[1:10]
xEntrez<-celAnnotationConvertor(geneList=x,initialIDs="Symbol",finalIDs="Entrez.gene")
#
library(org.Ce.eg.db)
x<-cbind(runif(10),runif(10))
rownames(x)<-names(as.list(org.Ce.egSYMBOL2EG))[1:10]
xEntrez<-celAnnotationConvertor(geneList=x,initialIDs="Symbol",finalIDs="Entrez.gene")
}
