\name{mammalAnnotationConvertor}
\alias{mammalAnnotationConvertor}
\title{
This function converts different identifier categories.
}
\description{
This function converts an initial named data vector to the same vector but with a different identifier category (this function also works on matrices with rows=genes (named)).
}
\usage{
mammalAnnotationConvertor(geneList, initialIDs = c(
"Ensembl.transcript", "Ensembl.prot", "Ensembl.gene", "Entrez.gene", 
"RefSeq", "Symbol", "GenBank"), finalIDs = c("Ensembl.transcript", 
"Ensembl.prot", "Ensembl.gene", "Entrez.gene", "RefSeq", "Symbol", 
"GenBank"), species = c("Homo_sapiens", "Rattus_norvegicus", 
"Mus_musculus"), keepMultipleMappings = TRUE)
}
\arguments{
  \item{geneList}{
A named data vector or data matrix (row names= gene IDs)
}
  \item{initialIDs}{
One of "Ensembl.transcript", "Ensembl.prot", "Ensembl.gene", "Entrez.gene", "RefSeq", "Symbol", "GenBank"
}
  \item{finalIDs}{
One of "Ensembl.transcript", "Ensembl.prot", "Ensembl.gene", "Entrez.gene", "RefSeq", "Symbol", "GenBank"
}
  \item{species}{
One of "Homo_sapiens","Rattus_norvegicus","Mus_musculus"
}
 \item{keepMultipleMappings}{
If TRUE the function is keeping the entries with multiple mappings (first mapping is kept), if  FALSE it is discarding them
}
}
\details{
This function removes the genes for which no mapping was found.This function relies on the org.Hs.eg.db/org.Mm.eg.db/org.Rn.eg.db packages and therefore only maps from any identifier to an Entrez gene id or from an Entrez gene ID to any identifier.
}
\value{
The same data vector-matrix but with another type of identifiers as names-row names.
}
\author{
Camille Terfve
}
\examples{
library(org.Hs.eg.db)
x<-runif(10)
names(x)<-names(as.list(org.Hs.egSYMBOL2EG))[1:10]
xEntrez<-mammalAnnotationConvertor(geneList=x,initialIDs="Symbol",
finalIDs="Entrez.gene", species="Homo_sapiens")
#
library(org.Hs.eg.db)
x<-cbind(runif(10),runif(10))
rownames(x)<-names(as.list(org.Hs.egSYMBOL2EG))[1:10]
xEntrez<-mammalAnnotationConvertor(geneList=x,initialIDs="Symbol",
finalIDs="Entrez.gene", species="Homo_sapiens")
}
