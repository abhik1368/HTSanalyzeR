\name{drosoAnnotationConvertor}
\alias{drosoAnnotationConvertor}
\title{
Conversion between different types of identifiers.
}
\description{
This function converts an initial named data vector (or data matrix with row names=identifiers) to the same vector but with a different identifier category.
}
\usage{
drosoAnnotationConvertor(geneList, initialIDs = c("Ensembl.transcript", 
"Ensembl.prot", "Ensembl.gene", "Entrez.gene", "RefSeq", "Symbol", "GenBank", 
"Flybase", "FlybaseCG", "FlybaseProt"), finalIDs = c("Ensembl.transcript", 
"Ensembl.prot", "Ensembl.gene", "Entrez.gene", "RefSeq", "Symbol", "GenBank", 
"Flybase", "FlybaseCG", "FlybaseProt"), keepMultipleMappings = TRUE, verbose=TRUE)
}
\arguments{
  \item{geneList}{
A named data vector (names= identifiers to be converted) or a matrix with identifiers to be converted as row names.
}
  \item{initialIDs}{
The current identifiers type. One of "Ensembl.transcript", "Ensembl.prot", "Ensembl.gene", "Entrez.gene", "RefSeq", "Symbol", "GenBank", "Flybase", "FlybaseCG", "FlybaseProt".
}
  \item{finalIDs}{
The final identifier type that is wanted. One of "Ensembl.transcript", "Ensembl.prot", "Ensembl.gene", "Entrez.gene", "RefSeq", "Symbol", "GenBank", "Flybase", "FlybaseCG", "FlybaseProt".
}
  \item{keepMultipleMappings}{
If TRUE, keeps the entries with multiple mappings (first mapping is kept). If FALSE, discarding them.
}
  \item{verbose}{
A logical value to decide to display detailed messages (when verbose=TRUE) or not (when verbose=FALSE)
}
}
\details{
The entries that could not be mapped to any identifiers are removed from the resulting data vector/matrix.
This function relies on the org.Dm.eg.dbpackage and therefore only maps 
\itemize{
	\item from any identifier to an Entrez gene id
	\item from an Entrez gene ID to any identifier	
	}
}
\value{
The same data vector/matrix but with names/row names converted.
}
\author{
Camille Terfve
}
\seealso{
mammalAnnotationConvertor, celAnnotationConvertor
}
\examples{
library(org.Dm.eg.db)
x<-runif(10)
names(x)<-names(as.list(org.Dm.egSYMBOL2EG))[1:10]
xEntrez<-drosoAnnotationConvertor(geneList=x,initialIDs="Symbol",
finalIDs="Entrez.gene")
#
library(org.Dm.eg.db)
x<-cbind(runif(10),runif(10))
rownames(x)<-names(as.list(org.Dm.egSYMBOL2EG))[1:10]
xEntrez<-drosoAnnotationConvertor(geneList=x,initialIDs="Symbol",
finalIDs="Entrez.gene")
}
