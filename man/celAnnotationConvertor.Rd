\name{celAnnotationConvertor}
\alias{celAnnotationConvertor}
\title{
Convert between different types of annotations for Caenorhabditis Elegans
}
\description{
This function converts an initial named data vector to the same vector 
but with a different identifier category. This function can also take a 
matrix, with rows=gene id's. This function removes the genes for which 
no mapping were found.
}
\usage{
celAnnotationConvertor(geneList, initialIDs = "Entrez.gene", finalIDs = 
"Entrez.gene", keepMultipleMappings = TRUE, verbose = TRUE)
}

\arguments{
  \item{geneList}{
a named integer or numeric vector, or a matrix with rows named by gene 
identifiers
}
  \item{initialIDs}{
a single character value specifying the type of initial identifiers for 
input 'geneList'. Current version can take one of the following types: 
"Ensembl.transcript", "Ensembl.prot", "Ensembl.gene", "Entrez.gene", 
"RefSeq", "Symbol", "GenBank" and "wormbase"
}
  \item{finalIDs}{
a single character value specifying the type of final identifiers to 
which users want to convert. Current version can take one of the following 
types: "Ensembl.transcript", "Ensembl.prot", "Ensembl.gene", "Entrez.gene", 
"RefSeq", "Symbol", "GenBank" and "wormbase"
}
  \item{keepMultipleMappings}{
a single logical value. If TRUE, the function keeps the entries with 
multiple mappings (first mapping is kept). If FALSE, the entries with 
multiple mappings will be discarded.
}
  \item{verbose}{
a single logical value deciding to display detailed messages (when 
verbose=TRUE) or not (when verbose=FALSE)
}
}
\details{
This function relies on the org.Ce.eg.db package and therefore only maps 
\itemize{
\item from any identifier to an Entrez gene id
\item from an Entrez gene ID to any identifier
}	
The entries that could not be mapped to any identifiers are removed from 
the resulting data vector/matrix.
}
\value{
the same vector/matrix but with names/rownames corresponding to a 
different type of identifiers.
}
\author{
Camille Terfve, Xin Wang
}
\seealso{
\code{\link[HTSanalyzeR:drosoAnnotationConvertor]{drosoAnnotationConvertor}}, \code{\link[HTSanalyzeR:mammalAnnotationConvertor]{mammalAnnotationConvertor}}, \code{\link[HTSanalyzeR:annotationConvertor]{annotationConvertor}}
}
\examples{
library(org.Ce.eg.db)
##example 1: convert a named vector
x <- runif(10)
names(x) <- names(as.list(org.Ce.egSYMBOL2EG))[1:10]
xEntrez <- celAnnotationConvertor(geneList=x, initialIDs="Symbol", 
finalIDs="Entrez.gene")
##example 2: convert a data matrix with row names as gene ids
x <- cbind(runif(10), runif(10))
rownames(x) <- names(as.list(org.Ce.egSYMBOL2EG))[1:10]
xEntrez <- celAnnotationConvertor(geneList=x, initialIDs="Symbol", 
finalIDs="Entrez.gene")
}
