\name{drosoAnnotationConvertor}
\alias{drosoAnnotationConvertor}
\title{
Convert between different types of gene identifiers for Drosophila Melanogaster
}
\description{
This function converts an initial data vector named by non-entrez ids to the same 
vector but with entrez ids, and vice versa. Genes for which no mapping were found 
will be removed. This function can also take a matrix, with gene identifiers as 
row names.
}
\usage{
drosoAnnotationConvertor(geneList, initialIDs = "Entrez.gene", finalIDs = 
"Entrez.gene", keepMultipleMappings = TRUE, verbose=TRUE)
}
\arguments{
  \item{geneList}{
a named integer or numeric vector, or a matrix with rows named by gene
identifiers
}
  \item{initialIDs}{
a single character value specifying the type of initial identifiers for
input 'geneList'. The current version can take one of the following types:
"Ensembl.transcript", "Ensembl.prot", "Ensembl.gene", "Entrez.gene",
"RefSeq", "Symbol", "GenBank", "Flybase", "FlybaseCG", "FlybaseProt".
}
  \item{finalIDs}{
a single character value specifying the type of final identifiers to
which users want to convert. The current version can take one of the 
following types: "Ensembl.transcript", "Ensembl.prot", "Ensembl.gene", 
"Entrez.gene", "RefSeq", "Symbol", "GenBank", "Flybase", "FlybaseCG", 
"FlybaseProt".
}
  \item{keepMultipleMappings}{
a single logical value. If TRUE, the function keeps the entries with
multiple mappings (first mapping is kept). If FALSE, the entries with
multiple mappings will be discarded.
}
  \item{verbose}{
a single logical value specifying to display detailed messages (when
verbose=TRUE) or not (when verbose=FALSE)
}
}
\details{
The entries that could not be mapped to any identifiers are removed from
the resulting data vector/matrix.
This function relies on the org.Dm.eg.dbpackage and therefore only maps 
\itemize{
	\item from any identifier to an Entrez gene id, or
	\item from an Entrez gene ID to any identifier	
	}
}
\value{
the same data vector/matrix but with names/row names converted.
}
\author{
Camille Terfve, Xin Wang
}
\seealso{
\code{\link[HTSanalyzeR:mammalAnnotationConvertor]{mammalAnnotationConvertor}}, \code{\link[HTSanalyzeR:celAnnotationConvertor]{celAnnotationConvertor}}, \code{\link[HTSanalyzeR:annotationConvertor]{annotationConvertor}}
}
\examples{
library(org.Dm.eg.db)
##example 1: convert a named vector
x <- runif(10)
names(x) <- names(as.list(org.Dm.egSYMBOL2EG))[1:10]
xEntrez <- drosoAnnotationConvertor(geneList=x, initialIDs="Symbol",
finalIDs="Entrez.gene")
##example 2: convert a data matrix with row names as gene ids
x <- cbind(runif(10), runif(10))
rownames(x) <- names(as.list(org.Dm.egSYMBOL2EG))[1:10]
xEntrez <- drosoAnnotationConvertor(geneList=x, initialIDs="Symbol",
finalIDs="Entrez.gene")
}
