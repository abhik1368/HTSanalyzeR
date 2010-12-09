\name{mammalAnnotationConvertor}
\alias{mammalAnnotationConvertor}
\title{
Convert between different types of identifiers for mammalian species
}
\description{
This function converts an initial named data vector to the same vector
but with a different identifier category (this function also works on
matrices with gene identifiers as row names).
}
\usage{
mammalAnnotationConvertor(geneList, initialIDs = "Entrez.gene", finalIDs = 
"Entrez.gene", species = "Hs", keepMultipleMappings = TRUE, verbose=TRUE)
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
"RefSeq", "Symbol" and "GenBank".
}
  \item{finalIDs}{
a single character value specifying the type of final identifiers to
which users want to convert. The current version can take one of the 
following types: "Ensembl.transcript", "Ensembl.prot", "Ensembl.gene", 
"Entrez.gene", "RefSeq", "Symbol" and "GenBank"
}
  \item{species}{
a single character value specifying the species: "Hs" ("Homo_sapiens"),
"Rn" ("Rattus_norvegicus") or "Mm" ("Mus_musculus")
}
 \item{keepMultipleMappings}{
a single logical value. If TRUE, the function keeps the entries with
multiple mappings (first mapping is kept). If FALSE, the entries with
multiple mappings will be discarded.
}
 \item{verbose}{
a single logical value indicating to display detailed messages (when
verbose=TRUE) or not (when verbose=FALSE)
}
}
\details{
This function removes the genes for which no mapping was found. This
function relies on the org.Hs.eg.db/org.Mm.eg.db/org.Rn.eg.db packages
and therefore only maps from any identifier to an Entrez gene ID or from
an Entrez gene ID to any identifier.
}
\value{
the same data vector/matrix but with another type of identifiers as
names/row names
}
\author{
Camille Terfve, Xin Wang
}
\seealso{
\code{\link[HTSanalyzeR:drosoAnnotationConvertor]{drosoAnnotationConvertor}}, \code{\link[HTSanalyzeR:celAnnotationConvertor]{celAnnotationConvertor}}, \code{\link[HTSanalyzeR:annotationConvertor]{annotationConvertor}}
}
\examples{
##example 1: convert a named vector
library(org.Hs.eg.db)
x <- runif(10)
names(x) <- names(as.list(org.Hs.egSYMBOL2EG))[1:10]
xEntrez <- mammalAnnotationConvertor(geneList=x, initialIDs="Symbol",
finalIDs="Entrez.gene", species="Hs")
##example 2: convert a data matrix with row names as gene ids
library(org.Hs.eg.db)
x <- cbind(runif(10), runif(10))
rownames(x) <- names(as.list(org.Hs.egSYMBOL2EG))[1:10]
xEntrez <- mammalAnnotationConvertor(geneList=x, initialIDs="Symbol",
finalIDs="Entrez.gene", species="Hs")
}
