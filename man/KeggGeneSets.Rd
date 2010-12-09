\name{KeggGeneSets}
\alias{KeggGeneSets}
\title{
Create a list of gene sets based on KEGG pathways terms
}
\description{
This function creates a list of gene sets based on KEGG pathways terms. 
It is species-specific, and returns a list of gene sets, each of which 
is a character vector of Entrez gene identifiers.
}
\usage{
KeggGeneSets(species="Dm")
}
\arguments{
  \item{species}{
a single character value specifying the species: "Dm" ("Drosophila_
melanogaster"), "Hs" ("Homo_sapiens"), "Rn" ("Rattus_norvegicus"), "Mm"
("Mus_musculus") or "Ce" ("Caenorhabditis_elegans"))
}
}
\details{
This function relies on the following packages: GSEABase, KEGG.db.
}
\value{
a list of gene sets, with names as KEGG pathway IDs. Each gene set 
is a character vector of Entrez gene identifiers.
}
\author{
Camille Terfve, Xin Wang
}
\seealso{
\code{\link[HTSanalyzeR:GOGeneSets]{GOGeneSets}}
}
\examples{
library(org.Dm.eg.db)
library(KEGG.db)
DM_KEGG<-KeggGeneSets(species = "Dm")
}
