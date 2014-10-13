\name{GOGeneSets}
\alias{GOGeneSets}
\title{
Create a list of gene sets based on GO terms
}
\description{
This function creates a list of gene sets based on GO terms. It is
species-specific, and returns a list of gene sets, each of which 
is a character vector of Entrez identifiers.
}
\usage{
GOGeneSets(species = "Dm", ontologies = "MF")
}
\arguments{
  \item{species}{
a single character value specifying a choice of species: "Dm" ("Drosophila_
melanogaster"), "Hs" ("Homo_sapiens"), "Rn" ("Rattus_norvegicus"), "Mm"
("Mus_musculus") or "Ce" ("Caenorhabditis_elegans"))
}
  \item{ontologies}{
a single character value or a character vector specifying an ontology or
multiple ontologies. The current version provides the following choices: 
"BP", "CC" and "MF"
}
}
\details{
This function relies on the following packages: GSEABase, GO.db, and either
org.Hs.eg.db, org.Mm.eg.db, org.Rn.eg.db, org.Ce.eg.db, org.Dm.eg.db.
}
\value{
a list of gene sets, with names as GO IDs. Each gene set is a character 
vector of Entrez identifiers.
}
\author{
Camille Terfve
}
\seealso{
\code{\link[HTSanalyzeR:KeggGeneSets]{KeggGeneSets}}
}
\examples{
library(GO.db)
library(org.Dm.eg.db)
Dm_GO_CC<-GOGeneSets(species="Dm",ontologies=c("CC"))
}
