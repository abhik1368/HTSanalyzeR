\name{GOGeneSets}
\alias{GOGeneSets}
\title{
Create gene set collections based on GO terms
}
\description{
This function creates gene set collections based on GO terms. It is species-specific, and returns a GeneSetCollection objects with the elements of the gene sets represented by Entrez Gene IDs.
}
\usage{
GOGeneSets(species = c("Drosophila_melanogaster", "Homo_sapiens", "Rattus_norvegicus", "Mus_musculus", "Caenorhabditis_elegans"),
 ontologies = c("BP", "MF", "CC"))
}
\arguments{
  \item{species}{
A choice of species (only one species of the following: human, mouse, rat,drosophila,celegans)
}
  \item{ontologies}{
A choice of ontology (multiple allowed): BP, CC, MF
}
}
\details{
This function relies on the following packages:GSEABase,GO.db,and either org.Hs.eg.db, org.Mm.eg.db, org.Rn.eg.db, org.Ce.eg.db, org.Dm.eg.db.
}
\value{
A GeneSetCollection object, with names=GO IDs and genes as Entrez gene IDs
}
\author{
Camille Terfve
}
\examples{
library(GO.db)
library(org.Dm.eg.db)
Dm.GO.CC<-GOGeneSets(species="Drosophila_melanogaster",ontologies=c("CC"))
}
