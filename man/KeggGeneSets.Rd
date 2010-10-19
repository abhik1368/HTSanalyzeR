\name{KeggGeneSets}
\alias{KeggGeneSets}
\title{
This function creates gene set collections based on Kegg pathways terms.
}
\description{
This function creates gene set collections based on Kegg pathways terms. It is species-specific, and returns a GeneSetCollection objects with the elements of the gene sets represented by Entrez Gene IDs.
}
\usage{
KeggGeneSets(species=c("Dm", "Hs", "Rn", "Mm", "Ce"))
}
\arguments{
  \item{species}{
A choice of species (only one species of the following: "Dm" ("Drosophila_melanogaster"), "Hs" ("Homo_sapiens"), 
"Rn" ("Rattus_norvegicus"), "Mm" ("Mus_musculus"), "Ce" ("Caenorhabditis_elegans"))
}
}
\details{
This function relies on the following packages:GSEABase,KEGG.db.
}
\value{
A GeneSetCollection object, with names=kegg pathway IDs and genes as entrez gene IDs
}
\author{
Camille Terfve
}
\seealso{
GOGeneSets
}
\examples{
library(org.Dm.eg.db)
library(KEGG.db)
kegg.droso<-KeggGeneSets(species = "Dm")
}
