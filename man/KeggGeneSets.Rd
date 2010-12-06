\name{KeggGeneSets}
\alias{KeggGeneSets}
\title{
Create gene set collections (lists) based on Kegg pathways terms
}
\description{
This function creates gene set collections (lists) based on Kegg pathways
terms. It is species-specific, and returns a list of gene set collections
(lists) with the elements of the gene sets represented by Entrez Gene IDs.
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
a collection (list) of gene sets, with names as kegg pathway IDs and
genes as entrez gene IDs
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
