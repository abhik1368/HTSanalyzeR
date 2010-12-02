\name{hyperGeoTest}
\alias{hyperGeoTest}
\title{
Performs hypergeometric tests for over-representation analysis
}
\description{
This function takes in a single gene set (GeneSet), a vector (GeneList)
of gene symbols for all tested genes, a vector of "hits" (hits), and a
p-value adjustment method. It outputs a vector containing the size of
the gene universe, the size of the gene set within this universe (i.e.
how many genes from the universe map to this gene set), the total number
of hits, the number of hits expected to occur in the gene set, the actual
hits observed in the gene set, and the pvalue from a hypergeometric test.
}
\usage{
hyperGeoTest(geneSet, universe, hits)
}
\arguments{
  \item{geneSet}{
a character vector specifying a gene set 
}
  \item{universe}{
a character vector of all gene identifiers (usually all genes tested in
a screen)
}
  \item{hits}{
a character vector of gene identifiers for those genes considered as hits
}
}
\value{
a numeric vector containing the size of the gene universe (named by
"Universe Size"), 
the size of the gene set within this universe (i.e. how many genes from
the universe map to this gene set) (named by "Gene Set Size"), 
the total number of hits (named by "Total Hits"), 
the number of hits expected to occur in the gene set (named by "Expected
Hits"), 
the actual hits observed in the gene set (named by "Observed Hits"), 
and the pvalue from a hypergeometric test (named by "Pvalue").
}
\author{
John C. Rose, Xin Wang
}
\seealso{
multiHyperGeoTest
}
\examples{
##example 1
gl <- runif(100, min=0, max=5)
gl <- order(gl)
names(gl) <- as.character(sample(x=seq(from=1, to=100, by=1), size=100,
replace=FALSE))
gs1 <- sample(names(gl), size=20, replace=FALSE)
gs2 <- sample(names(gl), size=20, replace=FALSE)
gsc <- list(subset1=gs1, subset2=gs2)
hypgeo <- hyperGeoTest(geneSet=gsc[["subset1"]], universe=names(gl), 
hits=names(gl)[which(abs(gl) > 2)])
##example 2 
library(org.Dm.eg.db)
library(KEGG.db)
##load phenotype vector (see the vignette for details about the
##preprocessing of this data set)
data("data4enrich")
##Prepare kegg gene sets
DM_KEGG<-KeggGeneSets(species="Dm")
##Do the tests
hypgeoResults <- hyperGeoTest(geneSet=DM_KEGG[[1]], universe=
names(data4enrich), hits=names(data4enrich)[which(abs(data4enrich) > 2)])
}
