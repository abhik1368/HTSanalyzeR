\name{multiHyperGeoTest}
\alias{multiHyperGeoTest}
\title{
 Hypergeometric test on gene set collections
}
\description{
This function performs hypergeometric tests for over-representation of
hits, on a list of gene sets. This function applies the hyperGeoTest
function to an entire list of gene sets and returns a data frame.
}
\usage{
multiHyperGeoTest(collectionOfGeneSets, universe, hits, minGeneSetSize =
15, pAdjustMethod = "BH", verbose = TRUE)
}
\arguments{
  \item{collectionOfGeneSets}{
a list of gene sets, each of which is a character vector of gene identifiers
}
  \item{universe}{
a character vector of all gene identifiers (usually all genes tested in
a screen)
}
  \item{hits}{
a character vector of gene identifiers for those genes considered as hits
}
  \item{minGeneSetSize}{
a single integer or numeric value specifying the minimum number of elements
in a gene set that must map to elements of the gene universe. Gene sets
with fewer genes than this number are removed from both hypergeometric
analysis and GSEA.
}
  \item{pAdjustMethod}{
a single character value specifying the p-value adjustment method to be
used (see 'p.adjust' for details)
}
  \item{verbose}{
a single logical value deciding to display detailed messages (when
verbose=TRUE) or not (when verbose=FALSE)
}
}
\value{
a dataframe containing the results of the hypergeometric test (one row
per gene set)
}
\author{
John C. Rose, Xin Wang
}
\seealso{
hyperGeoTest
}
\examples{
##example 1
gl <- runif(100, min=0, max=5)
gl <- gl[order(gl, decreasing=TRUE)]
names(gl) <- as.character(sample(x=seq(from=1, to=100, by=1), size=100,
replace=FALSE))
gs1 <- sample(names(gl), size=20, replace=FALSE)
gs2 <- sample(names(gl), size=20, replace=FALSE)
gsc <- list(subset1=gs1, subset2=gs2)
hypgeo<-multiHyperGeoTest(collectionOfGeneSets=gsc, universe=names(gl), 
hits=names(gl)[which(abs(gl) > 2)], minGeneSetSize = 2, pAdjustMethod ="BH")
##example 2
\dontrun{
library(org.Dm.eg.db)
library(KEGG.db)
##load phenotype vector (see the vignette for details about the
##preprocessing of this data set)
data("data4enrich")
DM_KEGG <- KeggGeneSets(species="Dm")
##Do multiple hypergeometric tests
hypgeoResults <- multiHyperGeoTest(collectionOfGeneSets=DM_KEGG,
universe=names(data4enrich), hits=names(data4enrich)[which(abs(
data4enrich) > 2)], minGeneSetSize = 15, pAdjustMethod = "BH")
}
}
