\name{gseaPlots}
\alias{gseaPlots}
\title{
Produce a collection of plots for GSEA
}
\description{
This function takes in the output of gseaScores and the gene list, and produce plots for GSEA.
}
\usage{
gseaPlots(runningScore, enrichmentScore, positions, geneList, output = "png", filepath=".", filename="test")
}
\arguments{
  \item{runningScore}{
a single numeric value specifying the enrichment score output from the function "gseaScores"
}
  \item{enrichmentScore}{
a numeric vector of running enrichment scores (only in mode "graph")
}
  \item{positions}{
a numeric vector specifying positions in the gene list of the genes in the gene set (only in mode "graph")
}
  \item{geneList}{
a numeric or integer vector of phenotypes in descending or ascending order with elements named by their EntrezIds (no duplicates nor NA values)
}
  \item{output}{
a single character value specifying the format of output image: "pdf" or "png"
}
  \item{filepath}{
a single character value specifying the directory where these plots will be stored
}
  \item{filename}{
a single character value specifying the name of the gene set for which these plots are produced
}
}
\details{
The plots are going to be produced and stored in the directory 'filepath' with the name 'filename.png' or 'filename.pdf'.
}
\value{
This function does not return anything, it just produces an image file containing the plots.
}
\references{
Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S. & Mesirov, J. P. (2005)  \emph{Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles.} Proc. Natl. Acad. Sci. USA 102, 15545-15550.
}
\author{
Camille Terfve, Xin Wang
}
\examples{
\dontrun{
library(org.Dm.eg.db)
library(KEGG.db)
##load phenotype vector (see the vignette for details about the preprocessing of this data set)
data("data4enrich")
kegg.droso<-KeggGeneSets(species="Dm")
test<-gseaScores(geneList=data4enrich,geneSet=kegg.droso[[1]],exponent=1,mode="graph")
gseaPlots(runningScore=test$runningScore, enrichmentScore=test$enrichmentScore, 
positions=test$positions,geneList=data4enrich, filepath=".", filename=names(kegg.droso)[1])
}
}
