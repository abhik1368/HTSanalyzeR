\name{gseaPlots}
\alias{gseaPlots}
\title{
Plot GSEA results for one gene set
}
\description{
This function takes in the output of gseaScores and the gene list, and
plot a figure of the ranked phenotypes, positions of the gene set and
the location of the enrichment score.  
}
\usage{
gseaPlots(runningScore, enrichmentScore, positions, geneList)
}
\arguments{
  \item{runningScore}{
a single numeric value specifying the enrichment score output from the
function "gseaScores"
}
  \item{enrichmentScore}{
a numeric vector of running enrichment scores (only in mode "graph")
}
  \item{positions}{
a numeric vector specifying positions in the gene list of the genes in
the gene set (only in mode "graph")
}
  \item{geneList}{
a numeric or integer vector of phenotypes in descending or ascending
order with elements named by their EntrezIds (no duplicates nor NA values)
}
}
\references{
Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L.,
Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander,
E. S. & Mesirov, J. P. (2005)  \emph{Gene set enrichment analysis: A
knowledge-based approach for interpreting genome-wide expression profiles.}
Proc. Natl. Acad. Sci. USA 102, 15545-15550.
}
\author{
Camille Terfve, Xin Wang
}
\seealso{
\code{\link[HTSanalyzeR:plotGSEA]{plotGSEA}}, \code{\link[HTSanalyzeR:viewGSEA]{viewGSEA}}
}
\examples{
\dontrun{
library(org.Dm.eg.db)
library(KEGG.db)
##load phenotype vector (see the vignette for details about the
##preprocessing of this data set)
data("KcViab_Data4Enrich")
DM_KEGG <- KeggGeneSets(species="Dm")
GSCscores <- gseaScores(geneList=KcViab_Data4Enrich, geneSet=DM_KEGG[[1]],
exponent=1, mode="graph")
gseaPlots(runningScore=GSCscores$runningScore, enrichmentScore=
GSCscores$enrichmentScore, positions=GSCscores$positions,
geneList=KcViab_Data4Enrich)
}
}
