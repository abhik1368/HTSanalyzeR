\name{viewEnrichMap}
\alias{viewEnrichMap}
\alias{viewEnrichMap,GSCA-method}
\title{
Plot on the screen a figure of the enrichment map for GSEA results
}
\description{
This is a generic function.

When implemented as the S4 method for objects of class \code{\link[HTSanalyzeR:GSCA]{GSCA}}, this function 
will plot and save an enrichment map for GSEA results.  

To use this function for objects of class \code{\link[HTSanalyzeR:GSCA]{GSCA}}:

viewEnrichMap(object, gscs, ntop=NULL, allSig=TRUE, gsNameType="id", 
displayEdgeLabel=TRUE, layout="layout.fruchterman.reingold")
}
\usage{
viewEnrichMap(object, ...)
}
\arguments{
	\item{object}{
an object. When this function is implemented as the S4 method of class \code{\link[HTSanalyzeR:GSCA]{GSCA}}, 
this argument is an object of class \code{\link[HTSanalyzeR:GSCA]{GSCA}}.
}
	\item{...}{
other arguments. (see below for the arguments supported by the method of
class \code{\link[HTSanalyzeR:GSCA]{GSCA}})
}
\describe{
	\item{gscs:}{
a character vector specifying the names of gene set collections of which
the top significant gene sets will be plot
}
	\item{ntop:}{
a single integer or numeric value specifying to plot how many gene sets
of top significance.  
}
	\item{allSig:}{
a single logical value. If 'TRUE', all significant gene sets (GSEA
adjusted p-value < 'pValueCutoff' of slot 'para') will be used; otherwise,
only top 'ntop' gene sets will be used.
}
	\item{gsNameType:}{
a single character value specifying the type of the gene set names that 
will be displayed as the names of nodes in the enrichment map. The type of
the gene set names should be one of the following: "id", "term" or "none".
}
	\item{displayEdgeLabel:}{
a single logical value specifying whether or not to display the labels of 
the edges in the enrichment map
}
	\item{layout:}{
a single character value specifying the layout of the enrichment map. 
(see help(layout) of the package \code{igraph})
}
}

}
\details{
The idea of this function is similar to the PLoS one paper by Merico et al. 

An enrichment map is a network to help better visualize and interpret the 
GSEA results. In an enrichment map, the nodes represent gene sets and the 
edges denote the Jaccard similarity coefficient of two gene sets. Nodes are 
colored by the sign of the enrichment scores (red:+, blue: -) and scaled 
according to the adjusted GSEA p-values (the darker the more significant). 
The size of nodes illustrates the size of gene sets, while the width of 
edges denotes the Jaccard coefficient.    
}
\value{
an object of \code{igraph} with all attributes about the enrichement map
}
\seealso{
\code{\link[HTSanalyzeR:plotEnrichMap]{plotEnrichMap}}
}
\references{
Merico, D. et al, Enrichment Map: A Network-Based Method for Gene-Set 
Enrichment Visualization and Interpretation, PloS one, 2010, e13984
}
\author{
Xin Wang \email{xw264@cam.ac.uk}
}

\examples{
\dontrun{
library(org.Dm.eg.db)
library(KEGG.db)
##load data for enrichment analyses
data("KcViab_GSCA")
##plot and save the enrichment map
viewEnrichMap(KcViab_GSCA, gscs=c("GO_MF"), allSig=TRUE, ntop=NULL, gsNameType="id", 
displayEdgeLabel=FALSE,layout="layout.fruchterman.reingold")
##append Gene set terms to results
KcViab_GSCA<-appendGSTerms(KcViab_GSCA, goGSCs=c("GO_BP","GO_MF","GO_CC"), 
keggGSCs=c("PW_KEGG"))
viewEnrichMap(KcViab_GSCA, gscs=c("GO_MF"), allSig=TRUE, ntop=NULL, gsNameType="term", 
displayEdgeLabel=FALSE,layout="layout.fruchterman.reingold")
}
}











