\name{plotEnrichMap}
\alias{plotEnrichMap}
\alias{plotEnrichMap,GSCA-method}
\title{
Plot and save an enrichment map for results of GSEA or hypergeometric tests
}
\description{
This is a generic function.

When implemented as the S4 method for objects of class \code{\link[HTSanalyzeR:GSCA]{GSCA}}, this function 
will plot and save an enrichment map for GSEA or Hypergeometric test results.  

To use this function for objects of class \code{\link[HTSanalyzeR:GSCA]{GSCA}}:

plotEnrichMap(object, resultName="GSEA.results", gscs, ntop=NULL, 
allSig=TRUE, gsNameType="id", displayEdgeLabel=TRUE, layout=
"layout.fruchterman.reingold", filepath=".", filename="test.png", 
output="png", ...)
}
\usage{
plotEnrichMap(object, ...)
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
	\item{resultName:}{
a single character value: 'HyperGeo.results' or 'GSEA.results'
} 
	\item{gscs:}{
a character vector specifying the names of gene set collections of which
the top significant gene sets will be plotted
}
	\item{ntop:}{
a single integer or numeric value indicating how many gene sets of top 
significance will be plotted.  
}
	\item{allSig:}{
a single logical value. If 'TRUE', all significant gene sets (adjusted 
p-value < 'pValueCutoff' of slot 'para') will be used; otherwise, only 
top 'ntop' gene sets will be used.
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
	\item{filepath:}{
a single character value specifying where to store the enrichment map.
}
	\item{output:}{
a single character value specifying the format of output image: "pdf" or
"png"
}
	\item{...}{
other arguments used by the function \code{png} or \code{pdf} such as
'width' and 'height'
}
}

}
\details{
See help(viewEnrichMap) for more details about the enrichment map for GSEA. 
}

\seealso{
\code{\link[HTSanalyzeR:viewEnrichMap]{viewEnrichMap}}
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
plotEnrichMap(KcViab_GSCA, gscs=c("GO_MF"), allSig=TRUE, ntop=NULL, gsNameType="id", 
displayEdgeLabel=FALSE,layout="layout.fruchterman.reingold", filepath="~", 
filename="GO_MF.pdf",output="pdf", width=8, height=8)
}
}











