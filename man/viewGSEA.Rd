\name{viewGSEA}
\alias{viewGSEA}
\alias{viewGSEA,GSCA-method}
\title{
Plot a figure of GSEA results for one gene set
}
\description{
This is a generic function.

When implemented as the S4 method for objects of class \code{\link[HTSanalyzeR:GSCA]{GSCA}}, this 
function plots a figure of the positions of the gene sets in the ranked
gene list and the location of the enrichment score. 

To use this function for objects of class \code{\link[HTSanalyzeR:GSCA]{GSCA}}:

viewGSEA(object, gscName, gsName)  
}
\usage{
viewGSEA(object, ...) 
}
\arguments{
	\item{object}{
an object. When this function is implemented as the S4 method of class 
\code{\link[HTSanalyzeR:GSCA]{GSCA}}, this argument is an object of class \code{\link[HTSanalyzeR:GSCA]{GSCA}}.
}
	\item{...}{
other arguments. (see below for the arguments supported by class \code{\link[HTSanalyzeR:GSCA]{GSCA}})
}
\describe{
	\item{gscName:}{
a single character value specifying the name of the gene set collection 
where the gene set is
}
	\item{gsName:}{
a single character value specifying the name of the gene set to be plotted
}
}
}
\details{
We suggest to print the names of top significant gene sets using the
function \code{\link[HTSanalyzeR:getTopGeneSets]{getTopGeneSets}} before
plotting the GSEA results.
}

\author{
Xin Wang \email{xw264@cam.ac.uk}
}

\seealso{
\code{\link[HTSanalyzeR:plotGSEA]{plotGSEA}}, \code{\link[HTSanalyzeR:gseaPlots]{gseaPlots}}
}

\examples{
\dontrun{
library(org.Dm.eg.db)
library(KEGG.db)
##load sample data
data("KcViab_GSCA")
##print summary of results
summarize(KcViab_GSCA, what="Result")
##print top significant gene sets in GO.BP
topPWKEGG<-getTopGeneSets(KcViab_GSCA, "GSEA.results", "PW_KEGG", allSig=TRUE)
##view a gene set
viewGSEA(KcViab_GSCA, "PW_KEGG", topPWKEGG[["PW_KEGG"]][1])
}
}











