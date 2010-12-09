\name{plotGSEA}
\alias{plotGSEA}
\alias{plotGSEA,GSCA-method}
\title{
Plot and save figures of GSEA results for top significant gene sets
}
\description{
This is a generic function.

When implemented as the S4 method for objects of class \code{\link[HTSanalyzeR:GSCA]{GSCA}}, this function 
plots figures of the positions of genes of the gene set in the ranked gene 
list and the location of the enrichment score for top significant gene sets.  

To use this function for objects of class \code{\link[HTSanalyzeR:GSCA]{GSCA}}:

plotGSEA(object, gscs, ntop=NULL, allSig=FALSE, filepath=".", output=
"png", ...)
}
\usage{
plotGSEA(object, ...)
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
a character vector specifying the names of gene set collections whose
top significant gene sets will be plotted
}
	\item{ntop:}{
a single integer or numeric value specifying how many gene sets of top 
significance will be plotted.  
}
	\item{allSig:}{
a single logical value. If 'TRUE', all significant gene sets (GSEA
adjusted p-value < 'pValueCutoff' of slot 'para') will be plotted; otherwise,
only top 'ntop' gene sets will be plotted.
}
	\item{filepath:}{
a single character value specifying where to store GSEA figures.
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
To make GSEA plots of top significance using this function, the user can
only choose one method: either assign an integer to the argument 'ntop'
or set the argument 'allSig' to 'TRUE'. Exceptions will occur if both
methods are used, or no method is used. Please also note that the argument
'ntop' is a cutoff for all gene set collections in the argument 'gscs'.

We suggest to perform \code{summarize(gsca, what="Result")} first to have 
an idea of how many significant gene sets there are, and then choose to 
plot them all or just the top ones.
}

\seealso{
\code{\link[HTSanalyzeR:viewGSEA]{viewGSEA}}, \code{\link[HTSanalyzeR:gseaPlots]{gseaPlots}}
}

\author{
Xin Wang \email{xw264@cam.ac.uk}
}

\examples{
\dontrun{
library(org.Dm.eg.db)
library(KEGG.db)
##load data for enrichment analyses
data("KcViab_Data4Enrich")
##select hits
hits <- names(KcViab_Data4Enrich)[which(abs(KcViab_Data4Enrich) > 2)]
##set up a list of gene set collections
PW_KEGG <- KeggGeneSets(species = "Dm")
gscList <- list(PW_KEGG = PW_KEGG)
##create an object of class 'GSCA'
gsca <- new("GSCA", listOfGeneSetCollections=gscList, geneList = 
KcViab_Data4Enrich, hits = hits)
##print summary of gsca
summarize(gsca)
##do preprocessing (KcViab_Data4Enrich has already been preprocessed)
gsca <- preprocess(gsca, species="Dm", initialIDs = "Entrez.gene", 
keepMultipleMappings = TRUE, duplicateRemoverMethod = "max", 
orderAbsValue = FALSE)
##print summary of gsca again
summarize(gsca)
##do hypergeometric tests and GSEA
gsca <- analyze(gsca, para = list(pValueCutoff = 0.05, pAdjustMethod 
= "BH", nPermutations = 1000, minGeneSetSize = 100, exponent = 1))
##print summary of results
summarize(gsca, what="Result")
##plot all significant gene sets
plotGSEA(gsca, gscs=c("PW_KEGG"), allSig=TRUE, filepath=".", output=
"pdf", width=8, height=8)
}
}











