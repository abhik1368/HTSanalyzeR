\name{viewGSEA}
\alias{viewGSEA}
\alias{viewGSEA,GSCA-method}
\title{
Plot on the screen a figure of GSEA results of one gene set
}
\description{
This is a generic function.

When implemented as the S4 method for objects of class \code{\link[HTSanalyzeR:GSCA]{GSCA}}, this 
function plots a figure of the positions of the gene sets in the ranked
phenotypes and the location of the enrichment score. 

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
a single character value specifying the name of the gene set to be plot
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
##do preprocessing (data4enrich has already been preprocessed)
gsca <- preprocess(gsca, species="Dm", initialIDs = "Entrez.gene", 
keepMultipleMappings = TRUE, duplicateRemoverMethod = "max", 
orderAbsValue = FALSE)
##print summary of gsca again
summarize(gsca)
##do hypergeometric tests and GSEA
gsca <- analyze(gsca, para = list(pValueCutoff = 0.05, pAdjustMethod 
= "BH", nPermutations = 1000, minGeneSetSize = 100,exponent = 1))
##print summary of results
summarize(gsca, what="Result")
##print top significant gene sets in GO.BP
topPWKEGG<-getTopGeneSets(gsca, "GSEA.results", "PW_KEGG", allSig=TRUE)
##view a gene set
viewGSEA(gsca, "PW_KEGG", topPWKEGG[["PW_KEGG"]][1])
}
}











