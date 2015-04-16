\name{summarize}
\alias{summarize}
\alias{summarize,GSCA-method}
\alias{summarize,NWA-method}
\alias{show}
\alias{show,GSCA-method}
\alias{show,NWA-method}
\title{
Print summary information for an object of class GSCA or NWA
}
\description{
This is a generic function.

When implemented as the S4 method for objects of class \code{\link[HTSanalyzeR:GSCA]{GSCA}} or \code{\link[HTSanalyzeR:NWA]{NWA}}, 
this function prints a summary of information about the slots of these
classes. 

To use this function for objects of class \code{\link[HTSanalyzeR:GSCA]{GSCA}} or \code{\link[HTSanalyzeR:NWA]{NWA}}:

summarize(object, what="ALL")
}
\usage{
summarize(object, ...)
}
\arguments{
	\item{object}{
an object. When this function is implemented as the S4 method of class
\code{\link[HTSanalyzeR:GSCA]{GSCA}} or \code{\link[HTSanalyzeR:NWA]{NWA}}, this argument is an object of class
\code{\link[HTSanalyzeR:GSCA]{GSCA}} or \code{\link[HTSanalyzeR:NWA]{NWA}}.
}
	\item{...}{
other arguments depending on class (see below for the arguments supported by 
the method of class \code{\link[HTSanalyzeR:GSCA]{GSCA}} and \code{\link[HTSanalyzeR:NWA]{NWA}})
}
\describe{
	\item{what:}{
a single character value or a character vector of key words specifying
what to print (see details below). 
}
}

}
\details{
For an object of class \code{\link[HTSanalyzeR:GSCA]{GSCA}}, the key words are 'GSC' 
(the slot 'listOfGeneSetCollections'), 'GeneList' (the slot 'geneList'), 'Hits'
(the slot 'hits'), 'Para' (the slot 'para'),  'Result' (the slot 'result')
and 'ALL' (all slots). 

For an object of class \code{\link[HTSanalyzeR:GSCA]{NWA}}, the key words include 
'Pval' (the slot 'pvalue'), 'Phenotype' (the slot 'phenotype'), 'Interactome' 
(the slot 'interactome'), 'Para' (the slot 'fdr'), 'Result' (the slot 'result') 
and 'ALL' (all slots).
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
}
}











