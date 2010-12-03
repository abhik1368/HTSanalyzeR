\name{summarize}
\alias{summarize}
\alias{summarize,GSCA-method}
\alias{summarize,NWA-method}
\alias{show}
\alias{show,GSCA-method}
\alias{show,NWA-method}
\title{
Print summary information about an object of class GSCA or NWA
}
\description{
This is a generic function.

When implemented as the S4 method for objects of class 'GSCA' or 'NWA', 
this function prints a summary of slots specified by users. 

To use this function for objects of class \code{GSCA} or \code{NWA}:

summarize(object, what="ALL")
}
\usage{
summarize(object, ...)
}
\arguments{
	\item{object}{
an object. When this function is implemented as the S4 method of class 'GSCA' 
or 'NWA', this argument is an object of class 'GSCA' or 'NWA'.
}
	\item{...}{
other arguments depending on class (see below for the arguments supported by 
class 'GSCA' and 'NWA')
}
\describe{
	\item{what:}{
a single character value or a character vector of key words specifying
what to print. 
}
}

}
\details{
For an object of class 'GSCA', key words are 'GSC' (slot
'listOfGeneSetCollections'), 'GeneList' (slot 'geneList'), 'Hits'
(slot 'hits'), 'Para' (slot 'para'),  'Result' (slot 'result') and 'ALL'
(all slots). 

For an object of class 'NWA', key words include 'Pval' (slot 'pvalue'), 
'Phenotype' (slot 'phenotype'), 'Interactome' (slot 'interactome'),
'Para' (slot 'fdr'), 'Result' (slot 'result') and 'ALL' (all slots).
}

\author{
Xin Wang \email{xw264@cam.ac.uk}
}

\examples{
\dontrun{
library(org.Dm.eg.db)
library(KEGG.db)
##load data for enrichment analyses
data("data4enrich")
##select hits
hits <- names(data4enrich)[which(abs(data4enrich) > 2)]
##set up a list of gene set collections
PW_KEGG <- KeggGeneSets(species = "Dm")
gscList <- list(PW_KEGG = PW_KEGG)
##create an object of class 'GSCA'
gsca <- new("GSCA", listOfGeneSetCollections=gscList, geneList =
data4enrich, hits = hits)
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
= "BH", nPermutations = 1000, minGeneSetSize = 100, exponent = 1))
##print summary of results
summarize(gsca, what="Result")
}
}











