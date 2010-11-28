\name{summary}
\alias{summary}
\title{
Print summary information about an object of class \code{GSCA} or \code{NWA}
}
\description{
This is an S4 method, which can be applied to an object of 
either class \code{GSCA} or \code{NWA}. It prints a summary of slots 
specified by users. 
}
\usage{
summary(object, what="ALL")
}
\arguments{
	\item{object}{
an object of class 'GSCA' or 'NWA'
}
	\item{what}{
a single character value or a character vector of key words specifying what to 
print. 
}
}
\details{
For an object of class 'GSCA', key words are 'GSC' (slot 'listOfGeneSetCollections'), 
'GeneList' (slot 'geneList'), 'Hits' (slot 'hits'), 'Para' (slot 'para'),  
'Result' (slot 'result') and 'ALL' (all slots). 

For an object of class 'NWA', key words include 'Pval' (slot 'pvalue'), 
'Phenotype' (slot 'phenotype'), 'Interactome' (slot 'interactome'), 'Para' 
(slot 'fdr'), 'Result' (slot 'result') and 'ALL' (all slots).
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
PW.KEGG <- KeggGeneSets(species = "Dm")
gsc.list <- list(PW.KEGG = PW.KEGG)
##create an object of class 'GSCA'
gsca <- new("GSCA", listOfGeneSetCollections=gsc.list,
geneList = data4enrich, hits = hits)
##print summary of gsca
summary(gsca)
##do preprocessing (data4enrich has already been preprocessed)
gsca <- preprocess(gsca, species="Dm", initialIDs = "Entrez.gene", 
keepMultipleMappings = TRUE, duplicateRemoverMethod = "max", 
orderAbsValue = FALSE)
##print summary of gsca again
summary(gsca)
##do hypergeometric tests and GSEA
gsca <- analyze(gsca, para = list(pValueCutoff = 0.05, pAdjustMethod 
= "BH", nPermutations = 1000, minGeneSetSize = 100,exponent = 1))
##print summary of results
summary(gsca, what="Result")
}
}











