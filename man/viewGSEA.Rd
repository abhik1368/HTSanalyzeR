\name{viewGSEA}
\alias{viewGSEA}
\title{
View the enrichment analysis result of a gene set
}
\description{
This is an S4 method for objects of class 'GSCA' to view a figure
including the ranked phenotypes, positions of the gene set and the
location of the maximum enrichment score.  
}
\usage{
viewGSEA(object, gsc, gs) 
}
\arguments{
	\item{object}{
an object of class 'GSCA'
}
	\item{gscName}{
a single character value specifying the name of the gene set collection 
where the gene set is
}
	\item{gsName}{
a single character value specifying the name of the gene set to be plot
}
}
\details{
We suggest to 'getTopGeneSets' first to have an idea of the names of top 
significant gene sets, and then choose to plot.
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
##print top significant gene sets in GO.BP
topPWKEGG<-getTopGeneSets(gsca, "GSEA.results", "PW_KEGG", allSig=TRUE)
##view a gene set
viewGSEA(gsca, "PW_KEGG", topPWKEGG[["PW_KEGG"]][1])
}
}











