\name{getTopGeneSets}
\alias{getTopGeneSets}
\title{
Select top significant gene sets from GSEA results
}
\description{
This is an S4 method for objects of class 'GSCA' to select top significant
gene sets from GSEA results for user-specified gene collections. If 'ntop'
is given, then top 'ntop' gene sets in gene set collections 'gscs' will be
selected and their names will be returned. If 'allSig=TRUE', then all
significant (adjusted p-value < 'pValueCutoff' see 'help("analyze")') will
be selected and their names will be returned.
}
\usage{
getTopGeneSets(object, resultName, gscs, ntop=NULL, allSig=FALSE)
}
\arguments{
	\item{object}{
an object of class 'GSCA'
}
	\item{reusltName}{
a single character value: 'HyperGeo.results' or 'GSEA.results'
} 
	\item{gscs}{
a character vector specifying the names of gene set collections of which
the top significant gene sets will be plot
}
	\item{ntop}{
a single integer or numeric value specifying to plot how many gene sets
of top significance.  
}
	\item{allSig}{
a single logical value. If 'TRUE', all significant gene sets (GSEA
adjusted p-value < 'pValueCutoff' of slot 'para') will be selected;
otherwise, only top 'ntop' gene sets will be selected.
}

}
\value{
a list including a character vector of names of top significant gene set 
for each gene set collection
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
}
}











