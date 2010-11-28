\name{plotGSEA}
\alias{plotGSEA}
\title{
Plot figures for Gene Set Enrichment Analysis
}
\description{
This is an S4 method for objects of class 'GSCA' to plot GSEA results.  
}
\usage{
plotGSEA(object, gscs, ntop=NULL, allSig=FALSE, filepath=".")
}
\arguments{
	\item{object}{
an object of class 'GSCA'
}
	\item{gscs}{
a list of gene set collections. See argument 'collectionOfGeneSets' in function 
'collectionGsea' for more details about a gene set collection.
}
	\item{ntop}{
a single integer or numeric value specifying to plot how many gene sets of 
top significance.  
}
	\item{allSig}{
a single logical value. If 'TRUE', all significant gene sets (GSEA adjusted 
p-value < 'pValueCutoff' of slot 'para') will be plot; otherwise, only top 
'ntop' gene sets will be plot.
}
	\item{filepath}{
a single character value specifying where to store GSEA figures.
}
}
\details{
To use this function to make GSEA plots of top significance, users can only 
choose one method, either assign 'ntop' a integer or set 'allSig = TRUE'. 
Exceptions will occur if both methods are used, or no method is used. Please 
also note that 'ntop' is a cutoff for all gene set collections in 'gscs'.

We suggest to summary(gsca, what="Result") first to have an idea of how many 
significant gene sets there are, and then choose to plot them all or just 
top few ones.
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
##plot all significant gene sets
plotGSEA(gsca, gscs=c("PW.KEGG"), allSig=TRUE)

}
}











