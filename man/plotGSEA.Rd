\name{plotGSEA}
\alias{plotGSEA}
\alias{plotGSEA,GSCA-method}
\title{
Plot and save figures of GSEA results of top significant gene sets
}
\description{
This is a generic function.

When implemented as the S4 method for objects of class 'GSCA', this function 
plots figures of the ranked phenotypes, positions of the gene set and the 
location of the maximum enrichment score for top significant gene sets.  

To use this function for objects of class 'GSCA':

plotGSEA(object, gscs, ntop=NULL, allSig=FALSE, filepath=".", output=
"png", ...)
}
\usage{
plotGSEA(object, ...)
}
\arguments{
	\item{object}{
an object. When this function is implemented as the S4 method of class 'GSCA', 
this argument is an object of class 'GSCA'.
}
	\item{...}{
other arguments. (see below for the arguments supported by class 'GSCA')
}
\describe{
	\item{gscs:}{
a character vector specifying the names of gene set collections of which
the top significant gene sets will be plot
}
	\item{ntop:}{
a single integer or numeric value specifying to plot how many gene sets
of top significance.  
}
	\item{allSig:}{
a single logical value. If 'TRUE', all significant gene sets (GSEA
adjusted p-value < 'pValueCutoff' of slot 'para') will be plot; otherwise,
only top 'ntop' gene sets will be plot.
}
	\item{filepath:}{
a single character value specifying where to store GSEA figures.
}
	\item{output:}{
a single character value specifying the format of output image: "pdf" or
"png"
}
	\item{...}{
other arguments used by function png() or pdf() such as width and height
}
}

}
\details{
To use this function to make GSEA plots of top significance, users can
only choose one method, either assign 'ntop' a integer or set 'allSig =
TRUE'. Exceptions will occur if both methods are used, or no method is
used. Please also note that 'ntop' is a cutoff for all gene set collections
in 'gscs'.

We suggest to summarize(gsca, what="Result") first to have an idea of how
many significant gene sets there are, and then choose to plot them all
or just top few ones.
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
gsca <- new("GSCA", listOfGeneSetCollections=gscList,
geneList = data4enrich, hits = hits)
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
##plot all significant gene sets
plotGSEA(gsca, gscs=c("PW_KEGG"), allSig=TRUE, filepath=".", output=
"pdf", width=8, height=8)
}
}











