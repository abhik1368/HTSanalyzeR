\name{appendGSTerms}
\alias{appendGSTerms}
\alias{appendGSTerms,GSCA-method}
\title{
Append gene set terms to GSCA results
}
\description{
This is a generic function.

When implemented as the S4 method for objects of class \code{\link[HTSanalyzeR:GSCA-class]{GSCA}},
this function finds corresponding annotation terms for KEGG and GO gene sets 
and inserts a column named "Gene.Set.Term" to each data frame in the GSCA results. 

To use this function for objects of class \code{\link[HTSanalyzeR:GSCA-class]{GSCA}}:

appendGSTerms(object, keggGSCs=NULL, goGSCs=NULL)
}
\usage{
appendGSTerms(object, ...)
}
\arguments{
	\item{object}{
an object. When this function is implemented as the S4 method of class 'GSCA', 
this argument is an object of class 'GSCA'.
}
	\item{...}{
other arguments depending on class (see below for the arguments supported by 
the method of class 'GSCA')
}
\describe{
  \item{keggGSCs:}{
a character vector of names of all KEGG gene set collections
}
  \item{goGSCs:}{
a character vector of names of all GO gene set collections
}
}

}
\details{
This function makes the GSCA results more readable by appending a column of 
terms for KEGG and GO gene sets. To do this, the user needs to specify the names 
of the gene set collections based on KEGG and GO, respectively. 

For each KEGG gene set, the species code in the KEGG id will be trimmed off, and 
then mapped to its corresponding annotation term using the function \code{mget} of 
the package \code{AnnotationDbi}.

For each GO gene set, the GO id will be mapped to corresponding GO term by 
the function \code{Term} of the package \code{GO.db}. 
}
\value{
In the end, this function will return an updated object of class \code{\link[HTSanalyzeR:GSCA-class]{GSCA}}.  
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
gsca <- new("GSCA", listOfGeneSetCollections=gscList,
geneList = KcViab_Data4Enrich, hits = hits)
##print gsca
gsca
##do preprocessing (KcViab_Data4Enrich has already been preprocessed)
gsca <- preprocess(gsca, species="Dm", initialIDs = "Entrez.gene", 
keepMultipleMappings = TRUE, duplicateRemoverMethod = "max", 
orderAbsValue = FALSE)
##print gsca again
gsca
##do hypergeometric tests and GSEA
gsca <- analyze(gsca, para = list(pValueCutoff = 0.05, pAdjustMethod 
= "BH", nPermutations = 1000, minGeneSetSize = 100,exponent = 1))
##append Kegg and GO gene set terms
gsca<-appendGSTerms(gsca, keggGSCs="PW_KEGG")
}
}










