\name{analyzeGeneSetCollections}
\alias{analyzeGeneSetCollections}
\title{
Hypergeometric tests and Gene Set Enrichment Analyses over a list of gene
set collections
}
\description{
This function takes a list of gene set collections, a named phenotype 
vector (with names of the phenotype vector as the universe), a vector 
of hits (gene names only) and returns the results of hypergeometric 
and gene set enrichment analyses for all of the gene set collections 
(with multiple hypothesis testing corrections).
}
\usage{
analyzeGeneSetCollections(listOfGeneSetCollections, geneList, hits, 
pAdjustMethod="BH", pValueCutoff=0.05, nPermutations=1000, 
minGeneSetSize=15, exponent=1, verbose=TRUE)
}
\arguments{
  \item{listOfGeneSetCollections}{
a list of gene set collections (a 'gene set collection' is a list of 
gene sets). Even if only one collection is being tested, it must be 
entered as an element of a 1-element list, e.g. 
\code{ListOfGeneSetCollections = list(YourOneGeneSetCollection)}. 
Naming the elements of listOfGeneSetCollections will result in these 
names being associated with the relevant dataframes in the output 
(meaningful names are advised)
}
  \item{geneList}{
a numeric or integer vector of phenotypes in descending or ascending 
order with elements named by their EntrezIds (no duplicates nor NA values)
}
  \item{hits}{
a character vector of the EntrezIds of hits, as determined by the user
}
  \item{pAdjustMethod}{
a single character value specifying the p-value adjustment method to be 
used (see 'p.adjust' for details)
}
  \item{pValueCutoff}{
a single numeric value specifying the cutoff for p-values considered 
significant
}
  \item{nPermutations}{
a single integer or numeric value specifying the number of permutations 
for deriving p-values in GSEA
}
  \item{minGeneSetSize}{
a single integer or numeric value specifying the minimum number of 
elements in a gene set that must map to elements of the gene universe. 
Gene sets with fewer than this number are removed from both hypergeometric 
analysis and GSEA.
}
  \item{exponent}{
a single integer or numeric value used in weighting phenotypes in GSEA 
(see the function \code{\link[HTSanalyzeR:gseaScores]{gseaScores}})
}
  \item{verbose}{
a single logical value specifying to display detailed messages (when 
verbose=TRUE) or not (when verbose=FALSE)
}
}
\details{
 All gene names must be EntrezIds in 'listOfGeneSetCollections', 
 'geneList', and 'hits'.
}
\value{
\item{HyperGeo.results}{a list of dataframes containing the results for 
all gene set collections in the input.}
\item{GSEA.results}{a similar list of dataframes containing the results 
from GSEA. As an example, to access the GSEA results for a gene set 
collection named "MyGeneSetCollection", one would enter:
output$GSEA.results$MyGeneSetCollection}
\item{Sig.pvals.in.both}{a list of dataframes containing the gene sets 
with p-values considered significant in both hypergeometric test and 
GSEA, before p-value correction. Each element of the list contains the 
results for one gene set collection.}
\item{Sig.adj.pvals.in.both}{a list of dataframes containing the gene 
sets with p-values considered significant in both hypergeometric test 
and GSEA, after p-value correction. Each element of the list contains 
the results for one gene set collection.}
}
\author{
John C. Rose, Xin Wang
}
\seealso{
\code{\link[HTSanalyzeR:analyze]{analyze}}
}
\examples{
\dontrun{
library(org.Dm.eg.db)
library(GO.db)
library(KEGG.db)
##load phenotype vector (see the vignette for details about the 
##preprocessing of this data set)
data("KcViab_Data4Enrich")
##Create a list of gene set collections for Drosophila melanogaster (Dm)
GO_MF <- GOGeneSets(species="Dm", ontologies="MF")
PW_KEGG <- KeggGeneSets(species="Dm")
ListGSC <- list(GO_MF=GO_MF, PW_KEGG=PW_KEGG)
##Conduct enrichment analyses
GSCAResults <- analyzeGeneSetCollections(
		listOfGeneSetCollections=ListGSC,
		geneList=KcViab_Data4Enrich,
		hits=names(KcViab_Data4Enrich)[which(abs(KcViab_Data4Enrich)>2)],
		pAdjustMethod="BH",
		nPermutations=1000,
		minGeneSetSize=200,
		exponent=1,
		verbose=TRUE
)
}
}
