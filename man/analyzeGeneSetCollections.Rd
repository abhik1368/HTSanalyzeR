\name{analyzeGeneSetCollections}
\alias{analyzeGeneSetCollections}
\title{
Enrichment analysis by hypergeometric test and GSEA
}
\description{
This function takes a list of GeneSetCollections, a named phenotype vector (with names(phenotype vector)=GeneUniverse), a list of hits (names only) and returns the results of hypergeometric and gene set enrichment analysis for all of the gene set collections (with multiple hypothesis testing correction).
}
\usage{
analyzeGeneSetCollections(ListOfGeneSetCollections, GeneList, hits, 
pAdjustMethod = "BH", p.value.cutoff = 0.05, npermutations = 1000, 
min.gene.set.size = 15, exponent = 1, whichSetIsKEGGIds = "none", 
whichSetIsGOIds = "none", verbose = TRUE)
}
\arguments{
  \item{ListOfGeneSetCollections}{
A list of gene set collections. Even if only one collection is being tested,it must be entered as an element of a 1-element list, e.g. 

\code{ListOfGeneSetCollections = list(YourOneGeneSetCollection)}. 
Naming the elements of ListOfGeneSetCollections will result in these names being associated with the relevant dataframes in the output (meaningful names are advised)
}
  \item{GeneList}{
A vector of numeric phenotypes in descending or ascending order with elements named by their EntrezIds (no duplicates nor NA values)
}
  \item{hits}{
A vector of the EntrezIds of hits, as determined by the user
}
  \item{pAdjustMethod}{
The p-value adjustment method to be used.
}
  \item{p.value.cutoff}{
The cutoff for p-values considered significant.
}
  \item{npermutations}{
The number of permutations for deriving pvalues in GSEA
}
  \item{min.gene.set.size}{
The minimum number of elements in a gene set that must map to elements of the gene universe. Gene sets with fewer than this number are removed from both hypergeometric analysis and GSEA.
}
  \item{exponent}{
Used in weighting phenotypes in GSEA (see "gseascores" function)
}
  \item{whichSetIsKEGGIds}{
 The index or a vector of indices of the GeneSetCollections with KEGGIds as setNames. The KEGG pathway names will be appended to these Ids in the output.
}
  \item{whichSetIsGOIds}{
The index or a vector of indices of the GeneSetCollections with GO Ids as setNames. The GO term names will be appended to these Ids in the output.
}
  \item{verbose}{
A logical value to decide to display detailed messages (when verbose=TRUE) or not (when verbose=FALSE)
}
}
\details{
 All gene names must be EntrezIds in the GeneSetCollections, GeneList, and the list of hits.
}
\value{
\item{HyperGeo.results}{a list of dataframes containing the results for all gene set collections in the input,and the last element of this list contains all of these results in a single dataframe.}
\item{GSEA.results}{a similar list of dataframes containing the results from GSEA. As an example, to access the GSEA results for a gene set collection named "MyGeneSetCollection", one would enter:output$GSEA.results$MyGeneSetCollection}
\item{Sig.pvals.in.both}{a list of dataframes containing the gene sets with p-values considered significant in both hypergeometric test and GSEA, before pvalue correction.  Each element of the list contains the results for one gene set collection, and the last element of this list contains all of these results in a single dataframe.}
\item{Sig.adj.pvals.in.both}{a list of dataframes containing the gene sets with p-values considered significant in both hypergeometric test and GSEA, after pvalue correction.  Each element of the list contains the results for one gene set collection, and the last element of this list contains all of these results in a single dataframe.}
}
\author{
John C. Rose
}
\examples{
\dontrun{
library(org.Dm.eg.db)
library(GO.db)
library(KEGG.db)
library(cellHTS2)
#Prepare the data
experimentName<-"KcViab"
dataPath<-system.file(experimentName, package="cellHTS2")
x<-readPlateList("Platelist.txt",name=experimentName,path=dataPath)
x<-configure(x,descripFile="Description.txt",confFile="Plateconf.txt",
logFile="Screenlog.txt",path=dataPath)
xn<-normalizePlates(x,scale="multiplicative",log=FALSE,method="median",
varianceAdjust="none")
xn<-annotate(xn,geneIDFile="GeneIDs_Dm_HFA_1.1.txt",path=dataPath)
xsc<-scoreReplicates(xn,sign="-",method="zscore")
xsc<-summarizeReplicates(xsc,summary="mean")
data4enrich<-as.vector(Data(xsc));
names(data4enrich)<-fData(xsc)[,"GeneID"] 
data4enrich<-data4enrich[which(!is.na(names(data4enrich)))]
data4enrichentrez<-drosoAnnotationConvertor(geneList=data4enrich,initialIDs="FlybaseCG")
data4enrichentrez<-data4enrichentrez[which(!is.na(data4enrichentrez))]
data4enrichentrez<-duplicateRemover(geneList=data4enrichentrez,method="max")
#Create the list of gene set collections for Drosophila melanogaster (Dm)
Dm.GO.MF<-GOGeneSets(species="Dm",ontologies=c("MF"))
Dm.GO.BP<-GOGeneSets(species="Dm",ontologies=c("BP"))
Dm.GO.CC<-GOGeneSets(species="Dm",ontologies=c("CC"))
kegg.droso<-KeggGeneSets(species="Dm");
gsc.list<-list(Dm.GO.MF=Dm.GO.MF,Dm.GO.BP=Dm.GO.BP,Dm.GO.CC=Dm.GO.CC,kegg.droso=kegg.droso)
#Conduct enrichment analyses
enrichment.analysis<-analyzeGeneSetCollections(
		ListOfGeneSetCollections=gsc.list,
		GeneList=data4enrichentrez,
		hits=names(data4enrichentrez)[which(abs(data4enrichentrez)>2)],
		pAdjustMethod="BH",
		npermutations=1000,
		min.gene.set.size=15,
		exponent=1,
		whichSetIsKEGGIds=4,
		whichSetIsGOIds=c(1,2,3),
		verbose=TRUE
	)
}
}
