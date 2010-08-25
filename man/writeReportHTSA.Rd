\name{writeReportHTSA}
\alias{writeReportHTSA}
\title{
Function to write the report of an analysis from HTSAnalyzer
}
\description{
This function writes an html report following a complete analysis of a dataset with the two wrapper functions geneSetAnalyze and networkAnalysis. Most of the parameters of this function can simply be copied and paste from the parameters of the two above functions.
}
\usage{
writeReportHTSA(experimentName, enrichmentAnalysis, cutoffHits, hits, 
listOfGeneSetCollections, geneList, geneListName, 
p.adj.method = "BH", nPermutations = 1000, min.gene.set.size = 15, 
exponent = 1, nwAnalysisOutput, nwAnalysisGraphFile = "EnrichedSubNw.png", 
controls = "neg", alternative = "two.sided", tests = c("T-test"), 
columns = c("t.test.pvalues.two.samples", "t.test.pvalues.one.sample"), 
species = c("Drosophila_melanogaster"), fdr = 0.001, genetic = FALSE, 
networkObject = NA, nGseaPlots = 10, whichSetIsKEGGIds=NA, whichSetIsGOIds=NA)
}
\arguments{
  \item{experimentName}{
The name of the experiment.  This is just for you own record.
}
  \item{enrichmentAnalysis}{
An object resulting from enrichment analysis with the function geneSetAnalyze.
}
  \item{cutoffHits}{
The cutoff that was used in the definition of the hits for the hypergeometric tests in the overrepresentation analysis.  This is only for you own record.
}
  \item{hits}{
The hits for the hypergeometric tests in the overrepresentation analysis (as a vector of gene names, expected to be entrez gene identifiers). 
}
  \item{listOfGeneSetCollections}{
The list of gene set collections that was used as input in the geneSetAnalyze function
}
  \item{geneList}{
The gene list that was used as input in the geneSetAnalyze function
}
  \item{geneListName}{
The name of the gene list that was used as input in the geneSetAnalyze function.  This is only for your own record.
}
  \item{p.adj.method}{
p value adjustment method for the overrepresentation analysis by geneSetAnalyze.   This is only for your own record.
}
  \item{nPermutations}{
Number of permutations used in the computation of the GSEA statistics
}
  \item{min.gene.set.size}{
Minimum gene set size for a gene set to be considered in the overrepresentation analysis.  This is just for your own record.
}
  \item{exponent}{
Exponent that was used in the GSEA. 
}
  \item{nwAnalysisOutput}{
The list that is obtained as a result of the networkAnalysis function
}
  \item{nwAnalysisGraphFile}{
The name of the image file containing a representation of the enriched subnetwork, produced by the networkAnalysis function.  If the name of the graph produced by the networkAnalysis function was not changed, this function will find the graph.
}
  \item{controls}{
Identifiers of the control samples in the cellHTS object, as input in the networkAnalysis function. This is for your own record.
}
  \item{alternative}{
The alternative used for the statistical tests in the networkAnalysis function. This is for your own record.
}
  \item{tests}{
The statistical tests performed in the networkAnalysis function. This is for your own record.
}
  \item{columns}{
The columns containing the statistical tests used in the networkAnalysis function. This is for your own record.
}
  \item{species}{
The species, as in the networkAnalysis function.
}
  \item{fdr}{
The FDR used in the networkAnalysis function for the scores calculation. This is for your own record.
}
  \item{genetic}{
FALSE/TRUE indicating if the genetic interaction data from The Biogrid dataset were kept in the network analysis.
}
  \item{networkObject}{
The interaction data set used in the networkAnalysis function. This is for your own record.
}
  \item{nGseaPlots}{
The number of plots to be produced for the GSEA analysis.  For each gene set collection, plots are produced for the "nplots" most significant p values.  
}
  \item{whichSetIsKEGGIds}{
Index or vector of indices of the gene set collections in the gene list that are KEGG pathways (as created by the KeggGeneSets function).
}
  \item{whichSetIsGOIds}{
Index or vector of indices of the gene set collections in the gene list that are GO terms (as created by the GOGeneSets function).
}
}
\details{
This function takes the outputs of the two wrapper functions either in the R environment or in the working directory, and writes a report in a subdirectory called HTSAnalyzerReport.  The GSEA plots will be produced in the working directory and will be fetched, along as the network representation, in that directory.  If you move those items the function might not work properly. Some parameters are described as "only for your own record", they are asked as input to thsi function so that immediately after analysis you can keep a trace directly in the report of the settings that you used for the analysis.

The following parameters can be directly copied from the input of "geneSetAnalyze": ListOfGeneSetCollections, GeneList, p.adj.method , npermutations , min.gene.set.size , exponent , whichSetIsKEGGIds, whichSetIsGOIds.

The following parameters can be directly copied from the input of "networkAnalysis":controls, alternative , tests , columns , species , fdr, genetic , biogridobject.

The following parameters are for your record only: Genelistname, experimentName, cutoffHits.

The following objects are generated by the two wrapper function "geneSetAnalyze" and "networkAnalysis": enrichmentAnalysis, nwAnalysisoutput, nwAnalysisgraphFile.

For the KEGG and GO gene sets, this function produces links to the relevant database entries.  These links will only work if the gene set names are in the following format: 
\itemize{
\item for KEGG 3 letters (species ID) followed by a set of numbers and then any number of non-digit characters
\item for GO "GO:" immediately followed by the GO term number (digits only) and any number of non-digits characters.
}
If the gene set names are not in these format, simply do not tag the gene set collections as "Kegg" or "GO".  This way, the GSEA plots will be produced normally but the links to web databases will not be produced.
}
\value{
Produce a set of html pages.  Start from the page called "index.html" to navigate those pages.
}
\author{
Camille Terfve
}
\examples{
\dontrun{
library(org.Dm.eg.db)
library(GO.db)
library(KEGG.db)
library(cellHTS2)
#Prepare the data for networks analysis
experimentName<-"KcViab"
dataPath<-system.file(experimentName, package="cellHTS2")
x<-readPlateList("Platelist.txt",name=experimentName,path=dataPath)
x<-configure(x,descripFile="Description.txt",confFile="Plateconf.txt",
logFile="Screenlog.txt",path=dataPath)
xn<-normalizePlates(x,scale="multiplicative",log=FALSE,method="median",
varianceAdjust="none")
xn<-annotate(xn,geneIDFile="GeneIDs_Dm_HFA_1.1.txt",path=dataPath)
#Prepare the data for gene set over-representation analysis
xsc<-scoreReplicates(xn,sign="-",method="zscore")
xsc<-summarizeReplicates(xsc,summary="mean")
data4enrich<-as.vector(Data(xsc));
names(data4enrich)<-fData(xsc)[,"GeneID"] 
data4enrich<-data4enrich[-which(is.na(names(data4enrich))==TRUE)]
data4enrichentrez<-drosoAnnotationConvertor(geneList=data4enrich,
initialIDs="FlybaseCG")
data4enrichentrez<-data4enrichentrez[-which(is.na(data4enrichentrez))]
data4enrichentrez<-duplicateRemover(geneList=data4enrichentrez,method="max")
#Make the list of gene sets
Dm.GO.MF<-GOGeneSets(species="Drosophila_melanogaster",ontologies=c("MF"))
Dm.GO.BP<-GOGeneSets(species="Drosophila_melanogaster",ontologies=c("BP"))
Dm.GO.CC<-GOGeneSets(species="Drosophila_melanogaster",ontologies=c("CC"))
kegg.droso<-KeggGeneSets(species="Drosophila_melanogaster");
gsc.list<-list(Dm.GO.MF=Dm.GO.MF,Dm.GO.BP=Dm.GO.BP,Dm.GO.CC=Dm.GO.CC,
kegg.droso=kegg.droso)
#Run the enrichment analysis
enrichment.analysis<-analyzeGeneSetCollections(
	ListOfGeneSetCollections=gsc.list,
	GeneList=data4enrichentrez,
	hits=names(data4enrichentrez)[which(abs(data4enrichentrez)>2)],
	pAdjMethod="BH",
	npermutations=1000,
	min.gene.set.size=15,
	exponent=1,
	whichSetIsKEGGIds=4,
	whichSetIsGOIds=c(1,2,3)
	)
#Run the networks analysis
test.module<-networkAnalysis(
	cellHTSobject=xn,
	annotationColumn="GeneID",
	controls="neg",
	alternative="two.sided",
	logged=FALSE,
	tests=c("T-test"),
	columns=c("t.test.pvalues.two.samples", "t.test.pvalues.one.sample"),
	species=c("Drosophila_melanogaster"),
	initialIDs="FlybaseCG",
	fdr=0.001,
	genetic=FALSE,
	biogridObject=NA,phenotypeVector=data4enrichentrez)	
#Write the report
writeReportHTSA(
	experimentName="cellHTS2DrosoData",
	enrichmentAnalysis=enrichment.analysis,
	cutoffHits=2,
	hits=names(data4enrichentrez)[which(abs(data4enrichentrez)>2)],
	listOfGeneSetCollections=gsc.list,
	geneList=data4enrichentrez,
	p.adj.method="BH",
	nPermutations=1000,
	min.gene.set.size=15,
	exponent=1,
	nwAnalysisOutput=test.module,
	nwAnalysisGraphFile="EnrichedSubNw.png",
	controls="neg",
	alternative="two.sided",
	tests=c("T-test"),
	columns=c("t.test.pvalues.two.samples","t.test.pvalues.one.sample"),
	species=c("Drosophila_melanogaster"),
	fdr=0.001,
	genetic=FALSE,
	networkObject=NA,
	nGseaPlots=10,
	geneListName="data4enrichentrez",
	whichSetIsKEGGIds=4,
	whichSetIsGOIds=c(1,2,3)
	)
}
}
