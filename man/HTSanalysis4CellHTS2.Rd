\name{HTSanalysis4CellHTS2}
\alias{HTSanalysis4CellHTS2}
\title{
Analysis pipeline for cellHTS2 objects
}
\description{
This function writes an html report following a complete analysis of a dataset with the two wrapper functions geneSetAnalyze and networkAnalysis. 
}
\usage{
HTSanalysis4CellHTS2(normCellHTSobject, scoreSign="-", scoreMethod="zscore", summarizeMethod="mean", annotationColumn="GeneID", 
species=c("Drosophila_melanogaster"), initialIDs="FlybaseCG", 
duplicateRemoverMethod="max", orderAbsValue=FALSE, listOfGeneSetCollections, 
cutoffHitsEnrichment=2, pAdjMethod="BH", nPermutations=1000, 
minGeneSetSize=15, exponent=1, whichSetIsKEGGIds, whichSetIsGOIds, 
nwStatsControls="neg", nwStatsAlternative="two.sided", 
nwStatsTests=c("T-test"), nwStatsColumns=c("t.test.pvalues.two.samples",
"t.test.pvalues.one.sample"),nwAnalysisFdr=0.001, nwAnalysisGenetic=FALSE, 
biogridObject=NA, nwAnalysisOrder=2, nGseaPlots=10)
}
\arguments{
  \item{normCellHTSobject}{
A normalized, configured and annotated cellHTS object.
}
  \item{scoreSign}{
A sign argument for the scoring function from cellHTS2 (see 'scoreReplicates').
}
  \item{scoreMethod}{
A method argument for the scoring function from cellHTS2 (see 'scoreReplicates').
}
  \item{summarizeMethod}{
A summary argument for the summarization function from cellHTS2 (see 'summarizeReplicates').
}
  \item{annotationColumn}{
The name of the column in the fData(normCellHTSobject) data frame from which the feature identifiers will be extracted
}
  \item{species}{
The species for which the data should be read. One of the following: "Drosophila_melanogaster", "Homo_sapiens", "Rattus_norvegicus", "Mus_musculus", "Caenorhabditis_elegans".
}
  \item{initialIDs}{
The type of identifiers in the column specified by annotationColumn these have to be one of: "Ensembl.transcript", "Ensembl.prot", "Ensembl.gene", "Entrez.gene", "RefSeq", "Symbol", "GenBank", "Flybase", "FlybaseCG", "FlybaseProt"
}
  \item{duplicateRemoverMethod}{
A choice of method to remove the duplicates (should the minimum, maximum of average observation for a same construct be kept). The minimum and maximum should be understood in terms of absolute values (i.e. min/max effect, no matter the sign). The fold change average method converts the fold changes to ratios, averages them and converts the average back to a fold change.
}
  \item{orderAbsValue}{
A choice of whether the values for the overrepresentation analysis should be converted to absolute value and then ordered (if TRUE),or ordered as they are (if FALSE).
}
  \item{listOfGeneSetCollections}{
A list of gene set collections. Even if only one GeneSetCollection is being tested,it must be entered as an element of a 1-element list, e.g. 

\code{ListOfGeneSetCollections = list(YourOneGeneSetCollection)}. 
Naming the elements of ListOfGeneSetCollections will result in these names being associated with the relevant dataframes in the output (meaningful names are advised)
}
  \item{cutoffHitsEnrichment}{
The cutoff that is used in the definition of the hits for the hypergeometric tests in the overrepresentation analysis. This cutoff is used in absolute value, since it is applied on scores, i.e. a cutoff of 2 when using z-scores means that we are selecting values that are two standard deviations away from the median of all samples. Therefore, the cutoff should be a positive number.
}
  \item{pAdjMethod}{
A p value adjustment method for multiple hypothesis testing in the overrepresentation analysis
}
  \item{nPermutations}{
Number of permutations used in the computation of the GSEA statistics
}
  \item{minGeneSetSize}{
Minimum gene set size for a gene set to be considered in the overrepresentation analysis. 
}
  \item{exponent}{
Exponent that is used in the GSEA. 
}
 \item{whichSetIsKEGGIds}{
Index or vector of indices of the gene set collections in the gene list that are KEGG pathways (as created by the KeggGeneSets function).
}
  \item{whichSetIsGOIds}{
Index or vector of indices of the gene set collections in the gene list that are GO terms (as created by the GOGeneSets function).
}
 \item{nwStatsControls}{
Identifiers of the control samples in the cellHTS object, to be used as a control population in the 2 samples tests (this HAS to be corresponding to how these control wells have been annotated in the column "controlStatus" of the fData(normCellHTSobject) dataframe) if nothing is specified, the function will look for negative controls labelled "neg"
}
  \item{nwStatsAlternative}{
The alternative used for the statistical tests in the networkAnalysis function. two.sided/less/greater (only one of those 3): which is the alternative hypothesis to be tested.
}
  \item{nwStatsTests}{
The statistical tests performed in the networkAnalysis function. Any combination of "T-test","MannWhitney","RankProduct".
}
  \item{nwStatsColumns}{
Any (relevant, i.e. that is produced in the tests) combination of "t.test.pvalues.two.samples", "t.test.pvalues.one.sample", "mannW.test.pvalues.one.sample", "mannW.test.pvalues.two.samples", "rank.product.pfp.greater", "rank.product.pfp.less"
}
  \item{nwAnalysisFdr}{
The FDR used in the networkAnalysis function for the scores calculation. 
}
  \item{nwAnalysisGenetic}{
FALSE/TRUE indicating if the genetic interaction data from The Biogrid dataset were kept in the network analysis.
}
  \item{biogridObject}{
The Biogrid data set used in the networkAnalysis function. If this argument is set to NA (default), the function will download from the BioGRID the complete data set of interactions for the species specified by the species argument.
}
  \item{nwAnalysisOrder}{
The order used in the networkAnalysis function for the scores calculation.
}
  \item{nGseaPlots}{
The number of plots to be produced for the GSEA analysis.  For each gene set collection, plots are produced for the "nplots" most significant p values.  
}
}
\details{
This function does a gene set overrepresentation analysis (hypergeometric test and gsea), a search for enriched subnetworks and writes a report in a subdirectory called HTSAnalyzerReport.  The GSEA plots will be produced in the working directory and will be fetched, along as the network representation, in that directory.  If you move those items the function might not work properly. 

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
#Prepare the data
experimentName<-"KcViab"
dataPath<-system.file(experimentName, package="cellHTS2")
x<-readPlateList("Platelist.txt",name=experimentName,path=dataPath)
x<-configure(x,descripFile="Description.txt",confFile="Plateconf.txt",
logFile="Screenlog.txt",path=dataPath)
xn<-normalizePlates(x,scale="multiplicative",log=FALSE,method="median",
varianceAdjust="none")
xn<-annotate(xn,geneIDFile="GeneIDs_Dm_HFA_1.1.txt",path=dataPath)
#Make the list of gene sets
Dm.GO.MF<-GOGeneSets(species="Drosophila_melanogaster",ontologies=c("MF"))
Dm.GO.BP<-GOGeneSets(species="Drosophila_melanogaster",ontologies=c("BP"))
Dm.GO.CC<-GOGeneSets(species="Drosophila_melanogaster",ontologies=c("CC"))
kegg.droso<-KeggGeneSets(species="Drosophila_melanogaster");
gsc.list<-list(Dm.GO.MF=Dm.GO.MF,Dm.GO.BP=Dm.GO.BP,Dm.GO.CC=Dm.GO.CC,
kegg.droso=kegg.droso)
cellHTS2DrosoData<-xn
HTSanalysis4CellHTS2(
 	normCellHTSobject=cellHTS2DrosoData,
 	annotationColumn="GeneID",
 	species=c("Drosophila_melanogaster"),
 	initialIDs="FlybaseCG",
 	listOfGeneSetCollections=gsc.list,
 	whichSetIsKEGGIds=2,
 	whichSetIsGOIds=1,
	)
}
}
