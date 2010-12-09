\name{HTSanalyzeR4cellHTS2}
\alias{HTSanalyzeR4cellHTS2}
\title{
An analysis pipeline for cellHTS2 objects
}
\description{
This function writes an html report following a complete analyses of a
dataset based on the two classes \code{\link[HTSanalyzeR:GSCA]{GSCA}} (Gene Set Collection Analysis)
and \code{\link[HTSanalyzeR:NWA]{NWA}} (NetWork Analysis) of this package. 
}
\usage{
HTSanalyzeR4cellHTS2(		
	normCellHTSobject,
	scoreSign = "-",
	scoreMethod = "zscore",
	summarizeMethod = "mean",
	annotationColumn = "GeneID",
	species = "Dm",
	initialIDs = "FlybaseCG",
	duplicateRemoverMethod = "max",
	orderAbsValue = FALSE,
	listOfGeneSetCollections,
	cutoffHitsEnrichment = 2,
	pValueCutoff = 0.05,
	pAdjustMethod = "BH",
	nPermutations = 1000,
	minGeneSetSize = 15,
	exponent = 1,
	keggGSCs,
	goGSCs,
	nwStatsControls = "neg",
	nwStatsAlternative = "two.sided",
	nwStatsTests = "T-test",
	nwStatsColumns = c("t.test.pvalues.two.samples", "t.test.pvalues.one.sample"),
	nwAnalysisFdr = 0.001,
	nwAnalysisGenetic = FALSE,
	interactionMatrix = NULL,
	nwAnalysisOrder = 2,
	ntop = NULL,
	allSig = TRUE,
	reportDir = "HTSanalyzerReport",
	verbose = TRUE
)
}
\arguments{
  \item{normCellHTSobject}{
a normalized, configured and annotated cellHTS object
}
  \item{scoreSign}{
a single character value specifying the 'sign' argument for the scoring
function from cellHTS2 (see 'scoreReplicates')
}
  \item{scoreMethod}{
a single character value specifying the 'method' argument for the scoring
function from cellHTS2 (see 'scoreReplicates')
}
  \item{summarizeMethod}{
a summary argument for the summarization function from cellHTS2 (see 
'summarizeReplicates')
}
  \item{annotationColumn}{
a single character value specifying the name of the column in the fData 
(cellHTSobject) dataframe from which the feature identifiers will be
extracted
}
	\item{species}{
a single character value specifying the species for which the data should
be read. The current version supports one of the following species: "Dm"
("Drosophila_melanogaster"), "Hs" ("Homo_sapiens"), "Rn" ("Rattus_norvegicus"),
"Mm" ("Mus_musculus"), "Ce" ("Caenorhabditis_elegans"). 
}
  \item{initialIDs}{
a single character value specifying the type of initial identifiers for
input 'geneList'. Current version can take one of the following types: 
"Ensembl.transcript", "Ensembl.prot", "Ensembl.gene", "Entrez.gene",
"RefSeq", "Symbol" and "GenBank" for all supported species; "Flybase",
"FlybaseCG" and "FlybaseProt" in addition for Drosophila Melanogaster; 
"wormbase" in addition for Caenorhabditis Elegans. 
}
   \item{duplicateRemoverMethod}{
a single character value specifying the method to remove the duplicates 
(should the minimum, maximum or average observation for a same construct
be kept). 
The current version provides "min" (minimum), "max" (maximum), "average" and
"fc.avg" (fold change average). The minimum and maximum should be understood
in terms of absolute values (i.e. min/max effect, no matter the sign). The
fold change average method converts the fold changes to ratios, averages
them and converts the average back to a fold change.
}
	\item{orderAbsValue}{
a single logical value determining whether the values should be converted 
to absolute value and then ordered (if TRUE), or ordered as they are (if
FALSE).	
}
  \item{listOfGeneSetCollections}{
a list of gene set collections (a 'gene set collection' is a list of gene
sets). Even if only one collection is being tested,it must be entered as
an element of a 1-element list, e.g. \code{ListOfGeneSetCollections =
list(YourOneGeneSetCollection)}. Naming the elements of listOfGeneSetCollections
will result in these names being associated with the relevant dataframes
in the output (meaningful names are advised)
}
  \item{cutoffHitsEnrichment}{
a single numeric or integer value specifying the cutoff that is used in
the definition of the hits for the hypergeometric tests in the over-
representation analysis. This cutoff is used in absolute value, since it
is applied on scores, i.e. a cutoff of 2 when using z-scores means that
we are selecting values that are two standard deviations away from the
median of all samples. Therefore, the cutoff should be a positive number.
}
  \item{pValueCutoff}{
a single numeric value specifying the cutoff for p-values considered
significant
}
  \item{pAdjustMethod}{
a single character value specifying the p-value adjustment method to be
used (see 'p.adjust' for details)
}
  \item{nPermutations}{
a single integer or numeric value specifying the number of permutations
for deriving p-values in GSEA
}
  \item{minGeneSetSize}{
a single integer or numeric value specifying the minimum number of elements
in a gene set that must map to elements of the gene universe. Gene sets
with fewer than this number are removed from both hypergeometric analysis
and GSEA.
}
  \item{exponent}{
a single integer or numeric value used in weighting phenotypes in GSEA
(see "gseaScores" function)
}
  \item{keggGSCs}{
a character vector of names of all KEGG gene set collections. This will
help create web links for KEGG terms. 
}
  \item{goGSCs}{
a character vector of names of all GO gene set collections. This will help 
create web links for GO terms. 
}
  \item{nwStatsControls}{
a single character value specifying the name of the controls to be used
as a control population in the two-sample tests (this HAS to be corresponding
to how these control wells have been annotated in the column "controlStatus"
of the fData(cellHTSobject) dataframe). If nothing is specified, the function 
will look for negative controls labelled "neg".
}
  \item{nwStatsAlternative}{
a single character value specifying the alternative hypothesis: "two.sided", 
"less" or "greater"
}
  \item{nwStatsTests}{
a single character value specifying the tests to be performed: "T-test", 
"MannWhitney" or "RankProduct". If nothing is specified, all three tests will 
be performed. Be aware that the Rank Product test is slower than the other 
two, and returns a percent false discovery (equivalent to a FDR, not a p-value).
}
  \item{nwStatsColumns}{
a character vector of any (relevant, i.e. that is produced in the tests) 
combination of "t.test.pvalues.two.samples", "t.test.pvalues.one.sample", 
"mannW.test.pvalues.one.sample", "mannW.test.pvalues.two.samples", 
"rank.product.pfp.greater", "rank.product.pfp.less"
}
  \item{nwAnalysisFdr}{
a single numeric value specifying the FDR used in the networkAnalysis 
function for the scores calculation
}
  \item{nwAnalysisGenetic}{
a single logical value indicating if the genetic interaction data of
the Biogrid dataset were kept in the network analysis
}
  \item{interactionMatrix}{
an interaction matrix including columns 'InteractionType', 'InteractorA'
and 'InteractorB'. If this matrix is available, the interactome can be
directly built based on it.
}
  \item{nwAnalysisOrder}{
the order used in the networkAnalysis function for the scores calculation
}
  \item{ntop}{
the number of plots to be produced for the GSEA analysis.  For each gene
set collection, plots are produced for the "nplots" most significant
p-values.  
}
  \item{allSig}{
a single logical value indicating whether or not to generate plots for
all significant gene sets. A gene set is significant if its corresponding
adjusted p-value is less than the \code{pValueCutoff} set in function
\code{analyze} (see function \code{analyze} for more details).   
 }
  \item{reportDir}{
a single character value specifying the directory to store reports
}
  \item{verbose}{
a single logical value indicating to display detailed messages (when verbose= 
TRUE) or not (when verbose=FALSE)
}
}
\details{
This pipeline function performs gene set over-representation, enrichment 
analyses and identifies enriched subnetworks, writes a report and saves results 
as 'RData' in the user-specified directory. 

For the KEGG and GO gene sets, this function produces links to the relevant 
database entries. To make the report more readable, the user can append a column 
named 'Gene.Set.Term' to each dataframe in the result slot of the object of 
class \code{\link[HTSanalyzeR:GSCA]{GSCA}} using the method \code{\link[HTSanalyzeR:appendGSTerms]{appendGSTerms}}.
To ensure that the function works properly, the gene set ID should comply with 
the following requirements:   
\itemize{
\item for KEGG 3 letters (species ID) followed by a set of numbers and then 
any number of non-digit characters
\item for GO "GO:" immediately followed by the GO term number (digits only) 
and any number of non-digits characters.
}
}
\value{
Produce a set of html pages and save \code{\link[HTSanalyzeR:GSCA]{GSCA}} and \code{\link[HTSanalyzeR:NWA]{NWA}} objects. The report
starts from the page called "index.html" to navigate those pages.
}
\author{
Camille Terfve, Xin Wang
}
\examples{
\dontrun{
library(org.Dm.eg.db)
library(GO.db)
library(KEGG.db)
library(cellHTS2)
library(BioNet)
library(igraph)
#prepare data from cellHTS2
experimentName <- "KcViab"
dataPath <- system.file(experimentName, package = "cellHTS2")
x <- readPlateList("Platelist.txt", name = experimentName, path=dataPath,
verbose=TRUE)
x <- configure(x, descripFile = "Description.txt", confFile = 
"Plateconf.txt",logFile = "Screenlog.txt", path = dataPath)
xn <- normalizePlates(x, scale = "multiplicative", log = FALSE, method =
"median", varianceAdjust = "none")
xn <- annotate(xn, geneIDFile = "GeneIDs_Dm_HFA_1.1.txt",path = dataPath)
cellHTS2DrosoData<-xn
#prepare a list of gene set collections
GO_MF <- GOGeneSets(species = "Dm",ontologies = c("MF"))
PW_KEGG <- KeggGeneSets(species = "Dm")
gscList <- list(GO_MF = GO_MF, PW_KEGG = PW_KEGG)
#note: set computer cluster here to do parallel computing
#options(cluster=makeCluster(4,"SOCK"))
HTSanalyzeR4cellHTS2(
	normCellHTSobject=cellHTS2DrosoData,
	annotationColumn="GeneID",
	species=c("Dm"),
	initialIDs="FlybaseCG",
	listOfGeneSetCollections=gscList,
	cutoffHitsEnrichment=2,
	minGeneSetSize=200,
	keggGSCs="PW_KEGG",
	goGSCs=c("GO_MF"),
	allSig=TRUE,
	reportDir="HTSanalyzerReport",
	verbose=TRUE
)
#note: stop the cluster if you created
#if(is(getOption("cluster"),"cluster")) {
#	stopCluster(getOption("cluster"))
#	options(cluster=NULL)
#}
##browse the index page
browseURL(file.path(getwd(), "HTSanalyzerReport", "index.html"))
}
}
