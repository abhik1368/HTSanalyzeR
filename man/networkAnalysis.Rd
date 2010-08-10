\name{networkAnalysis}
\alias{networkAnalysis}
\title{
Wrapper function around the network analysis process in this package.
}
\description{
This function finds subnetworks enriched for genes with significant phenotypes.  It is a wrapper around the whole network analysis process in this package.
}
\usage{
networkAnalysis(cellHTSobject, annotationColumn = "GeneID", controls = "neg", 
alternative = c("two.sided", "less", "greater"), logged = FALSE, 
tests = c("T-test"), 
columns = c("t.test.pvalues.two.samples", "t.test.pvalues.one.sample"), 
species = c("Drosophila_melanogaster"), initialIDs = "FlybaseCG", 
fdr = 0.001, genetic = FALSE, phenotypeVector = NA, biogridObject = NA, 
order = 2, link)
}
\arguments{
  \item{cellHTSobject}{
An object of class cellHTS.
}
  \item{annotationColumn}{
The name of the column in the fData(cellHTSobject) data frame from which the feature identifiers will be extracted
}
  \item{controls}{
The name of the controls to be used as a control population in the 2 samples tests (this HAS to be corresponding to how these control wells have been annotated in the column "controlStatus" of the fData(cellHTSobject) dataframe) if nothing is specified, the function will look for negative controls labelled "neg"
}
  \item{alternative}{
two.sided/less/greater (only one of those 3): which is the alternative hypothesis to be tested
}
  \item{logged}{
TRUE or FALSE: has the data been logged during the normalization process
}
  \item{tests}{
Any combination of "T-test","MannWhitney","RankProduct"
}
  \item{columns}{
Any (relevant, i.e. that is produced in the tests) combination of "t.test.pvalues.two.samples", "t.test.pvalues.one.sample", "mannW.test.pvalues.one.sample", "mannW.test.pvalues.two.samples", "rank.product.pfp.greater", "rank.product.pfp.less"
}
  \item{species}{
One of the following: "Drosophila_melanogaster", "Homo_sapiens", "Rattus_norvegicus", "Mus_musculus", "Caenorhabditis_elegans"
}
  \item{initialIDs}{
The type of identifiers in the column specified by annotationColumn these have to be one of: "Ensembl.transcript", "Ensembl.prot", "Ensembl.gene", "Entrez.gene", "RefSeq", "Symbol", "GenBank", "Flybase", "FlybaseCG", "FlybaseProt"

}
  \item{fdr}{
The fdr for the scoring of nodes (see Dittrich et al., 2008)
}
  \item{genetic}{
FALSE or TRUE: should genetic interactions be kept? if nothing is specified, they will be removed from the data
}
  \item{phenotypeVector}{
A named (entrez ID) vector of phenotypes, if available if a phenotype vector is given, the these will be used to colour the nodes of the network the nodes that are not in the screen (no phenotype value) will be square nodes, the other ones will be circles
}
  \item{biogridObject}{
Do you already have a biogrid cellHTSobject (as created by the function "biogrid") if nothing is specified this is set to NA and teh biogrid cellHTSobject will be build. If you do give a biogrid cellHTSobject there, trhe parameters "link" and "species" are not mandatory
}
  \item{order}{
The order of the statistic for the aggregation of pvalues (see Dittrich et al., 2008)
}
  \item{link}{
The link to the biogrid data in tab2 format (on the biogrid website) if nothing is provided the link is the one valid on June 7th 2010.	
}
}
\details{
This function is made to work on the tab2 format from the Biogrid (i.e. the first line is a header containing the columns names, and the identifiers used are Entrez identifiers).
}
\value{
A matrix with one row for each interaction, and three columns: InteractorA, InteractorB (both given by their Entrez Identifiers) and InteractionType (physical or genetic).
}
\references{
Beisser D, Klau GW, Dandekar T, Muller T, Dittrich MT. BioNet: an R-Package for the functional analysis of biological networks. Bioinformatics. 2010 Apr 15;26(8):1129-30.

Dittrich MT, Klau GW, Rosenwald A., Dandekar T and Muller T. \emph{Identifying functional modules in protein-protein interaction networks: an integrated exact approach.} Bioinformatics 2008 24(13):i223-i231.
}
\author{
Camille Terfve
}
\seealso{
enricheSubNw
}
\examples{
\dontrun{
library(org.Dm.eg.db)
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
#
xsc<-scoreReplicates(xn,sign="-",method="zscore")
xsc<-summarizeReplicates(xsc,summary="mean")
data4enrich<-as.vector(Data(xsc));
names(data4enrich)<-fData(xsc)[,"GeneID"] 
data4enrich<-data4enrich[-which(is.na(names(data4enrich))==TRUE)]
data4enrichentrez<-drosoAnnotationConvertor(geneList=data4enrich,
initialIDs="FlybaseCG")
data4enrichentrez<-data4enrichentrez[-which(is.na(data4enrichentrez))]
data4enrichentrez<-duplicateRemover(geneList=data4enrichentrez,method="max")
#
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
}
}
