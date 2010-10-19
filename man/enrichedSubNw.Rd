\name{enrichedSubNw}
\alias{enrichedSubNw}
\title{
Searching for significant subnetworks.
}
\description{
This function takes a matrix of p-values generated for a set of genes and searches for enriched sub-networks in an interaction data set from The Biogrid.
}
\usage{
enrichedSubNw(pvaluesMatrix, columns = c("t.test.pvalues.two.samples", 
"t.test.pvalues.one.sample"), species = c("Dm"), 
fdr = 0.001, genetic = FALSE, biogridObject = NA, link , order = 2,
biogridDataDir="Data", verbose=TRUE)
}
\arguments{
  \item{pvaluesMatrix}{
A matrix of pvalues with rows corresponding to genes and columns to p.values (both named). See details.
}
  \item{columns}{
The columns of the pvaluesMatrix to be used in the computation
}
  \item{species}{
One of the following: "Dm" ("Drosophila_melanogaster"), "Hs" ("Homo_sapiens"), "Rn" ("Rattus_norvegicus"), "Mm" ("Mus_musculus"), 
"Ce" ("Caenorhabditis_elegans") (see the "biogrid" function).
}
  \item{fdr}{
The fdr for the scoring of nodes (see the BioNet package and documentation)
}
  \item{genetic}{
FALSE-TRUE: should genetic interactions be kept? If nothing is specified, they will be removed from the data
}
  \item{biogridObject}{
Do you already have a biogrid object (as created by the function "biogrid"). If nothing is specified this is set to NA and the biogrid object will be build. If you do give a biogrid object there, the parameters "link" and "species" are not mandatory
}
  \item{link}{
The link to the biogrid data (on the biogrid website) (see the "biogrid" function)
}
  \item{order}{
The order of the statistic for aggregation of pvalues (see Dittrich et al., 2008).
}
  \item{dataDirectory}{
The directory to store downloaded file.
}
  \item{verbose}{
A logical value to decide to display detailed messages (when verbose=TRUE) or not (when verbose=FALSE)
}
}
\details{
pvaluesMatrix: If you have a vector of p-values just do: 

\code{pvaluesMatrix<-as.matrix(vectorPvalues,ncol=1)}

\code{rownames(pvaluesMatrix)<-names(vectorPvalues)} and input \code{columns=1}
Be careful, the matrix MUST have row names, and these MUST correspond to Entrez Identifiers. Use the annotationconvertor functions if this is not the case.
}
\value{
A graphNEL object consisting of the enriched subnetwork.
}
\references{
Beisser D, Klau GW, Dandekar T, Muller T, Dittrich MT. \emph{BioNet: an R-Package for the functional analysis of biological networks.} Bioinformatics. 2010 Apr 15;26(8):1129-30.

Dittrich MT, Klau GW, Rosenwald A., Dandekar T and Muller T. \emph{Identifying functional modules in protein-protein interaction networks: an integrated exact approach.} Bioinformatics 2008 24(13):i223-i231.
}
\author{
Camille Terfve
}
\seealso{
networkAnalysis, biogridDataDownload
}
\examples{
\dontrun{
library(cellHTS2)
library(Org.Dm.eg.db)
library(BioNet)
experimentName<-"KcViab"
dataPath<-system.file(experimentName, package="cellHTS2")
x<-readPlateList("Platelist.txt",name=experimentName,path=dataPath)
x<-configure(x,descripFile="Description.txt",confFile="Plateconf.txt",
logFile="Screenlog.txt",path=dataPath)
xn<-normalizePlates(x,scale="multiplicative",log=FALSE,method="median",
varianceAdjust="none")
xn<-annotate(xn,geneIDFile="GeneIDs_Dm_HFA_1.1.txt",path=dataPath)
test.stats<-cellHTS2OutputStatTests(cellHTSobject=xn,alternative="two.sided",
tests="T-test")
test.stats.entrez<-drosoAnnotationConvertor(geneList=test.stats,
initialIDs="FlybaseCG",finalIDs="Entrez.gene")
subnw<-enrichedSubNw(pvaluesMatrix=test.stats.entrez,
columns=c("t.test.pvalues.two.samples","t.test.pvalues.one.sample"),
species="Dm",fdr=0.001,genetic=FALSE,
biogridObject=NA,order=2, biogridDataDir="Data",verbose=TRUE)
}
}
