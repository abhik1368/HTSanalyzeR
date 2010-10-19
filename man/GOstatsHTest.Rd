\name{GOstatsHTest}
\alias{GOstatsHTest}
\title{
Implements the hyperGTest function from GOStats.
}
\description{
This function is a slightly simplified version of the hyperGTest function from GOStats.  This function produces a GOHypergParams object and then implements the hyperGTest function from GOStats.
}
\usage{
GOstatsHTest(hits, GeneList, annotation, ontology = c("BP", "CC", "MF"), 
pvalueCutoff, conditional = FALSE, testDirection = "over")
}
\arguments{
  \item{hits}{
Vector of gene symbols -- hits genes
}
  \item{GeneList}{
Vector of gene symbols -- all tested genes (GeneList)
}
  \item{annotation}{
The annotation data package to be used
}
  \item{ontology}{
A string indicating the desired ontology to  test (Cellular Component, Biological Process, or Molecular Function)
}
  \item{pvalueCutoff}{
A cutoff for p values called significant
}
  \item{conditional}{
If conditional is TRUE, then the test uses GOstats' conditional algorithm which accounts for the topology of the GOterms being tested.
}
  \item{testDirection}{
Test direction ("over","under") determines whether the test identifies over or underrepresented GO terms (see GOStats vignette, pg 6).
}
}
\value{
A GOHyperGResult instance.
}
\references{
S. Falcon  and R. Gentleman. \emph{Using GOstats to test gene lists for GO term association.} Bioinformatics 2007 23(2):257-258.
}
\author{
John C. Rose
}
\examples{
\dontrun{
library(GOstats)
library(org.Dm.eg.db)
library(cellHTS2)
experimentName<-"KcViab"
dataPath<-system.file(experimentName, package="cellHTS2")
x<-readPlateList("Platelist.txt",name=experimentName,path=dataPath)
x<-configure(x,descripFile="Description.txt",confFile="Plateconf.txt",logFile="Screenlog.txt",path=dataPath)
xn<-normalizePlates(x,scale="multiplicative",log=FALSE,method="median",varianceAdjust="none")
xn<-annotate(xn,geneIDFile="GeneIDs_Dm_HFA_1.1.txt",path=dataPath)
xsc<-scoreReplicates(xn,sign="-",method="zscore")
xsc<-summarizeReplicates(xsc,summary="mean")
data4enrich<-as.vector(Data(xsc));
names(data4enrich)<-fData(xsc)[,"GeneID"] 
data4enrich<-data4enrich[which(!is.na(names(data4enrich)))]
data4enrichentrez<-drosoAnnotationConvertor(geneList=data4enrich,initialIDs="FlybaseCG")
data4enrichentrez<-data4enrichentrez[which(!is.na(data4enrichentrez))]
data4enrichentrez<-duplicateRemover(geneList=data4enrichentrez,method="max")
gost<-GOstatsHTest(hits=names(data4enrichentrez)[which(abs(data4enrichentrez)>2)],
GeneList=names(data4enrichentrez),pvalueCutoff=0.05,annotation="org.Dm.eg.db",
ontology="CC")
}
}
