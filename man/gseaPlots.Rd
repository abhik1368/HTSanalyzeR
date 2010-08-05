\name{gseaPlots}
\alias{gseaPlots}
\title{
Produce a collection of plots for GSEA.
}
\description{
Function that takes the output of gseaScores and the gene list, and does the plots for GSEA.
}
\usage{
gseaPlots(Running.Score, Enrichment.Score, Positions, geneList, 
output = "png", name)
}
\arguments{
  \item{Running.Score}{
Output of the function "gseascores"
}
  \item{Enrichment.Score}{
Output of the function "gseascores"
}
  \item{Positions}{
Output of the function "gseascores"
}
  \item{geneList}{
A gene list as in "gseascores"
}
  \item{output}{
"pdf" or "png"
}
  \item{name}{
Name of the gene set for which these plots are produced
}
}
\details{
The plots are going to be produced in the current working directory.
}
\value{
This function does not return anything, it just produces an image file containing the plots.
}
\references{
Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S. & Mesirov, J. P. (2005)  \emph{Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles.} Proc. Natl. Acad. Sci. USA 102, 15545-15550.
}
\author{
Camille Terfve
}
\examples{
\dontrun{
library(KEGG.db)
library(org.Dm.eg.db)
library(cellHTS2)
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
data4enrich<-data4enrich[-which(is.na(names(data4enrich))==TRUE)]
data4enrichentrez<-drosoAnnotationConvertor(geneList=data4enrich,
initialIDs="FlybaseCG")
data4enrichentrez<-data4enrichentrez[-which(is.na(data4enrichentrez))]
data4enrichentrez<-duplicateRemover(geneList=data4enrichentrez,method="max")
kegg.droso<-KeggGeneSets(species="Drosophila_melanogaster")
test<-gseaScores(geneList=data4enrichentrez,geneSet=geneIds(kegg.droso[[1]]),
exponent=1,mode="graph")
gseaPlots(Running.Score=test$Running.Score,Enrichment.Score=test$Enrichment.Score,
Positions=test$Positions,geneList=data4enrichentrez,name=names(kegg.droso)[1])
}
}
