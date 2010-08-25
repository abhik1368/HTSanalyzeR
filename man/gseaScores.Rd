\name{gseaScores}
\alias{gseaScores}
\title{
Compute enrichment scores for GSEA
}
\description{
Function to compute enrichment scores for GSEA, running score and position of hits for a gene set.
}
\usage{
gseaScores(geneList, geneSet, exponent, mode = c("graph", "score"))
}
\arguments{
  \item{geneList}{
An ordered list of genes as vector with names = gene IDs (no NAs, no duplicates)
}
  \item{geneSet}{
A gene set (no names, just a vector of characters corresponding to the IDs)
}
  \item{exponent}{
Exponent to weight the phenotypes in the GSEA score calculation
}
  \item{mode}{
"score" or "graph": should only a score be returned, or all the necessary elements to make a plot (see gseaplots)
}
}
\details{
The type of identifiers used in the gene sets and gene list must obviously match.
}
\value{
\item{Enrichment.Score}{the enrichment score}
\item{Running.Score}{running enrichment score (only in mode "graph")}
\item{Positions}{positions in the gene list of the genes in the gene set (only in mode "graph")}
}
\references{
Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S. & Mesirov, J. P. (2005)  \emph{Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles.} Proc. Natl. Acad. Sci. USA 102, 15545-15550.
}
\author{
Camille Terfve
}
\examples{
gl<-runif(100,min=0,max=5)
gl<-gl[order(gl,decreasing=TRUE)]
names(gl)<-as.character(sample(x=seq(from=1,to=100,by=1),size=100,replace=FALSE))
gs<-sample(names(gl),size=20,replace=FALSE)
gsea<-gseaScores(geneList=gl,geneSet=gs,mode="score",exponent=1)
#
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
}
}
