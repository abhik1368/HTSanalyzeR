\name{multiHyperGeoTest}
\alias{multiHyperGeoTest}
\title{
 Hypergeometric test on gene set collections.
}
\description{
This function performs hypergeometric test for over-representation of hits,  on gene set collections. This function applies the phyp.geo function to an entire gene set collection and returns a data frame.
}
\usage{
multiHyperGeoTest(GeneSet, GeneList, hits, min.gene.set.size = 15, 
pAdjustMethod = c("holm", "hochberg", "hommel", "bonferroni", 
"BH", "BY", "fdr", "none"), verbose = TRUE)
}
\arguments{
  \item{GeneSet}{
GeneSetCollection object
}
  \item{GeneList}{
Named of a data vector (names=genes identifiers)
}
  \item{hits}{
Names of hits
}
  \item{min.gene.set.size}{
Minimum size for a gene set to be included in the analysis
}
  \item{pAdjustMethod}{
A multiple testing p value adjustment method
}
  \item{verbose}{
A logical value to decide to display detailed messages (when verbose=TRUE) or not (when verbose=FALSE)
}
}
\value{
A data frame containing the results of the hypergeometric test (one row per gene set)
}
\author{
John C. Rose
}
\seealso{
hyperGeoTest
}
\examples{
gl<-runif(100,min=0,max=5)
gl<-order(gl)
names(gl)<-as.character(sample(x=seq(from=1,to=100,by=1),size=100,replace=FALSE))
gs1<-GeneSet(sample(names(gl),size=20,replace=FALSE), setName="subset1")
gs2<-GeneSet(sample(names(gl),size=20,replace=FALSE), setName="subset2")
gsc<-GeneSetCollection(c(gs1,gs2))
hypgeo<-multiHyperGeoTest(GeneSet=gsc, GeneList=names(gl), 
hits=names(gl)[which(abs(gl) > 2)], min.gene.set.size = 2, 
pAdjustMethod ="BH")
##
\dontrun{
library(org.Dm.eg.db)
library(KEGG.db)
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
data4enrich<-data4enrich[which(!is.na(names(data4enrich)))]
data4enrichentrez<-drosoAnnotationConvertor(geneList=data4enrich,
initialIDs="FlybaseCG")
data4enrichentrez<-data4enrichentrez[which(!is.na(data4enrichentrez))]
data4enrichentrez<-duplicateRemover(geneList=data4enrichentrez,method="max")
kegg.droso<-KeggGeneSets(species="Dm")
hypgeo.kegg<-multiHyperGeoTest(GeneSet=kegg.droso, GeneList=names(data4enrichentrez), 
hits=names(data4enrichentrez)[which(abs(data4enrichentrez) > 2)], 
min.gene.set.size = 15, pAdjustMethod ="BH")
}
}
