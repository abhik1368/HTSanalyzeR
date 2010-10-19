\name{hyperGeoTest}
\alias{hyperGeoTest}
\title{
Performs hypergeometric tests for overrepresentation analysis
}
\description{
This function takes in a single gene set (GeneSet), a vector (GeneList) of gene symbols for all tested genes, a vector of "hits" (hits), and a p-value adjustment method. It outputs a vector containing the size of the gene universe, the size of the gene set within this universe (i.e. how many genes from the universe map to this gene set), the total number of hits, the number of hits expected to occur in the gene set, the actual hits observed in the gene set, and the pvalue from a hypergeometric test.
}
\usage{
hyperGeoTest(GeneSet, GeneList, hits)
}
\arguments{
  \item{GeneSet}{
A vector of names of genes in a gene set
}
  \item{GeneList}{
A vector of gene identifiers (usually all genes tested in a screen)
}
  \item{hits}{
A vector of gene identifiers for those genes considered as hits
}
}
\value{
A vector containing the size of the gene universe, the size of the gene set within this universe (i.e. how many genes from the universe map to this gene set), the total number of hits, the number of hits expected to occur in the gene set, the actual hits observed in the gene set, and the pvalue from a hypergeometric test.
}
\author{
John C. Rose
}
\seealso{
multiHyperGeoTest
}
\examples{
gl<-runif(100,min=0,max=5)
gl<-order(gl)
names(gl)<-as.character(sample(x=seq(from=1,to=100,by=1),size=100,replace=FALSE))
gs1<-sample(names(gl),size=20,replace=FALSE)
gs2<-sample(names(gl),size=20,replace=FALSE)
gsc<-list(subset1=gs1,subset2=gs2)
hypgeo<-hyperGeoTest(GeneSet=gsc[["subset1"]], GeneList=names(gl), 
hits=names(gl)[which(abs(gl) > 2)])
##
#Prepare the data
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
#Remove NA's and duplicates
data4enrichentrez<-data4enrichentrez[which(!is.na(data4enrichentrez))]
data4enrichentrez<-duplicateRemover(geneList=data4enrichentrez,method="max")
#Build a gene set collection
kegg.droso<-KeggGeneSets(species="Dm")
#Do the tests
hypgeo.kegg<-hyperGeoTest(GeneSet=kegg.droso[[1]], 
GeneList=names(data4enrichentrez), 
hits=names(data4enrichentrez)[which(abs(data4enrichentrez) > 2)])
}
