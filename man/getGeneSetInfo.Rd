\name{getGeneSetInfo}
\alias{getGeneSetInfo}
\title{
GSEA associated information about a gene set
}
\description{
This function returns a data frame of the genes in a selected gene set, their gene symbols, numeric phenotypes, phenotype rank, and whether they are a member of the leading edge subset (indicated by a 1 in the pertinent column).
}
\usage{
getGeneSetInfo(GeneList, GeneSet, exponent = 1,species=c("Dm","Hs","Rn","Mm","Ce"))
}
\arguments{
  \item{GeneList}{
A named vector of phenotypes, ordered, with no duplicates nor NA values
}
  \item{GeneSet}{
A gene set in GeneSet format
}
  \item{exponent}{
An exponent to weight the phenotypes in a GSEA
}
  \item{species}{
One of "Dm" ("Drosophila_melanogaster"), "Hs" ("Homo_sapiens"), "Rn" ("Rattus_norvegicus"), "Mm" ("Mus_musculus"), "Ce" ("Caenorhabditis_elegans")
}
}
\details{
This function assumes that the initial identifiers of both the gene list and the gene set are Entrez gene identifiers.
}
\value{
A data frame of the genes in a selected gene set, their gene symbols, numeric phenotypes, phenotype rank, and whether they are a member of the leading edge subset (indicated by a 1 in the pertinent column). 
}

\author{
John C. Rose and Camille Terfve
}
\examples{
library(org.Dm.eg.db)
gl<-runif(100,min=0,max=5)
gl<-gl[order(gl,decreasing=TRUE)]
names(gl)<-names(as.list(org.Dm.egSYMBOL))[1:100]
gs1<-sample(names(gl),size=20,replace=FALSE)
gs2<-sample(names(gl),size=20,replace=FALSE)
gsc<-list(subset1=gs1,subset2=gs2)
info<-getGeneSetInfo(GeneList=gl,GeneSet=gsc[["subset1"]],exponent=1,species="Dm")
####
\dontrun{
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
data4enrich<-data4enrich[which(!is.na(names(data4enrich)))]
data4enrichentrez<-drosoAnnotationConvertor(geneList=data4enrich,
initialIDs="FlybaseCG")
data4enrichentrez<-data4enrichentrez[which(!is.na(data4enrichentrez))]
data4enrichentrez<-duplicateRemover(geneList=data4enrichentrez,method="max")
kegg.droso<-KeggGeneSets(species="Dm")
info<-getGeneSetInfo(GeneList=data4enrichentrez,GeneSet=kegg.droso[[1]],exponent=1,species="Dm")
}
}
