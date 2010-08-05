\name{HTSanalyzeR-package}
\alias{HTSanalyzeR-package}
\alias{HTSanalyzeR_package}
\docType{package}

\title{
This package does enrichment and network analysis on high throughput screens such as RNAi screens.
}

\description{
This package does enrichment (hypergeometric tests and GSEA) and network analysis (using the package BioNet) on high throughput screens  such as RNAi screens.
}

\details{
	\tabular{ll}{
Package: \tab HTSAnalyzeR\cr
Type: \tab Package\cr
Version: \tab 1.1.6\cr
Date: \tab 2010-07-13\cr
License: \tab Artistic-2.0\cr
LazyLoad: \tab yes\cr
	}
}

\author{
John C. Rose, Camille D.A. Terfve , Xin Wang and Florian Markowetz
Maintainer: Xin Wang <Xin.Wang@cancer.org.uk>
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
xsc<-scoreReplicates(xn,sign="-",method="zscore")
xsc<-summarizeReplicates(xsc,summary="mean")
data4enrich<-as.vector(Data(xsc));
names(data4enrich)<-fData(xsc)[,"GeneID"] 
data4enrich<-data4enrich[-which(is.na(names(data4enrich))==TRUE)]
data4enrichentrez<-drosoAnnotationConvertor(geneList=data4enrich,
initialIDs="FlybaseCG")
data4enrichentrez<-data4enrichentrez[-which(is.na(data4enrichentrez))]
data4enrichentrez<-duplicateRemover(geneList=data4enrichentrez,method="max")
#Make the list of gene set collections
Dm.GO.MF<-GOGeneSets(species="Drosophila_melanogaster",ontologies=c("MF"))
Dm.GO.BP<-GOGeneSets(species="Drosophila_melanogaster",ontologies=c("BP"))
Dm.GO.CC<-GOGeneSets(species="Drosophila_melanogaster",ontologies=c("CC"))
kegg.droso<-KeggGeneSets(species="Drosophila_melanogaster");
gsc.list<-list(Dm.GO.MF=Dm.GO.MF,Dm.GO.BP=Dm.GO.BP,Dm.GO.CC=Dm.GO.CC,
kegg.droso=kegg.droso)
#Run the enrichment analysis
enrichment.analysis<-analyzeGeneSetCollections(
	ListOfGeneSetCollections=gsc.list,
	GeneList=data4enrichentrez,
	hits=names(data4enrichentrez)[which(abs(data4enrichentrez)>2)],
	pAdjustMethod="BH",
	npermutations=1000,
	min.gene.set.size=15,
	exponent=1,
	whichSetIsKEGGIds=4,
	whichSetIsGOIds=c(1,2,3)
	)
#Run the networks analysis
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
#Write the report
writeReportHTSA(
	experimentName="cellHTS2DrosoData",
	enrichmentAnalysis=enrichment.analysis,
	cutoffHits=2,
	hits=names(data4enrichentrez)[which(abs(data4enrichentrez)>2)],
	listOfGeneSetCollections=gsc.list,
	geneList=data4enrichentrez,
	p.adj.method="BH",
	nPermutations=1000,
	min.gene.set.size=15,
	exponent=1,
	nwAnalysisOutput=test.module,
	nwAnalysisGraphFile="EnrichedSubNw.png",
	controls="neg",
	alternative="two.sided",
	tests=c("T-test"),
	columns=c("t.test.pvalues.two.samples","t.test.pvalues.one.sample"),
	species=c("Drosophila_melanogaster"),
	fdr=0.001,
	genetic=FALSE,
	biogridObject=NA,
	nGseaPlots=10,
	geneListName="data4enrichentrez",
	whichSetIsKEGGIds=4,
	whichSetIsGOIds=c(1,2,3)
	)
	}
}
