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
Version: \tab 1.2.4\cr
Date: \tab 2010-10-19\cr
License: \tab Artistic-2.0\cr
LazyLoad: \tab yes\cr
	}
}

\author{
Camille D.A. Terfve , Xin Wang, John C. Rose, and Florian Markowetz
Maintainer: Xin Wang <Xin.Wang@cancer.org.uk>
}

\examples{
\dontrun{
library(org.Dm.eg.db)
library(GO.db)
library(KEGG.db)
library(cellHTS2)

verbose=TRUE
#Prepare the data
cat("-Performing pre-processing ... \n")
experimentName <- "KcViab"
dataPath <- system.file(experimentName, package = "cellHTS2")
x <- readPlateList("Platelist.txt", name = experimentName,path = dataPath,verbose=verbose)
x <- configure(x, descripFile = "Description.txt", confFile = "Plateconf.txt",logFile = "Screenlog.txt", path = dataPath)
xn <- normalizePlates(x, scale = "multiplicative", log = FALSE,method = "median", varianceAdjust = "none")
xn <- annotate(xn, geneIDFile = "GeneIDs_Dm_HFA_1.1.txt",path = dataPath)
xsc <- scoreReplicates(xn, sign = "-", method = "zscore")
xsc <- summarizeReplicates(xsc, summary = "mean")	
	
data4enrich <- as.vector(Data(xsc))
names(data4enrich) <- fData(xsc)[, "GeneID"]
data4enrich <- data4enrich[which(!is.na(names(data4enrich)))]
data4enrichentrez <- drosoAnnotationConvertor(geneList = data4enrich,initialIDs = "FlybaseCG",verbose=verbose)
data4enrichentrez <- data4enrichentrez[which(!is.na(data4enrichentrez))]
data4enrichentrez <- duplicateRemover(geneList = data4enrichentrez,method = "max")
data.hits <- names(data4enrichentrez)[which(abs(data4enrichentrez) > 2)]
Dm.GO.MF <- GOGeneSets(species = "Dm",ontologies = c("MF"))
Dm.GO.BP <- GOGeneSets(species = "Dm",ontologies = c("BP"))
Dm.GO.CC <- GOGeneSets(species = "Dm",ontologies = c("CC"))
Dm.PW.KEGG <- KeggGeneSets(species = "Dm")
	
gsc.list <- list(Dm.GO.MF = Dm.GO.MF, Dm.GO.BP = Dm.GO.BP, Dm.GO.CC = Dm.GO.CC, Dm.PW.KEGG = Dm.PW.KEGG)
	
gsc.solo.list <- list(Dm.GO.MF = Dm.GO.MF)
	
cat("-Pre-processing complete \n\n")
#if R pakage 'snow' is available, users can choose to create a cluster object called "htsa.cluster" for much faster speed 
#options(htsa.cluster=makeCluster(4,"SOCK"))
#gene set enrichment analyses
enrichment.analysis <- analyzeGeneSetCollections(ListOfGeneSetCollections = gsc.list,
GeneList = data4enrichentrez, hits = data.hits, p.value.cutoff = 0.05,
pAdjustMethod = "BH", npermutations = 1000, min.gene.set.size = 100,
exponent = 1, whichSetIsKEGGIds = 4, whichSetIsGOIds = c(1, 2, 3),verbose=verbose)

if(is(getOption("htsa.cluster"),"cluster")) 
	stopCluster(getoption("htsa.cluster"))
if(!("BioNet" %in% search())) library(BioNet) 
#network analyses		
test.module <- networkAnalysis(cellHTSobject = xn, annotationColumn = "GeneID",
controls = "neg", alternative = "two.sided", logged = FALSE,
tests = c("T-test"), columns = c("t.test.pvalues.two.samples",
"t.test.pvalues.one.sample"), species = c("Dm"),
initialIDs = "FlybaseCG", fdr = 0.001, genetic = FALSE,
biogridObject = NA, order = 2, reportdir="HTSanalyzerReport",verbose=verbose)
#Report
writeReportHTSA(experimentName = "cellHTS2DrosoData",
enrichmentAnalysis = enrichment.analysis, cutoffHits = 2,
hits = names(data4enrichentrez)[which(abs(data4enrichentrez) >2)], listOfGeneSetCollections = gsc.list, geneList = data4enrichentrez,
p.adj.method = "BH", nPermutations = 1000, min.gene.set.size = 100,
exponent = 1, nwAnalysisOutput = test.module, 
controls = "neg", alternative = "two.sided", tests = c("T-test"),
columns = c("t.test.pvalues.two.samples", "t.test.pvalues.one.sample"),
species = c("Dm"), fdr = 0.001,
genetic = FALSE, networkObject = NA, nGseaPlots = 2,
geneListName = "data4enrichentrez", whichSetIsKEGGIds = 4,
whichSetIsGOIds = c(1, 2, 3),reportdir="HTSanalyzerReport")
	}
}
