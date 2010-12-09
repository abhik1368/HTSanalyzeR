\name{GSCA-class}
\alias{GSCA}
\alias{GSCA-class}
\alias{writeHits}
\alias{writeHits,GSCA-method}
\docType{class}
\title{
An S4 class for Gene Set Collection Analyses on high-throughput screens
}
\description{
This S4 class includes a series of methods to do gene set enrichment 
analysis and hypergeometric tests for high-throughput screens. 
}
\section{Objects from the Class}{ 
Objects of class \code{GSCA} can be created from \code{
new("GSCA", listOfGeneSetCollections, geneList, hits)} (see
the examples below)
}
\section{Slots}{
	\describe{	
		\item{\code{listOfGeneSetCollections}:}{
		a list of gene set collections (a 'gene set collection' is a list 
		of gene sets). 
		}
		
		\item{\code{geneList}:}{
		a numeric or integer vector of phenotypes named by gene identifiers.
		}
		
		\item{\code{hits}:}{
		a character vector of the gene identifiers (used as hits in the hypergeometric tests).
		}
		
		\item{\code{para}:}{
		a list of parameters for hypergeometric tests and GSEA. These 
		parameters are \code{pValueCutoff}, \code{pAdjustMethod}, \code{nPermutations}, 
		\code{minGeneSetSize} and \code{exponent} (see function
		\code{\link[HTSanalyzeR:analyzeGeneSetCollections]{analyzeGeneSetCollections}} for detailed descriptions about
		these parameters). 
		}
		
		\item{\code{result}:}{
		a list of results (see the returned values in the function 
		\code{\link[HTSanalyzeR:analyzeGeneSetCollections]{analyzeGeneSetCollections}}).
		}
		
		\item{\code{summary}:}{
		a list of summary information for \code{listOfGeneSetCollections},
		\code{geneList}, \code{hits}, \code{para}, and \code{result}.
		}
		
		\item{\code{preprocessed}:}{
		a single logical value specifying whether or not the input data has been 
		preprocessed.
		}
	}
}
\section{Methods}{
	An overview of methods with class-specific functionality:
	More detailed introduction can be found in help for each specific function.
	\describe{	
		\item{\code{preprocess}}{ 
		do preprocessing on input vectors of phenotypes and hits including: 
		a) removing NAs in the \code{geneList} and \code{hits}; 
		b) invoking function \code{\link[HTSanalyzeR:duplicateRemover]{duplicateRemover}} to process
		duplicated phenotypes (see \code{\link[HTSanalyzeR:duplicateRemover]{duplicateRemover}} for more details); 
		c) invoking function \code{\link[HTSanalyzeR:annotationConvertor]{annotationConvertor}} to convert annotations; 
		d) ranking phenotypes in a decreasing order.
		}
		\item{\code{analyze}}{
		perform hypergeometric tests and Gene Set Enrichment Analysis based on
		input parameter list \code{para}.
		}
		\item{\code{appendGSTerms}}{
		append gene set terms to GSCA results
		}
		\item{\code{summarize}}{
		print summary information about \code{listOfGeneSetCollections},
		\code{geneList}, \code{hits}, \code{para}, and \code{result}.
		}
		\item{\code{getTopGeneSets}}{
		select top significant gene sets from \code{object@results$`resultName`} 
		by setting \code{ntop} or \code{allSig}.
		}
		\item{\code{writeHits}}{
		write observed hits in gene sets for hypergeometric tests.
		}
		\item{\code{viewGSEA}}{
		view a figure of GSEA results for a gene set in a gene set collection.
		}
		\item{\code{plotGSEA}}{
		plot and save figures of GSEA results for top significant gene sets in a gene set collection.
		}
		\item{\code{viewEnrichMap}}{
		plot an enrichment map for GSEA or hypergeometric test results
		}
		\item{\code{plotEnrichMap}}{
		plot and save an enrichment map for GSEA or hypergeometric test results
		}
		\item{\code{report}}{
		generate html reports.
		}
	}
}
\author{Xin Wang \email{xw264@cam.ac.uk}}
\seealso{
  \code{\link[HTSanalyzeR:preprocess]{preprocess}}
  \code{\link[HTSanalyzeR:analyze]{analyze}}
  \code{\link[HTSanalyzeR:appendGSTerms]{appendGSTerms}}
  \code{\link[HTSanalyzeR:summarize]{summarize}}
  \code{\link[HTSanalyzeR:getTopGeneSets]{getTopGeneSets}}
  \code{\link[HTSanalyzeR:writeHits]{writeHits}}
  \code{\link[HTSanalyzeR:viewGSEA]{viewGSEA}}
  \code{\link[HTSanalyzeR:plotGSEA]{plotGSEA}}
  \code{\link[HTSanalyzeR:viewEnrichMap]{viewEnrichMap}}
  \code{\link[HTSanalyzeR:plotEnrichMap]{plotEnrichMap}}
  \code{\link[HTSanalyzeR:report]{report}}
}
\examples{
\dontrun{
library(org.Dm.eg.db)
library(KEGG.db)
library(AnnotationDbi)
library(igraph)
##load data for enrichment analyses
data("KcViab_Data4Enrich")
##select hits
hits <- names(KcViab_Data4Enrich)[which(abs(KcViab_Data4Enrich) > 2)]
##set up a list of gene set collections
PW_KEGG <- KeggGeneSets(species = "Dm")
gscs <- list(PW_KEGG = PW_KEGG)
##create an object of class 'GSCA'
gsca <- new("GSCA", listOfGeneSetCollections=gscs, geneList = 
KcViab_Data4Enrich, hits = hits)
##do preprocessing (KcViab_Data4Enrich has already been preprocessed)
gsca <- preprocess(gsca, species="Dm", initialIDs = "Entrez.gene", 
keepMultipleMappings = TRUE, duplicateRemoverMethod = "max", 
orderAbsValue = FALSE)
##do hypergeometric tests and GSEA
gsca <- analyze(gsca, para = list(pValueCutoff = 0.05, pAdjustMethod 
= "BH", nPermutations = 1000, minGeneSetSize = 60, exponent = 1))
##print summary information
summarize(gsca)
##get all significant gene sets in "PW_KEGG"
sigGSs<-getTopGeneSets(gsca, "GSEA.results", "PW_KEGG", allSig=TRUE)
##view a GSEA figure
viewGSEA(gsca, gscName="PW_KEGG", gsName=sigGSs[["PW_KEGG"]][1])
dev.off()
##append gene set terms to results
gsca<-appendGSTerms(gsca, keggGSCs="PW_KEGG")
##view an enrichment map for GSEA results
eb<-viewEnrichMap(gsca, gscs="PW_KEGG", allSig=TRUE, gsNameType="term", 
displayEdgeLabel=FALSE, layout="layout.fruchterman.reingold")
##write html reports
report(object = gsca, experimentName = "GSCATest", species = "Dm", 
allSig = TRUE, keggGSCs = "PW_KEGG", reportDir="GSCATestReport")
##browse the index page
browseURL(file.path(getwd(), "GSCATestReport", "index.html"))
}
}
\keyword{classes}
