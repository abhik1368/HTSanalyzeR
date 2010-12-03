\name{GSCA-class}
\alias{GSCA}
\alias{GSCA-class}
\alias{writeHits}
\alias{writeHits,GSCA-method}
\docType{class}
\title{
A class for Gene Set Collection Analyses
}
\description{
This class includes methods to do gene set enrichment analysis 
and hypergeometric tests on a vector of phenotypes from experiments
such as RNA interference screening. 
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
		a numeric or integer vector of phenotypes.
		}
		
		\item{\code{hits}:}{
		a character vector of the gene identifiers.
		}
		
		\item{\code{para}:}{
		a list of parameters for hypergeometric tests and GSEA. These 
		parameters are \code{pValueCutoff}, \code{pAdjustMethod}, \code{nPermutations}, 
		\code{minGeneSetSize} and \code{exponent} (see function
		\code{analyzeGeneSetCollections} for detailed descriptions about
		these parameters). 
		}
		
		\item{\code{result}:}{
		a list of results (see the returned values in function 
		\code{analyzeGeneSetCollections}).
		}
		
		\item{\code{summary}:}{
		a list of summary information for \code{listOfGeneSetCollections},
		\code{geneList}, \code{hits}, \code{para}, and \code{result}.
		}
		
		\item{\code{preprocessed}:}{
		a logical value specifying whether or not input data has been 
		preprocessed.
		}
	}
}
\section{Methods}{
	An overview of methods with class-specific functionality:
	More detailed introduction can be found in help for each specific function.
	\describe{	
		\item{\code{preprocess}}{ 
		do preprocessing on input gene list and hit list including: 
		a) removing NAs in \code{geneList} and \code{hits}; 
		b) invoking function \code{duplicateRemover} to process
		duplicated phenotypes (see \code{duplicateRemover} for more details); 
		c) invoking function \code{annotationConvertor} to convert annotations; 
		d) ranking phenotypes in a decreasing order.
		}
		\item{\code{analyze}}{
		do hypergeometric tests and Gene Set Enrichment Analysis based on
		input parameter list \code{para}.
		}
		\item{\code{summarize}}{
		print summary information of \code{listOfGeneSetCollections},
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
		plot a figure of GSEA results of a gene set in a gene set collection.
		}
		\item{\code{plotGSEA}}{
		plot and save figures for Gene Set Enrichment Analysis.
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
  \code{\link[HTSanalyzeR:summarize]{summarize}}
  \code{\link[HTSanalyzeR:getTopGeneSets]{getTopGeneSets}}
  \code{\link[HTSanalyzeR:writeHits]{writeHits}}
  \code{\link[HTSanalyzeR:plotGSEA]{plotGSEA}}
  \code{\link[HTSanalyzeR:report]{report}}
}
\examples{
\dontrun{
library(org.Dm.eg.db)
library(KEGG.db)
##load data for enrichment analyses
data("data4enrich")
##select hits
hits <- names(data4enrich)[which(abs(data4enrich) > 2)]
##set up a list of gene set collections
PW_KEGG <- KeggGeneSets(species = "Dm")
gscs <- list(PW_KEGG = PW_KEGG)
##create an object of class 'GSCA'
gsca <- new("GSCA", listOfGeneSetCollections=gscs, geneList = data4enrich,
hits = hits)
##do preprocessing (data4enrich has already been preprocessed)
gsca <- preprocess(gsca, species="Dm", initialIDs = "Entrez.gene", 
keepMultipleMappings = TRUE, duplicateRemoverMethod = "max", 
orderAbsValue = FALSE)
##do hypergeometric tests and GSEA
gsca <- analyze(gsca, para = list(pValueCutoff = 0.05, pAdjustMethod 
= "BH", nPermutations = 1000, minGeneSetSize = 100, exponent = 1))
##print summary information
summarize(gsca)
##get all significant gene sets in "PW_KEGG"
sigGSs<-getTopGeneSets(gsca, "GSEA.results", "PW_KEGG", allSig=TRUE)
##view a GSEA figure
viewGSEA(gsca, gscName="PW_KEGG", gsName=sigGSs[["PW_KEGG"]][1])
##write html reports
report(object = gsca, experimentName = "GSCATest", species = "Dm", 
allSig = TRUE, keggGSCs = "PW_KEGG", reportDir="GSCATestReport")
##browse the index page
browseURL(file.path(getwd(), "GSCATestReport", "index.html"))
}
}
\keyword{classes}
