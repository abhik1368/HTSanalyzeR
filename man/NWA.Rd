\name{NWA-class}
\alias{NWA}
\alias{NWA-class}
\docType{class}
\title{
A class for NetWork Analyses
}
\description{
This class includes methods to do network analysis on a vector of p-values 
from experiments such as RNA interference screening. 
}
\section{Objects from the Class}{ 
Objects of class \code{NWA} can be created from \code{
new("NWA", pvalues, phenotypes=NULL, interactome=NULL)} (see
the examples below)
}
\section{Slots}{
	\describe{	
		\item{\code{pvalues}:}{
		a numeric vector of p-values.
		}
		
		\item{\code{phenotypes}:}{
		a numeric or integer vector of phenotypes.
		}
		
		\item{\code{interactome}:}{
		an object of class \code{graphNEL}.
		}
		
		\item{\code{fdr}:}{
		one parameter for BioNet to score nodes in the interactome.
		}
		
		\item{\code{result}:}{
		a list of subnetwork module identified by 'BioNet' and a vector of 
		labels for nodes of the subnetwork module.
		}
		
		\item{\code{summary}:}{
		a list of summary information for \code{p-values}, \code{phenotypes}, 
		\code{interactome} and \code{result}.
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
		a) removing NAs in \code{p-values} and \code{phenotypes}; 
		b) invoking function \code{duplicateRemover} to process
		duplicated phenotypes and p-values (see \code{duplicateRemover} for more 
		details); 
		c) invoking function \code{annotationConvertor} to convert annotations; 
		}
		\item{\code{analyze}}{
		invoke function \code{networkAnalysis} to identify sub-networks based on
		input parameter list \code{para}.
		}
		\item{\code{summary}}{
		print summary information of \code{p-values}, \code{phenotypes}, \code{interactome} and \code{result}.
		}
		\item{\code{interactome}}{
		build an interactome for the network analysis.
		}
		\item{\code{plotSubNet}}{
		plot figures for identified subnetworks.
		}
		\item{\code{report}}{
		generate html reports.
		}
	}
}
\author{Xin Wang \email{xw264@cam.ac.uk}}
\seealso{
  \code{\link[NWA:preprocess]{preprocess}}
  \code{\link[NWA:analyze]{analyze}}
  \code{\link[NWA:summary]{summary}}
  \code{\link[NWA:interactome]{interactome}}
  \code{\link[NWA:plotSubNet]{plotSubNet}}
  \code{\link[NWA:report]{report}}
}
\examples{
##load p-values and phenotypes
data("pvalues","data4enrich")
##load Biogrid interactome for Drosophila Melanogaster
data("Biogrid.Dm.interactome")
##create a NWA (NetWork Analysis) object
nwa<-new("NWA",pvalues=pvalues,phenotypes=data4enrich, 
interactome=Biogrid.Dm.interactome)
##preprocessing
nwa<-preprocess(nwa,species="Dm",initialIDs="Entrez.gene", 
keepMultipleMappings=TRUE, duplicateRemoverMethod="max")
##To create an interactome:
##nwa<-interactome(nwa,species="Dm",reportDir="HTSanalyzerReport",genetic=FALSE)
##do network analysis
nwa<-analyze(nwa, fdr=0.001, species="Dm")
##plot identified subnetwork
plotSubNet(nwa, filepath=".", filename="subnetwork.png")
##report as html pages
report(object=nwa, experimentName="NWATest", species="Dm", 
allSig=TRUE, keggGSCs="PW.KEGG", goGSCs=c("GO.BP","GO.MF","GO.CC"), 
reportDir="NWATestReport")
}
\keyword{classes}
