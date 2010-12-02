\name{networkAnalysis}
\alias{networkAnalysis}
\title{
Identify enriched subnetworks 
}
\description{
This function finds subnetworks enriched for genes with significant
phenotypes based on the package 'BioNet'. 
}
\usage{
networkAnalysis(pvalues,graph,fdr=0.001,verbose=TRUE)
}
\arguments{
  \item{pvalues}{
a numeric vector of p-values   
 } 
  \item{fdr}{
a single numeric value specifying the false discovery for the scoring of
nodes (see BioNet::scoreNodes and Dittrich et al., 2008 for details)
}
  \item{verbose}{
a single logical value deciding to display detailed messages (when
verbose=TRUE) or not (when verbose=FALSE)
}
}
\details{
This function takes in a vector of p-values and a graph standing for the
interactome to identify the maximum scoring subnetwork based on the BioNet
package.
}
\value{
a subnetwork module of class 'graphNEL' 
}
\references{
Beisser D, Klau GW, Dandekar T, Muller T, Dittrich MT. BioNet: an R-Package
for the functional analysis of biological networks. Bioinformatics. 2010
Apr 15;26(8):1129-30.

Dittrich MT, Klau GW, Rosenwald A., Dandekar T and Muller T. \emph{Identifying
functional modules in protein-protein interaction networks: an integrated
exact approach.} Bioinformatics 2008 24(13):i223-i231.
}
\author{
Camille Terfve, Xin Wang
}
\seealso{
networkPlot, NWA-class
}
\examples{
\dontrun{
library(BioNet)
##load pvalues (see the vignette for details about the preprocessing of
##this data set)
data("pvalues")
##load interactome
data("Biogrid.Dm.interactome")
##identify subnetworks
enrichedSubNet <- networkAnalysis(pvalues=pvalues,
graph=Biogrid.Dm.interactome, fdr=0.001, verbose=TRUE)	
}
}
