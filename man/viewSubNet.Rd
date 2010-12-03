\name{viewSubNet}
\alias{viewSubNet}
\alias{viewSubNet,NWA-method}
\title{
Plot to screen a figure of the enriched subnetwork
}
\description{
This is an generic function. 

When implemented as the S4 method for class 'NWA', function \code{viewSubNet}
invokes function 'networkPlot' to plot to screen the subnetwork identified 
by 'BioNet'.

}
\usage{
viewSubNet(object)
}
\arguments{
	\item{object}{
an object. When implemented as S4 methods of class 'NWA', this argument 
is an object of class 'GSCA'.
}
}
\details{
After the 'analyze' step for an object of class 'NWA', users can generate
the enriched subnetwork identified by package 'BioNet'. If slot 'phenotype'
has been input when initializing the object, this function will send it to
function 'networkPlot' as the argument 'phenotypeVector' to highlight nodes
in different colors. If the argument 'species' of function 'analyze' has
been assigned, labels of nodes of this subnetwork will be mapped to gene
symbols corresponding to the species; otherwise, Entrez identifiers will
be used as these labels.    
}

\author{
Xin Wang \email{xw264@cam.ac.uk}
}
\seealso{
\code{networkPlot}, \code{analyze}, \code{NWA}
}
\examples{
\dontrun{
library(BioNet)
##load p-values and phenotypes
data("pvalues","data4enrich")
##load Biogrid interactome for Drosophila Melanogaster
data("Biogrid.Dm.interactome")
##create a NWA (NetWork Analysis) object
nwa <- new("NWA", pvalues=pvalues, phenotypes=data4enrich, 
interactome=Biogrid.Dm.interactome)
##preprocessing
nwa <- preprocess(nwa, species="Dm", initialIDs="Entrez.gene", 
keepMultipleMappings=TRUE, duplicateRemoverMethod="max")
##To create an interactome:
##nwa<-interactome(nwa, species="Dm", reportDir="HTSanalyzerReport",
##genetic=FALSE)
##do network analysis
nwa <- analyze(nwa, fdr=0.001, species="Dm")
graphics.off()
##plot to screen the identified subnetwork
viewSubNet(nwa)
}
}











