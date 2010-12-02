\name{plotSubNet}
\alias{plotSubNet}
\alias{viewSubNet}
\title{
Plot and/or save a figures of the enriched subnetwork
}
\description{
viewSubNet:

an S4 method for objects of class 'NWA', which invokes function 'networkPlot' 
to plot to screen the subnetwork identified by 'BioNet'.

plotSubNet:

an S4 method for objects of class 'NWA', which invokes function 'networkPlot' 
to plot and save the subnetwork identified by 'BioNet'. 
}
\usage{

viewSubNet(object)

plotSubNet(object, filepath, filename, output, ...)
}
\arguments{
  \item{object}{
an object of class 'NWA'
}
  \item{filepath}{
a single character value specifying the directory where the figure will
be stored
}
  \item{filename}{
a single character value specifying the name of the figure to save
}
  \item{output}{
a single character value specifying the format of output image: "pdf"
or "png"
}
  \item{...}{
all other arguments for 'pdf' or 'png' such as 'width' and 'height'
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
##plot identified subnetwork
viewSubNet(nwa)
plotSubNet(nwa, filepath=".", filename="subnetwork.pdf", 
output="pdf", width=8, height=8)
}
}











