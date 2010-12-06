\name{plotSubNet}
\alias{plotSubNet}
\alias{plotSubNet,NWA-method}
\title{
Plot and save a figure of the enriched subnetwork
}
\description{
This is an generic function. 

When implemented as the S4 method for class \code{\link[HTSanalyzeR:NWA]{NWA}}, this function 
invokes the function \code{\link[HTSanalyzeR:networkPlot]{networkPlot}} to plot and save the subnetwork identified 
by the 'BioNet' package. 

To use this function for objects of class \code{NWA}:

plotSubNet(object, filepath, filename, output, ...)
}
\usage{
plotSubNet(object, ...)
}
\arguments{
	\item{object}{
an object. When implemented as S4 methods of class \code{\link[HTSanalyzeR:NWA]{NWA}}, this argument 
is an object of class \code{\link[HTSanalyzeR:NWA]{NWA}}.
}
	\item{...}{
other arguments. (see below for the arguments supported by the method of
class \code{\link[HTSanalyzeR:NWA]{NWA}})
}
\describe{
	\item{filepath:}{
a single character value specifying the directory where the figure will
be stored
}
	\item{filename:}{
a single character value specifying the name of the figure to save
}
	\item{output:}{
a single character value specifying the format of output image: "pdf"
or "png"
}
	\item{...:}{
all other arguments for the function \code{pdf} or \code{png} such as
the argument 'width' and 'height'
}
}
}
\details{
After the analyses step for an object of class 'NWA', users can generate
the enriched subnetwork identified by the 'BioNet' package. If the slot
'phenotype' has been input when initializing the object, this function
will send it to the function \code{\link[HTSanalyzeR:networkPlot]{networkPlot}} as the argument \code{phenotypeVector}
to highlight nodes in different colors. If the argument \code{species} of the
function \code{\link[HTSanalyzeR:analyze]{analyze}} has been assigned, labels of
nodes of this subnetwork will be mapped to gene symbols corresponding to
the species; otherwise, Entrez identifiers will be used as the labels.    
}

\author{
Xin Wang \email{xw264@cam.ac.uk}
}
\seealso{
\code{networkPlot}, \code{viewSubNet}, \code{analyze}
}
\examples{
\dontrun{
library(BioNet)
##load p-values and phenotypes
data("KcViab_Data4Enrich","KcViab_PVals")
##load Biogrid interactome for Drosophila Melanogaster
data("Biogrid_DM_Interactome")
##create a NWA (NetWork Analysis) object
nwa <- new("NWA", pvalues=KcViab_PVals, phenotypes=KcViab_Data4Enrich, 
interactome=Biogrid_DM_Interactome)
##preprocessing
nwa <- preprocess(nwa, species="Dm", initialIDs="Entrez.gene", 
keepMultipleMappings=TRUE, duplicateRemoverMethod="max")
##To create an interactome:
##nwa<-interactome(nwa, species="Dm", reportDir="HTSanalyzerReport",
##genetic=FALSE)
##do network analysis
nwa <- analyze(nwa, fdr=0.001, species="Dm")
graphics.off()
##plot and save the identified subnetwork
plotSubNet(nwa, filepath=".", filename="subnetwork.pdf", 
output="pdf", width=8, height=8)
}
}











