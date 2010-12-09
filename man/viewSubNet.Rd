\name{viewSubNet}
\alias{viewSubNet}
\alias{viewSubNet,NWA-method}
\title{
Plot a figure of the enriched subnetwork
}
\description{
This is an generic function. 

When implemented as the S4 method of class \code{\link[HTSanalyzeR:NWA]{NWA}}, this function 
invokes the function \code{\link[HTSanalyzeR:networkPlot]{networkPlot}} to plot the subnetwork identified 
by the 'BioNet' package.

}
\usage{
viewSubNet(object)
}
\arguments{
	\item{object}{
an object. When implemented as S4 methods of class \code{\link[HTSanalyzeR:NWA]{NWA}}, this argument 
is an object of class 'GSCA'.
}
}
\details{
After the analyses step for an object of class \code{\link[HTSanalyzeR:NWA]{NWA}}, the user can generate
the enriched subnetwork identified by the package 'BioNet'. If the slot
'phenotype' was inputted during the initialization of the object, this 
function will send it to the function \code{\link[HTSanalyzeR:networkPlot]{networkPlot}} as the argument 
'phenotypeVector' to highlight nodes in different colors. If the argument 
'species' of the function \code{\link[HTSanalyzeR:analyze]{analyze}} has been assigned,
the labels of the nodes of this subnetwork will be mapped to gene symbols
corresponding to the species; otherwise, the Entrez identifiers will be
used as the labels.    
}

\author{
Xin Wang \email{xw264@cam.ac.uk}
}
\seealso{
\code{networkPlot}, \code{analyze}, \code{plotSubNet}
}
\examples{
\dontrun{
##load sample data
data("KcViab_NWA")
##plot to screen the identified subnetwork
viewSubNet(KcViab_NWA)
}
}











