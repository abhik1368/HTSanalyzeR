\name{duplicateRemover}
\alias{duplicateRemover}
\title{
Remove the duplicates in a gene list
}
\description{
This function gets rid of the duplicates in a gene list. It is used to 
prepare a gene list for enrichment analysis.
}
\usage{
duplicateRemover(geneList, method = "max")
}
\arguments{
  \item{geneList}{
a single named numeric or integer vector
}
  \item{method}{
a single character value specifying the method to remove the duplicates 
(should the minimum, maximum or average observation for a same construct 
be kept). Current version provides "min" (minimum), "max" (maximum), 
"average" and "fc.avg" (fold change average). The minimum and maximum 
should be understood in terms of absolute values (i.e. min/max effect, 
no matter the sign). The fold change average method converts the fold 
changes to ratios, averages them and converts the average back to a 
fold change.
}
}
\value{
a gene list (increasing) with duplicates removed, as a single named vector
}
\author{
Camille Terfve, John C. Rose and Xin Wang
}
\examples{
x<-c(5,1,3,-2,6)
names(x)<-c("gene1","gene3","gene7","gene3","gene4")
xprocessed<-duplicateRemover(geneList=x,method="max")
}
