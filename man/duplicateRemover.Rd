\name{duplicateRemover}
\alias{duplicateRemover}
\title{
Remove the duplicates in a gene list.
}
\description{
This function gets rid of the duplicates in a gene list. It is used to prepare a gene list for enrichment analysis.  It also orders the gene list by phenotype, as needed for the gene set enrichment analysis functions.
}
\usage{
duplicateRemover(geneList, method = c("max", "min", "average", 
"fold.change.average"), absValue = FALSE)
}
\arguments{
  \item{geneList}{
A gene list, as a single named vector
}
  \item{method}{
A choice of method to remove the duplicates (should the minimum, maximum of average observation for a same construct be kept). The minimum and maximum should be understood in terms of absolute values (i.e. min/max effect, no matter the sign). The fold change average method converts the fold changes to ratios, averages them and converts the average back to a fold change.
}
  \item{absValue}{
A choice of whether the values should be converted to absolute value and then ordered (if TRUE),or ordered as they are (if FALSE).
}
}
\value{
An ordered gene list (increasing) with duplicates removed, as a single named vector
If absValue is TRUE:a gene list (absolute values) ordered by absolute value (increasing) with duplicates removed, as a single named vector.
}
\author{
Camille Terfve and John C. Rose
}
\examples{
x<-c(5,1,3,-2,6)
names(x)<-c("gene1","gene3","gene7","gene3","gene4")
xprocessed<-duplicateRemover(geneList=x,method="max",absValue=FALSE)
}
