\name{pairwiseGseaPlot}
\alias{pairwiseGseaPlot}
\title{
Produce a collection of plots for GSEA, by pairs of phenotypes.
}
\description{
Function that produces plots from the GSEA analysis, but plotting in parallel for two phenotypes, for one gene set.
}
\usage{
pairwiseGseaPlot(gl1,gl2,geneSet,exponent=1,output="png",filepath,filename)
}
\arguments{
  \item{gl1}{
a named vector where names are gene identifiers of the same type as the ones in the gene set collection, and values are the measurement on phenotype 2 corresponding to those genes. This vector MUST be ordered
}
  \item{gl2}{
a named vector where names are gene identifiers of the same type as the ones in the gene set collection, and values are the measurement on phenotype 2 corresponding to those genes this vector MUST be ordered
}
  \item{geneSet}{
a vector of identifiers (see the help on function gseaScores)
}
  \item{exponent}{
exponent for the GSEA, keep it at one (see the help on function gseaScores)
}
  \item{output}{
either "png" or "pdf"
}
  \item{filepath}{
The directory where these plots will be stored
}
  \item{filename}{
Name of the gene set for which these plots are produced
}
}
\details{
The plot is going to be produced in the current working directory.
}
\value{
This function does not return anything, it just produces an image file containing the plots.
}
\references{
Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S. & Mesirov, J. P. (2005)  \emph{Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles.} Proc. Natl. Acad. Sci. USA 102, 15545-15550.
}
\author{
Camille Terfve
}
\examples{
\dontrun{
gl1<-runif(100,min=-5,max=5)
gl1<-gl1[order(gl1,decreasing=TRUE)]
gl2<-runif(100,min=-5,max=5)
gl2<-gl2[order(gl2,decreasing=TRUE)]
names(gl1)<-as.character(sample(x=seq(from=1,to=100,by=1),size=100,replace=FALSE))
names(gl2)<-as.character(sample(x=seq(from=1,to=100,by=1),size=100,replace=FALSE))
gs1<-sample(names(gl1),size=20,replace=FALSE)
gs2<-sample(names(gl1),size=20,replace=FALSE)
gsc<-list(subset1=gs1,subset2=gs2)
pairwiseGseaPlot(gl1=gl1,gl2=gl2,geneSet=gsc[["subset1"]],filepath=".", filename="geneSet")
}
}
