\name{InterProGeneSets}
\alias{InterProGeneSets}
\title{
Create gene set collections based on Interpro annotations.
}
\description{
This functions returns a gene set collection where gene sets are Interpro IDs-descriptions. The gene sets are limited to those with interproIDs with genes in the GeneList.
}
\usage{
InterProGeneSets(GeneList,species=c("human","mouse","rat","drosophila","celegans"))
}
\arguments{
  \item{GeneList}{
A numeric vector named with Entrez Gene IDs
}
  \item{species}{
The species for which the data should be read. One of the following: "human", "mouse", "rat", "drosophila", "celegans"
}	
}
\details{
This function uses the package biomaRt, it should therefore be loaded prior to using this function, and this function cannot be used without internet access.
}
\value{
A GeneSetCollection object.
}
\references{
Durinck et al. \emph{BioMart and Bioconductor: a powerful link between biological databases and microarray data analysis.} Bioinformatics 2005 21(16):3439-3440.
}
\author{
Jack C. Rose
}
\examples{
\dontrun{
library(biomaRt)
library(org.Hs.eg.db)
gl<-seq(1:100)
x<-as.list(org.Hs.egENSEMBL)
xx<-names(x)
names(gl)<-sample(xx,size=100)
gsc<-InterProGeneSets(GeneList=gl,species="human")
}
}
