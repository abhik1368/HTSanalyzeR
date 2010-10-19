\name{pairwisePhenoMannWith}
\alias{pairwisePhenoMannWith}
\title{
Performs a Mann-Whitney test for shift in location of genes from gene sets, on a pair of phenotypes.
}
\description{
This function performs a Mann-Whitney test for shift in location of genes from gene sets, on a pair of phenotypes. It looks for gene sets that are represented towards the 2 different ends of two ranked lists of genes, i.e. whose phenotype distribution is located around two different values in the two phenotypes list, rather than spread on the whole list in both lists.
}
\usage{
pairwisePhenoMannWith(gl1,gl2,gsc,cutoff=15,MHT="BH")
}
\arguments{
  \item{gl1}{
a named vector where names are gene identifiers of the same type as the ones in the gene set collection, and values are the measurement on phenotype 1 corresponding to those genes
}
  \item{gl2}{
a named vector where names are gene identifiers of the same type as the ones in the gene set collection, and values are the measurement on phenotype 2 corresponding to those genes
}
  \item{gsc}{
a list of gene sets 
}
  \item{cutoff}{
minimum number of genes in one gene set (see the help on function collectionGsea)
}
  \item{MHT}{
a method for multiple hypothesis testing correction (one of the p.adjust methods)
}
}
\value{
A table with a row for each gene set, containing the p value for the Mann-Whitney test, and the adjusted
p value.  The table is ordered by the p value column.
}
\author{
Camille Terfve
}
\examples{
gl1<-runif(100,min=-5,max=5)
gl2<-runif(100,min=-5,max=5)
names(gl1)<-as.character(sample(x=seq(from=1,to=100,by=1),size=100,replace=FALSE))
names(gl2)<-as.character(sample(x=seq(from=1,to=100,by=1),size=100,replace=FALSE))
gs1<-sample(names(gl1),size=20,replace=FALSE)
gs2<-sample(names(gl1),size=20,replace=FALSE)
gsc<-list(subset1=gs1,subset2=gs2)
gl1gl2.gsc<-pairwisePhenoMannWith(gl1=gl1,gl2=gl2,gsc=gsc,cutoff=2)
}
