\name{pairwiseGsea}
\alias{pairwiseGsea}
\title{
Performs GSEA on a pair of phenotypes.
}
\description{
This function performs pairwise GSEA: it looks for gene sets that are specifically over-represented towards the 2 different ends of two ranked lists of genes, in a gene set collection.
}
\usage{
pairwiseGsea(gl1,gl2,gsc,exponent=1,nperm=1000,cutoff=15,MHT="BH")
}
\arguments{
  \item{gl1}{
a named vector where names are gene identifiers of the same type as the ones in the gene set collection, and values are the measurement on phenotype 1 corresponding to those genes this vector MUST be ordered (decreasing of increasing)
}
  \item{gl2}{
a named vector where names are gene identifiers of the same type as the ones in the gene set collection, and values are the measurement on phenotype 2 corresponding to those genes this vector MUST be ordered
}
  \item{gsc}{
a gene set collection in geneSetCollection format (see the help on function collectionGsea)
}
  \item{exponent}{
exponent for the GSEA, keep it at one (see the help on function collectionGsea)
}
  \item{nperm}{
number of permutations for the GSEA (see the help on function collectionGsea)
}
  \item{cutoff}{
minimum number of genes in one gene set (see the help on function collectionGsea)
}
  \item{MHT}{
a method for multiple hypothesis testing correction (one of the p.adjust methods)
}
}
\value{
A table with a row for each gene set, containing the p values for the GSEA, and the observed scores for each of the phenotypes independently.  The table is ordered by the p value column.
}
\details{
phenotypes 1 and 2 must be measured on the same genes, i.e. the two vectors  gl1 and gl2 must have the same length and their names must match, but the two vectors must be ordered separately, i.e.one phenotype vector is ordered based on the values of that phenotype only
}
\references{
Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S. & Mesirov, J. P. (2005) \emph{ Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles.} Proc. Natl. Acad. Sci. USA 102, 15545-15550.
}
\author{
Camille Terfve
}
\examples{
gl1<-runif(100,min=-5,max=5)
gl1<-gl1[order(gl1,decreasing=TRUE)]
gl2<-runif(100,min=-5,max=5)
gl2<-gl2[order(gl2,decreasing=TRUE)]
names(gl1)<-as.character(sample(x=seq(from=1,to=100,by=1),size=100,replace=FALSE))
names(gl2)<-as.character(sample(x=seq(from=1,to=100,by=1),size=100,replace=FALSE))
gs1<-sample(names(gl1),size=20,replace=FALSE)
gs2<-sample(names(gl1),size=20,replace=FALSE)
gsc<-list(subset1=gs1,subset2=gs2)
gl1gl2.gsc<-pairwiseGsea(gl1=gl1,gl2=gl2,gsc=gsc)
}
