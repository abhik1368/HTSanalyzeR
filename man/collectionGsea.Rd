\name{collectionGsea}
\alias{collectionGsea}
\title{
Compute observed and permutation-based scores for a collection of Gene Sets.
}
\description{
This function computes observed and permutation-based scores associated with a gene set enrichment analysis for a collection of Gene Sets.
}
\usage{
collectionGsea(collectionOfGeneSets, geneList, exponent, permutationsMatrix, cutoffGeneSetSize)
}
\arguments{
  \item{collectionOfGeneSets}{
Collection of gene sets in "GeneSetCollection" format
}
  \item{geneList}{
geneList is a named data vector, ordered, as for gseaScores.  The data vector cannot contain any duplicates nor NAs.
}
  \item{exponent}{
The exponent of the GSEA method.
}
  \item{permutationsMatrix}{
A matrix of permutationsMatrix from "permutationsMatrixGenerator".
}
  \item{cutoffGeneSetSize}{
Minimum size required for a gene set to be considered.
}
}
\value{
\item{Observed.scores}{The observed score for the given gene set (a named vector)}
\item{Permutation.scores}{The scores for the permutationsMatrix (one column for each permutation and a row for each gene set)}
}
\references{
Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S. & Mesirov, J. P. (2005) \emph{ Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles.} Proc. Natl. Acad. Sci. USA 102, 15545-15550.
}
\author{
Camille Terfve
}
\examples{
gl<-runif(100,min=0,max=5)
gl<-order(gl)
names(gl)<-as.character(sample(x=seq(from=1,to=100,by=1),size=100,replace=FALSE))
gs1<-GeneSet(sample(names(gl),size=20,replace=FALSE), setName="subset1")
gs2<-GeneSet(sample(names(gl),size=20,replace=FALSE), setName="subset2")
gsc<-GeneSetCollection(c(gs1,gs2))
perm<-permutationsGeneratorGsea(nPermutations=1000,geneList=gl)
scores.collection<-collectionGsea(collectionOfGeneSets=gsc,geneList=gl,
exponent=1,permutationsMatrix=perm,cutoffGeneSetSize=5)
pvalues.collection<-permutationPvalueCollectionGsea(
permScores=scores.collection$Permutation.scores,
dataScores=scores.collection$Observed.scores)
#
\dontrun{
library(KEGG.db)
library(org.Dm.eg.db)
library(cellHTS2)
experimentName<-"KcViab"
dataPath<-system.file(experimentName, package="cellHTS2")
x<-readPlateList("Platelist.txt",name=experimentName,path=dataPath)
x<-configure(x,descripFile="Description.txt",confFile="Plateconf.txt",
logFile="Screenlog.txt",path=dataPath)
xn<-normalizePlates(x,scale="multiplicative",log=FALSE,method="median",
varianceAdjust="none")
xn<-annotate(xn,geneIDFile="GeneIDs_Dm_HFA_1.1.txt",path=dataPath)
xsc<-scoreReplicates(xn,sign="-",method="zscore")
xsc<-summarizeReplicates(xsc,summary="mean")
data4enrich<-as.vector(Data(xsc));
names(data4enrich)<-fData(xsc)[,"GeneID"] 
data4enrich<-data4enrich[-which(is.na(names(data4enrich))==TRUE)]
data4enrichentrez<-drosoAnnotationConvertor(geneList=data4enrich,
initialIDs="FlybaseCG")
data4enrichentrez<-data4enrichentrez[-which(is.na(data4enrichentrez))]
data4enrichentrez<-duplicateRemover(geneList=data4enrichentrez,method="max")
kegg.droso<-KeggGeneSets(species="Drosophila_melanogaster")
perm<-permutationsGeneratorGsea(nPermutations=1000,
geneList=data4enrichentrez)
scores.collection<-collectionGsea(collectionOfGeneSets=kegg.droso,
geneList=data4enrichentrez,exponent=1,permutationsMatrix=perm,cutoffGeneSetSize=15)
pvalues.collection<-permutationPvalueCollectionGsea(
permScores=scores.collection$Permutation.scores,
dataScores=scores.collection$Observed.scores)
}
}
