\name{permutationPvalueCollectionGsea}
\alias{permutationPvalueCollectionGsea}
\title{
This function computes the GSEA p-values for a collection of gene sets
}
\description{
Compute the nominal p-value associated with a GSEA for a collection of gene sets, from the outputs of collectionGsea
}
\usage{
permutationPvalueCollectionGsea(permScores, dataScores)
}
\arguments{
  \item{permScores}{
The output of collectionGsea: matrix of permutation-based scores (rows=gene sets, named)
}
  \item{dataScores}{
The output of collectionGsea: vector of observed scores (named)
}
}
\value{
A vector of p-values (named), one for each gene set
}
\references{
Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S. & Mesirov, J. P. (2005) \emph{ Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles.} Proc. Natl. Acad. Sci. USA 102, 15545-15550.
}
\author{
Camille Terfve
}
\seealso{
collectionGsea
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
data4enrich<-data4enrich[which(!is.na(names(data4enrich))==TRUE)]
data4enrichentrez<-drosoAnnotationConvertor(geneList=data4enrich,
initialIDs="FlybaseCG")
data4enrichentrez<-data4enrichentrez[which(!is.na(data4enrichentrez))]
data4enrichentrez<-duplicateRemover(geneList=data4enrichentrez,method="max")
kegg.droso<-KeggGeneSets(species="Dm")
scores.collection<-collectionGsea(collectionOfGeneSets=kegg.droso,
geneList=data4enrichentrez,exponent=1,npermutations=1000,cutoffGeneSetSize=15)
pvalues.collection<-permutationPvalueCollectionGsea(
permScores=scores.collection$Permutation.scores,
dataScores=scores.collection$Observed.scores)
}
}
