\name{gseaScores}
\alias{gseaScores}
\alias{gseaScoresBatch}
\alias{gseaScoresBatchParallel}
\title{
Compute enrichment scores for GSEA (Gene Set Enrichment Analysis)
}
\description{
'gseaScores' computes the enrichment score, running sum scores and
positions of hits for GSEA on one gene set.

'gseaScoresBatch' computes enrichment scores for both input 'geneList'
and its permutations for GSEA on one gene set.

'gseaScoresBatchParallel' computes enrichment scores for both input
'geneList' and their permutations for GSEA on multiple gene sets in 
parallel.
}
\usage{
gseaScores(geneList, geneSet, exponent=1, mode = "score")

gseaScoresBatch(geneList, geneNames.perm, geneSet, exponent=1,
nPermutations=1000)

gseaScoresBatchParallel(geneList, geneNames.perm, collectionOfGeneSets,
exponent=1, nPermutations=1000)
}
\arguments{
  \item{geneList}{
a numeric or integer vector of phenotypes in descending or ascending order
with elements named by their gene identifiers (no duplicates nor NA values)
}
  \item{geneNames.perm}{
a character matrix including permuted gene identifiers of input 'geneList'.
The first column of this matrix should be the gene identifiers of 'geneList',
while other columns are vectors of permuted gene identifiers.   
}
  \item{geneSet}{
a character vector specifying a gene set (no names, just a vector of
characters corresponding to the IDs)
}
  \item{collectionOfGeneSets}{
a list of gene sets. Each gene set in the list is a character vector of gene
identifiers
}
  \item{exponent}{
a single integer or numeric value used in weighting phenotypes in GSEA
}
  \item{nPermutations}{
a single integer or numeric value specifying the number of permutations 
for deriving p-values in GSEA
}
  \item{mode}{
a single character value specifying to return only a score (if set as
"score"), or all the necessary elements to make a plot (if "graph")
(see "gseaPlot" for more details).
}
}
\details{
The type of identifiers used in the gene sets and gene list must obviously
match.
}
\value{
* gseaScores will return:
  \item{enrichmentScore}{
a single numeric value of the enrichment score
  }
  \item{runningScore}{
a numeric vector of running sum scores (only in mode "graph")
  }
  \item{positions}{
a numeric vector of positions in the ranked phenotype vector of the genes 
in the gene set (only in mode "graph")
}

* gseaScoresBatch will return a list consisting of:
  \item{scoresObserved}{
a single numeric value of the enrichment score for input 'geneList'
}
  \item{scoresperm}{
a numeric vector of enrichment scores for permutation tests
}

* gseaScoresBatchParallel will return a matrix, in which one column is
the 'gseaScoresBatch' result for each gene set.
}
\references{
Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L.,
Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander,
E. S. & Mesirov, J. P. (2005)  \emph{Gene set enrichment analysis: A
knowledge-based approach for interpreting genome-wide expression profiles.}
Proc. Natl. Acad. Sci. USA 102, 15545-15550.
}
\author{
Xin Wang, Camille Terfve
}
\examples{
##example 1
gl <- runif(100, min=0, max=5)
gl <- gl[order(gl, decreasing=TRUE)]
names(gl) <- as.character(sample(x=seq(from=1, to=100, by=1), size=100,
replace=FALSE))
gs <- sample(names(gl), size=20, replace=FALSE)
##GSEA for input phenotype vector on a single gene set
gsea <- gseaScores(geneList=gl, geneSet=gs, mode="score", exponent=1)
##GSEA for both input phenotype vector and permutation tests on a single gene set
nPermutations <- 100
glPerm <- sapply(1:nPermutations, function(n) names(gl)[sample(1:length(gl),
length(gl), replace=FALSE)])
glPerm <- cbind(names(gl), glPerm)
gseaBatch<-gseaScoresBatch(geneList=gl, geneNames.perm=glPerm, geneSet=gs,
nPermutations=100, exponent=1)
##example 2
\dontrun{
library(KEGG.db)
library(org.Dm.eg.db)
library(snow)
##load phenotype vector (see the vignette for details about the
##preprocessing of this data set)
data("KcViab_Data4Enrich")
DM_KEGG <- KeggGeneSets(species="Dm")
##GSEA for input gene list on a single gene set
test <- gseaScores(geneList=KcViab_Data4Enrich, geneSet=DM_KEGG[[1]], 
exponent=1, mode="graph")
##GSEA for both input gene list and permutation tests on multiple gene
##sets in parallel
nPermutations <- 100
glPerm <- sapply(1:nPermutations, function(n) names(KcViab_Data4Enrich)[
sample(1:length(KcViab_Data4Enrich), length(KcViab_Data4Enrich), 
replace=FALSE)])
glPerm<-cbind(names(KcViab_Data4Enrich), glPerm)
options(cluster=makeCluster(4,"SOCK"))
gseaBatchPar <- gseaScoresBatchParallel(geneList=KcViab_Data4Enrich,
geneNames.perm=glPerm, collectionOfGeneSets=DM_KEGG[1:10], nPermutations=100,
exponent=1)
if(is(getOption("cluster"),"cluster")) {
	stopCluster(getOption("cluster"))
	options(cluster=NULL)
}
}
}
