\name{cellHTS2OutputStatTests}
\alias{cellHTS2OutputStatTests}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute statistical tests on a cellHTS cellHTSobject.
}
\description{
This function takes a normalized, configured and annotated cellHTS cellHTSobject and computes statistical tests on it for the significance of a set of observations for each condition tested in a hight throughput screen.
}
\usage{
cellHTS2OutputStatTests(cellHTSobject, annotationColumn = "GeneID",
 controls = "neg", alternative = c("two.sided", "less", "greater"),
 logged = FALSE, tests = c("T-test", "MannWhitney", "RankProduct"))
}
%- maybe also 'usage' for other cellHTSobjects documented here.
\arguments{
  \item{cellHTSobject}{
An object of class cellHTS
}
  \item{annotationColumn}{
The name of the column in the fData(cellHTSobject) data frame from which the feature identifiers will be extracted
}
  \item{controls}{
the name of the controls to be used as a control population in the 2 samples tests (this HAS to be corresponding to how these control wellshave been annotated in the column "controlStatus" of the fData(cellHTSobject) dataframe). If nothing is specified, the function will look for negative controls labelled "neg".
}
  \item{alternative}{
two.sided-less-greater (only one of those 3): which is the alternative hypothesis to be tested
}
  \item{logged}{
TRUE-FALSE: has the data been logged during the normalization process
}
  \item{tests}{
"T-test", "MannWhitney", "RankProduct": the tests to be performed, if nothing is specified, all 3 tests will be perfomed. Be aware that the rank product test is slower than the other 2, and returns a percent false discovery (equivalent to a FDR, not a p-value) this one requires the package "RankProduct".
}
}
\details{
The tests are computed taking into account only the wells labelled "sample" in the column "controlStatus" of the fData(cellHTSobject).

The two samples test compare the set of observation for one construct to the values obtained for a population considered as "control".  The one sample tests compare the set of observation for one construct to the median of all values obtained across all constructs labelled as "sample".  This type of test assumes that most constructs are expected to show a negligible effect, it is therefore not advised to used this type of tests when the constructs tested have been pre-screened for being associated with a phenotype.

Please be aware that both types of tests are less reliable when the number of replicates for each construct is low.
}
\value{
A matrix with 2 columns for each type of test (two samples and one sample test) except the Rank Product (no alternative), and a row for each construct (row names corresponding to the identifiers given by the "annotationcolumn" entry).
}
\references{
Michael Boutros, Ligia P. Bras L and Wolfgang Huber. \emph{Analysis of cell-based RNAi screens.} Genome Biology 7:7 R66 (2006)."
}
\author{
Camille Terfve
}
\examples{
\dontrun{
library(cellHTS2)
experimentName<-"KcViab"
dataPath<-system.file(experimentName, package="cellHTS2")
x<-readPlateList("Platelist.txt",name=experimentName,path=dataPath)
x<-configure(x,descripFile="Description.txt",confFile="Plateconf.txt",
logFile="Screenlog.txt",path=dataPath)
xn<-normalizePlates(x,scale="multiplicative",log=FALSE,method="median",
varianceAdjust="none")
xn<-annotate(xn,geneIDFile="GeneIDs_Dm_HFA_1.1.txt",path=dataPath)
test.stats<-cellHTS2OutputStatTests(cellHTSobject=xn,alternative="two.sided",tests="T-test")
}
}
