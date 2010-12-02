\name{cellHTS2OutputStatTests}
\alias{cellHTS2OutputStatTests}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute statistical tests on a cellHTS object
}
\description{
This function takes a normalized, configured and annotated cellHTS
cellHTSobject and computes statistical tests on it for the significance
of a set of observations for each condition tested in a high-throughput
screen.
}
\usage{
cellHTS2OutputStatTests(cellHTSobject, annotationColumn = "GeneID",
controls = "neg", alternative = "two.sided", logged = FALSE, tests =
"T-test")
}
%- maybe also 'usage' for other cellHTSobjects documented here.
\arguments{
  \item{cellHTSobject}{
an object of class cellHTS
}
  \item{annotationColumn}{
a single character value specifying the name of the column in the
fData(cellHTSobject) data frame from which the feature identifiers will
be extracted
}
  \item{controls}{
a single character value specifying the name of the controls to be used 
as a control population in the 2 samples tests (this HAS to be 
corresponding to how these control wellshave been annotated in the 
column "controlStatus" of the fData(cellHTSobject) dataframe). If nothing 
is specified, the function will look for negative controls labelled "neg".
}
  \item{alternative}{
a single character value specifying the alternative hypothesis:"two.sided",
"less" or "greater"
}
  \item{logged}{
a single logical value specifying whether or not the data has been logged
during the normalization process
}
  \item{tests}{
a single character value specifying the tests to be performed: "T-test",
"MannWhitney" or "RankProduct". If nothing is specified, all 3 tests will
be performed. Be aware that the rank product test is slower than the other
2, and returns a percent false discovery (equivalent to a FDR, not a
p-value) this one requires the package "RankProduct".
}
}
\details{
The tests are computed taking into account only the wells labelled "sample"
in the column "controlStatus" of the fData(cellHTSobject).

The two samples test compare the set of observation for one construct to
the values obtained for a population considered as "control".  The one
sample tests compare the set of observation for one construct to the median
of all values obtained across all constructs labelled as "sample".  This
type of test assumes that most constructs are expected to show a negligible
effect, it is therefore not advised to used this type of tests when the
constructs tested have been pre-screened for being associated with a
phenotype.

Please be aware that both types of tests are less reliable when the number
of replicates for each construct is low.
}
\value{
a matrix with 2 columns for each type of test (two samples and one sample
test) except the Rank Product (no alternative), and a row for each construct
(row names corresponding to the identifiers given by the "annotationcolumn"
entry).
}
\references{
Michael Boutros, Ligia P. Bras L and Wolfgang Huber. \emph{Analysis of
cell-based RNAi screens.} Genome Biology 7:7 R66 (2006)."
}
\author{
Camille Terfve, Xin Wang
}
\examples{
\dontrun{
library(cellHTS2)
##load normalized cellHTS object (see the vignette for details about the
##preprocessing of this data set)
data("KcViab.normalized")
test.stats <- cellHTS2OutputStatTests(cellHTSobject=KcViab.normalized,
alternative="two.sided", tests="T-test")
}
}
