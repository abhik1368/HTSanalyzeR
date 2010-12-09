\name{report}
\alias{report}
\alias{report,GSCA-method}
\alias{report,NWA-method}
\title{
Write html reports for enrichment or network analyses
}
\description{
This is a generic function.

When implemented as the method of class \code{\link[HTSanalyzeR:GSCA]{GSCA}} or \code{\link[HTSanalyzeR:NWA]{NWA}}, 
this function produces reports for either the Gene Set Collection 
Analysis or the Network Analysis. 

To use \code{report} for objects of class \code{\link[HTSanalyzeR:GSCA]{GSCA}} or \code{\link[HTSanalyzeR:NWA]{NWA}}:

report(object, experimentName="Unknown", species=NULL, ntop=NULL,
allSig=FALSE, keggGSCs=NULL, goGSCs=NULL, reportDir="HTSanalyzerReport")
}
\usage{
report(object, ...)
}
\arguments{
  \item{object}{
an object. When implemented as S4 methods of class \code{\link[HTSanalyzeR:GSCA]{GSCA}} or \code{\link[HTSanalyzeR:NWA]{NWA}}, this 
argument is either an object of class \code{\link[HTSanalyzeR:GSCA]{GSCA}} or \code{\link[HTSanalyzeR:NWA]{NWA}}.
}
  \item{...}{
other arguments. (see below for the arguments supported by the method of
class \code{\link[HTSanalyzeR:GSCA]{GSCA}} and/or \code{\link[HTSanalyzeR:NWA]{NWA}})
}
\describe{
  \item{experimentName:}{
a single character value specifying the name of the experiment (just
for you own record)
} 
  \item{species:}{
a single character value specifying the species for which the data should
be read. Current version support one of the following species: "Dm"
("Drosophila_melanogaster"), "Hs" ("Homo_sapiens"), "Rn" ("Rattus_
norvegicus"), "Mm" ("Mus_musculus"), "Ce" ("Caenorhabditis_elegans").
}
  \item{ntop:}{
a single integer value specifying the number of plots to be produced for
the GSEA analysis. For each gene set collection, plots are produced for
the 'ntop' most significant gene sets.  
}
  \item{allSig:}{
a single logical value determining whether or not to generate plots for
all significant gene sets. A gene set is significant if its corresponding
adjusted p-value is less than the \code{pValueCutoff} set in the function
\code{\link[HTSanalyzeR:analyze]{analyze}}. (see help(analyze) for more details)   
 }
  \item{keggGSCs:}{
a character vector of names of all KEGG gene set collections. This will
help create web links for KEGG terms. 
}
  \item{goGSCs:}{
a character vector of names of all GO gene set collections. This will
help create web links for GO terms. 
}
  \item{reportDir:}{
a single character value specifying the directory to store reports
}
}

}
\details{
This function take in the objects of the two wrapper classes (\code{\link[HTSanalyzeR:GSCA]{GSCA}} 
and \code{\link[HTSanalyzeR:NWA]{NWA}}) and writes a report into the user-specified directory.
An index html file containing a summary of all results and hyperlinked tabs to
more detailed results will be generated in the root directory. The other
html files will be stored in a subdirectory called 'html'. All images
including GSEA plots and the subnetwork figure will be produced in a sub-
directory called 'image'. All documents or text files such as the files
containing significant gene sets of the hypergeometric test results will
be stored in a subdirectory called 'doc'. 
}
\seealso{
\code{\link[HTSanalyzeR:reportAll]{reportAll}}, \code{\link[HTSanalyzeR:writeReportHTSA]{writeReportHTSA}}
}
\author{
Xin Wang, Camille Terfve
}
\examples{
\dontrun{
##(see the vignette for details about the preprocessing of this data set)
library(KEGG.db)
library(GO.db)
library(AnnotationDbi)
##report for a GSCA object
data("KcViab_GSCA")
##append gene set terms to results
KcViab_GSCA<-appendGSTerms(KcViab_GSCA, keggGSCs="PW_KEGG", 
goGSCs=c("GO_BP","GO_MF","GO_CC"))
report(object=KcViab_GSCA, experimentName="KcViab", species="Dm",
allSig=TRUE, keggGSCs="PW_KEGG", goGSCs=c("GO_BP","GO_MF","GO_CC"),
reportDir="HTSanalyzerGSCAReport")
browseURL(file.path(getwd(), "HTSanalyzerGSCAReport", "index.html"))
##report for a NWA object
data("KcViab_NWA")	
report(object=KcViab_NWA, experimentName="KcViab", species="Dm",
allSig=TRUE, keggGSCs="PW_KEGG", goGSCs=c("GO_BP","GO_MF","GO_CC"),
reportDir="HTSanalyzerNWReport")
browseURL(file.path(getwd(), "HTSanalyzerNWReport", "index.html"))
}
}
