\name{writeReportHTSA}
\alias{writeReportHTSA}
\title{
Write html reports for enrichment and/or network analyses
}
\description{
This function writes an html report following a complete analysis
of a dataset with the objects of class \code{\link[HTSanalyzeR:GSCA]{GSCA}} and/or \code{\link[HTSanalyzeR:NWA]{NWA}}. 
}
\usage{
writeReportHTSA(gsca, nwa, experimentName="Unknown", species=NULL,
ntop=NULL, allSig=FALSE, keggGSCs=NULL, goGSCs=NULL, reportDir=
"HTSanalyzerReport")
}
\arguments{
  \item{gsca}{
an object of class \code{\link[HTSanalyzeR:GSCA]{GSCA}} 
}
  \item{nwa}{
an object of class \code{\link[HTSanalyzeR:NWA]{NWA}}
}
  \item{experimentName}{
a single character value specifying the name of the experiment (just
for you own record)
} 
  \item{species}{
a single character value specifying the species for which the data should
be read. Current version support one of the following species: "Dm"
("Drosophila_melanogaster"), "Hs" ("Homo_sapiens"), "Rn" ("Rattus_
norvegicus"), "Mm" ("Mus_musculus"), "Ce" ("Caenorhabditis_elegans").
}
  \item{ntop}{
a single integer value specifying the number of plots to be produced for
the GSEA analysis. For each gene set collection, plots are produced for
the 'ntop' most significant gene sets.  
}
  \item{allSig}{
a single logical value determining whether or not to generate plots for
all significant gene sets. A gene set is significant if its corresponding
adjusted p-value is less than the \code{pValueCutoff} set in function
\code{\link[HTSanalyzeR:analyze]{analyze}}.    
 }
  \item{keggGSCs}{
a character vector of names of all KEGG gene set collections. This will
help create web links for KEGG terms. 
}
  \item{goGSCs}{
a character vector of names of all GO gene set collections. This will
help create web links for GO terms. 
}
  \item{reportDir}{
a single character value specifying the directory to store reports
}
}
\details{
This function takes in the objects of the two wrapper classes (\code{\link[HTSanalyzeR:GSCA]{GSCA}}
and \code{\link[HTSanalyzeR:NWA]{NWA}}) and writes a report into a user-specified directory.
An index html file containing a summary of all results and hyperlinked tabs
linking to more detailed results will be generated in the root directory.
The other html files will be stored in a subdirectory called 'html'. All
images including the GSEA plots and subnetwork figure will be produced in
a subdirectory called 'image'. All documents or text files such as the files
containing significant gene sets of the hypergeometric test results will
be stored in a subdirectory called 'doc'. 
}
\author{
Xin Wang, Camille Terfve
}
\seealso{
\code{\link[HTSanalyzeR:report]{report}}, \code{\link[HTSanalyzeR:reportAll]{reportAll}}
}
\examples{
\dontrun{
##(see the vignette for details about the preprocessing of this data set)
data("KcViab_GSCA")
data("KcViab_NWA")	
##report both analyses
writeReportHTSA(gsca=KcViab_GSCA, nwa=KcViab_NWA, experimentName="KcViab",
species="Dm", allSig=TRUE, keggGSCs="PW_KEGG", goGSCs=c("GO_BP","GO_MF",
"GO_CC"), reportDir="HTSanalyzerReport")
browseURL(file.path(getwd(), "HTSanalyzerReport", "index.html"))
}
}
