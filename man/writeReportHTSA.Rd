\name{writeReportHTSA}
\alias{writeReportHTSA}
\alias{reportAll}
\alias{report}
\title{
Write html reports for enrichment and/or network analyses
}
\description{
\code{writeReportHTSA} writes an html report following a complete analysis of a dataset with the objects of class 'GSCA' and 'NWA'. 

\code{reportAll} is the S4 method version of \code{writeReportHTSA}.

\code{report} is designed specifically for each class. It can take either an object of \code{GSCA} or \code{NWA}, and produce reports for either the Gene Set Collection Analysis or the Network Analysis. 
}
\usage{
writeReportHTSA(gsca, nwa, experimentName="Unknown", species=NULL, ntop=NULL, allSig=FALSE, keggGSCs=NULL, goGSCs=NULL, reportDir="HTSanalyzerReport")

reportAll(gsca, nwa, experimentName="Unknown", species=NULL, ntop=NULL, allSig=FALSE, keggGSCs=NULL, goGSCs=NULL, reportDir="HTSanalyzerReport")

report(object, experimentName="Unknown", species=NULL, ntop=NULL, allSig=FALSE, keggGSCs=NULL, goGSCs=NULL, reportDir="HTSanalyzerReport")
}
\arguments{
  \item{object}{
either an object of class \code{GSCA} or \code{NWA}
}
  \item{gsca}{
an object of class \code{GSCA} (see ?GSCA for more details)
}
  \item{nwa}{
an object of class \code{NWA} (see ?NWA for more details)
}
  \item{experimentName}{
a single character value specifying the name of the experiment (just for you own record)
} 
  \item{species}{
a single character value specifying the species for which the data should be read. Current version support one of the following species: "Dm" ("Drosophila_melanogaster"), "Hs" ("Homo_sapiens"), "Rn" ("Rattus_norvegicus"), "Mm" ("Mus_musculus"), "Ce" ("Caenorhabditis_elegans").
}
  \item{ntop}{
a single integer value specifying the number of plots to be produced for the GSEA analysis.  For each gene set collection, plots are produced for the 'ntop' most significant p values.  
}
  \item{allSig}{
a single logical value determining whether or not to generate plots for all significant gene sets. A gene set is significant if its corresponding adjusted p-value is less than the \code{pValueCutoff} set in function \code{analyze}. (see function \code{analyze} for more details)   
 }
  \item{keggGSCs}{
a character vector of names of all KEGG gene set collections. This will help create web links for KEGG terms. 
}
  \item{goGSCs}{
a character vector of names of all GO gene set collections. This will help create web links for GO terms. 
}
  \item{reportDir}{
a single character value specifying the directory to store reports
}
}
\details{
These functions take the objects of the two wrapper classes (\code{GSCA} and \code{NWA}) and write a report in user-specified directory. An index html file containing a summary of all results and hyperlink tables to more detailed results will be generated in the root directory. The other html files will be stored in a subdirectory called 'html'. All images including GSEA plots and subnetwork figure will be produced in a subdirectory called 'image'. All documents or text files such as the files containing significant gene sets of the hypergeometric test results will be stored in a subdirectory called 'doc'. 
}
\author{
Xin Wang, Camille Terfve
}
\examples{
\dontrun{
##(see the vignette for details about the preprocessing of this data set)
##report for a GSCA object
data("KcViab.gsca")
report(object=KcViab.gsca, experimentName="KcViab", species="Dm", allSig=TRUE, 
keggGSCs="PW.KEGG", goGSCs=c("GO.BP","GO.MF","GO.CC"), reportDir="HTSanalyzerGSCAReport")
##report for a NWA object
data("KcViab.nwa")	
report(object=KcViab.nwa, experimentName="KcViab", species="Dm", allSig=TRUE, 
keggGSCs="PW.KEGG", goGSCs=c("GO.BP","GO.MF","GO.CC"), reportDir="HTSanalyzerNWReport")
##report both analyses
reportAll(gsca=KcViab.gsca, nwa=KcViab.nwa, experimentName="KcViab", species="Dm", allSig=TRUE, 
keggGSCs="PW.KEGG", goGSCs=c("GO.BP","GO.MF","GO.CC"), reportDir="HTSanalyzerReport")	
}
}
