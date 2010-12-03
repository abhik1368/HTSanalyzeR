\name{report}
\alias{report}
\alias{report,GSCA-method}
\alias{report,NWA-method}
\title{
Write html reports for enrichment or network analyses
}
\description{
This is a generic function.

When implemented as the method of class \code{GSCA} or \code{NWA}, 
\code{report} produces reports for either the Gene Set Collection 
Analysis or the Network Analysis. 

To use \code{report} for objects of class \code{GSCA} or \code{NWA}:

report(object, experimentName="Unknown", species=NULL, ntop=NULL,
allSig=FALSE, keggGSCs=NULL, goGSCs=NULL, reportDir="HTSanalyzerReport")
}
\usage{
report(object, ...)
}
\arguments{
  \item{object}{
an object. When implemented as S4 methods of class 'NWA' or 'GSCA', this 
argument is either an object of class \code{GSCA} or \code{NWA}.
}
  \item{...}{
other arguments. (see below for the arguments supported by class 'NWA' 
and/or 'GSCA')
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
the GSEA analysis.  For each gene set collection, plots are produced for
the 'ntop' most significant p-values.  
}
  \item{allSig:}{
a single logical value determining whether or not to generate plots for
all significant gene sets. A gene set is significant if its corresponding
adjusted p-value is less than the \code{pValueCutoff} set in function
\code{analyze}. (see function \code{analyze} for more details)   
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
These functions take the objects of the two wrapper classes (\code{GSCA}
and \code{NWA}) and write a report in user-specified directory. An index
html file containing a summary of all results and hyperlink tables to
more detailed results will be generated in the root directory. The other
html files will be stored in a subdirectory called 'html'. All images
including GSEA plots and subnetwork figure will be produced in a sub-
directory called 'image'. All documents or text files such as the files
containing significant gene sets of the hypergeometric test results will
be stored in a subdirectory called 'doc'. 
}
\author{
Xin Wang, Camille Terfve
}
\examples{
\dontrun{
##(see the vignette for details about the preprocessing of this data set)
##report for a GSCA object
data("KcViab.gsca")
report(object=KcViab.gsca, experimentName="KcViab", species="Dm",
allSig=TRUE, keggGSCs="PW.KEGG", goGSCs=c("GO.BP","GO.MF","GO.CC"),
reportDir="HTSanalyzerGSCAReport")
browseURL(file.path(getwd(), "HTSanalyzerGSCAReport", "index.html"))
##report for a NWA object
data("KcViab.nwa")	
report(object=KcViab.nwa, experimentName="KcViab", species="Dm",
allSig=TRUE, keggGSCs="PW.KEGG", goGSCs=c("GO.BP","GO.MF","GO.CC"),
reportDir="HTSanalyzerNWReport")
browseURL(file.path(getwd(), "HTSanalyzerNWReport", "index.html"))
}
}
