\name{preprocess}
\alias{preprocess}
\alias{preprocess,GSCA-method}
\alias{preprocess,NWA-method}
\title{
A preprocessing method for objects of class GSCA or NWA
}
\description{
This is a generic function.

When implemented as the S4 method for objects of class \code{\link[HTSanalyzeR:GSCA]{GSCA}} or \code{\link[HTSanalyzeR:NWA]{NWA}}, 
this function filters out invalid data, removes duplicated genes, 
converts annotations to Entrez identifiers, etc.

To use this function for objects of class \code{\link[HTSanalyzeR:GSCA]{GSCA}}:

preprocess(object, species="Dm", initialIDs="FlybaseCG", keepMultipleMappings
=TRUE, duplicateRemoverMethod="max", orderAbsValue=FALSE, verbose=TRUE)

To use this function for objects of \code{\link[HTSanalyzeR:NWA]{NWA}}:

preprocess(object, species="Dm", initialIDs="FlybaseCG", keepMultipleMappings
=TRUE, duplicateRemoverMethod="max", verbose=TRUE)
}
\usage{
preprocess(object, ...)
}
\arguments{
	\item{object}{
an object. When this function is implemented as the S4 method of class \code{\link[HTSanalyzeR:GSCA]{GSCA}} 
or \code{\link[HTSanalyzeR:NWA]{NWA}}, this argument is an object of class 'GSCA' or \code{\link[HTSanalyzeR:NWA]{NWA}}.
}
	\item{...}{
other arguments depending on class (see below for the arguments supported by 
class \code{\link[HTSanalyzeR:GSCA]{GSCA}} and/or \code{\link[HTSanalyzeR:NWA]{NWA}})
}
	\describe{
	\item{species:}{
a single character value specifying the species for which the data should
be read. The current version supports one of the following species: "Dm"
("Drosophila_melanogaster"), "Hs" ("Homo_sapiens"), "Rn" ("Rattus_
norvegicus"), "Mm" ("Mus_musculus"), "Ce" ("Caenorhabditis_elegans").
This is an optional argument here. If it is provided, then the labels of
nodes of the identified subnetwork will be mapped from Entrez IDs to gene
symbols; otherwise, Entrez IDs will be used as labels for those nodes.
}
  \item{initialIDs:}{
a single character value specifying the type of initial identifiers for
input 'geneList'. Current version can take one of the following types: 
"Ensembl.transcript", "Ensembl.prot", "Ensembl.gene", "Entrez.gene",
"RefSeq", "Symbol" and "GenBank" for all supported species; "Flybase",
"FlybaseCG" and "FlybaseProt" in addition for Drosophila Melanogaster; 
"wormbase" in addition for Caenorhabditis Elegans. 
}
  \item{keepMultipleMappings:}{
a single logical value. If 'TRUE', the function keeps the entries with
multiple mappings (first mapping is kept). If 'FALSE', the entries with
multiple mappings will be discarded.
}
	\item{duplicateRemoverMethod:}{
a single character value specifying the method to remove the duplicates 
(should the minimum, maximum or average observation for a same construct
be kept). Current version provides "min" (minimum), "max" (maximum),
"average" and "fc.avg" (fold change average). The minimum and maximum
should be understood in terms of absolute values (i.e. min/max effect,
no matter the sign). The fold change average method converts the fold
changes to ratios, averages them and converts the average back to a fold
change.
}
	\item{orderAbsValue:}{
a single logical value indicating whether the values should be converted 
to absolute values and then ordered (if TRUE), or ordered as they are (if
FALSE).	This argument is only for class \code{GSCA}.
}
	\item{verbose:}{
a single logical value suggesting to display detailed messages (when 
verbose=TRUE) or not (when verbose=FALSE)
}
}
}
\details{

This function will do the following preprocessing steps:
\describe{	
	\item{1:}{ filter out p-values (the slot \code{pvalues} of class \code{NWA}), 
	phenotypes (the slot \code{phenotypes} of class \code{NWA}) and data for 
	enrichment (the slot \code{geneList} of class \code{GSCA}) with NA values 
	or without valid names, and invalid gene names (the slot \code{hits} of class 
	\code{GSCA});
}
	\item{2:}{ invoke function \code{duplicateRemover} to remove duplicated genes 
	in the slot \code{pvalues}, \code{phenotypes} of class \code{NWA}, and the
	slot \code{geneList} and \code{hits} of class \code{GSCA};
}
	\item{3:}{ invoke function \code{annotationConvertor} to convert annotations from
	\code{initialIDs} to Entrez identifiers. Please note that the slot \code{hits} and 
	the names of the slot \code{geneList} of class \code{GSCA}, the names of the slot 
	\code{pvalues} and the names of the slot \code{phenotypes} of class \code{NWA} must 
	have the same type of gene annotation specified by \code{initialIDs};
}	
	\item{4:}{ order the data for enrichment decreasingly for objects of 
	class \code{GSCA}. 
}
}
See the function \code{\link[HTSanalyzeR:duplicateRemover]{duplicateRemover}} for more details about how to remove duplicated 
genes.

See the function \code{\link[HTSanalyzeR:annotationConvertor]{annotationConvertor}} for more details about how to convert 
annotations.
}
\value{
In the end, this function will return an updated object of class \code{\link[HTSanalyzeR:GSCA]{GSCA}}
or \code{\link[HTSanalyzeR:NWA]{NWA}}.  
}
\author{
Xin Wang \email{xw264@cam.ac.uk}
}
\seealso{
\code{\link[HTSanalyzeR:duplicateRemover]{duplicateRemover}}, \code{\link[HTSanalyzeR:annotationConvertor]{annotationConvertor}}
}
\examples{
\dontrun{
library(org.Dm.eg.db)
library(KEGG.db)
##load data for enrichment analyses
data("KcViab_Data4Enrich")
##select hits
hits <- names(KcViab_Data4Enrich)[which(abs(KcViab_Data4Enrich) > 2)]
##set up a list of gene set collections
PW_KEGG <- KeggGeneSets(species = "Dm")
gscList <- list(PW_KEGG = PW_KEGG)
##create an object of class 'GSCA'
gsca <- new("GSCA", listOfGeneSetCollections=gscList, geneList =
KcViab_Data4Enrich, hits = hits)
##print gsca
summarize(gsca, what = c("GeneList", "Hits"))
##do preprocessing (KcViab_Data4Enrich has already been preprocessed)
gsca <- preprocess(gsca, species="Dm", initialIDs = "Entrez.gene", 
keepMultipleMappings = TRUE, duplicateRemoverMethod = "max", 
orderAbsValue = FALSE)
##print updated object
summarize(gsca, what = c("GeneList", "Hits"))
}
}










