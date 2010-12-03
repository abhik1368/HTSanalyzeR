\name{analyze}
\alias{analyze}
\alias{analyze,GSCA-method}
\alias{analyze,NWA-method}
\title{
Gene Set Collection Analysis or NetWork Analysis
}
\description{
This is a generic function.

When implemented as the S4 method for objects of class 'GSCA', this 
function will invoke function \code{analyzeGeneSetCollections} to do 
hypergeometic tests and GSEA. 

When implemented on an object of class \code{NWA}, it will call 
function \code{networkAnalysis} to do subnetwork identification. 

To use this function for objects of class \code{GSCA}:

analyze(object, para = list(pValueCutoff = 0.05, pAdjustMethod = "BH", 
nPermutations = 1000, minGeneSetSize = 15, exponent = 1), verbose = TRUE)

To use this function for objects of class \code{NWA}:

analyze(object, fdr=0.001, species, verbose=TRUE)
}
\usage{
analyze(object, ...)
}
\arguments{
	\item{object}{
an object. When this function is implemented as the S4 method of class 'GSCA' 
or 'NWA', this argument is an object of class 'GSCA' or 'NWA'.
}
	\item{...}{
other arguments depending on class (see below for the arguments supported by 
class 'GSCA' and/or 'NWA')
}
\describe{
	\item{para:}{
a list of parameters for GSEA and hypergeometric tests including:
para$pValueCutoff: a single numeric value specifying the cutoff for p-values 
considered significant	
}
	\item{para$pAdjustMethod:}{
a single character value specifying the p-value adjustment method 
to be used (see 'p.adjust' for details)
}
	\item{para$nPermutations:}{
a single integer or numeric value specifying the number of permutations 
for deriving p-values in GSEA
}
	\item{para$minGeneSetSize:}{
a single integer or numeric value specifying the minimum number of 
elements in a gene set that must map to elements of the gene universe. 
Gene sets with fewer than this number are removed from both hypergeometric 
analysis and GSEA.
}
	\item{para$exponent:}{
a single integer or numeric value used in weighting phenotypes in GSEA 
(see 'gseaScores' function)
}
	\item{fdr:}{
a single numeric value specifying the false discovery for the scoring 
of nodes (see BioNet::scoreNodes and Dittrich et al., 2008 for details)	
}
	\item{species:}{
a single character value specifying the species for which the data 
should be read. The current version support one of the following 
species: "Dm" ("Drosophila_melanogaster"), "Hs" ("Homo_sapiens"), 
"Rn" ("Rattus_norvegicus"), "Mm" ("Mus_musculus"), "Ce" 
("Caenorhabditis_elegans"). This is an optional argument here. If it is 
provided, then the labels of nodes of the identified subnetwork will be
mapped from Entrez IDs to gene symbols; otherwise, Entrez IDs will be
used as labels for those nodes.
}
	\item{verbose:}{
a single logical value deciding to display detailed messages (when 
verbose=TRUE) or not (when verbose=FALSE)
}
}

}
\details{
For objects of class \code{GSCA}:

This function will store the results from function \code{analyzeGeneSetCollections} 
in slot \code{result}, and update information about these results to 
slot \code{summary} of class \code{GSCA}. 

See function \code{analyzeGeneSetCollections} for the detailed information
about the returned results.

For objects of class \code{NWA}:

This function will store the subnetwork module identified by BioNet (if 
\code{species} is given, labels of nodes will also be mapped from Entrez IDs 
to gene symbols), and update information about these results to slot 
\code{summary} of class \code{NWA}.

See function \code{networkAnalysis} for the detailed information about the
returned results.
}
\value{
In the end, this function will return an updated object of class \code{GSCA} 
or \code{NWA}.  
}
\references{
Beisser D, Klau GW, Dandekar T, Muller T, Dittrich MT. \emph{BioNet: an R-Package 
for the functional analysis of biological networks.} Bioinformatics. 2010 
Apr 15;26(8):1129-30.

Dittrich MT, Klau GW, Rosenwald A., Dandekar T and Muller T. 
\emph{Identifying functional modules in protein-protein interaction 
networks: an integrated exact approach.} Bioinformatics 2008 24(13):
i223-i231.
}
\author{
Xin Wang \email{xw264@cam.ac.uk}
}
\seealso{
\code{analyzeGeneSetCollection}, \code{networkAnalysis}
}
\examples{
\dontrun{
library(org.Dm.eg.db)
library(KEGG.db)
##load data for enrichment analyses
data("data4enrich")
##select hits
hits <- names(data4enrich)[which(abs(data4enrich) > 2)]
##set up a list of gene set collections
PW.KEGG <- KeggGeneSets(species = "Dm")
gsc.list <- list(PW.KEGG = PW.KEGG)
##create an object of class 'GSCA'
gsca <- new("GSCA", listOfGeneSetCollections=gsc.list,
geneList = data4enrich, hits = hits)
##print gsca
gsca
##do preprocessing (data4enrich has already been preprocessed)
gsca <- preprocess(gsca, species="Dm", initialIDs = "Entrez.gene", 
keepMultipleMappings = TRUE, duplicateRemoverMethod = "max", 
orderAbsValue = FALSE)
##print gsca again
gsca
##do hypergeometric tests and GSEA
gsca <- analyze(gsca, para = list(pValueCutoff = 0.05, pAdjustMethod 
= "BH", nPermutations = 1000, minGeneSetSize = 100,exponent = 1))
##updated object
gsca
summary(gsca)
}
}










