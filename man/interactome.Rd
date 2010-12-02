\name{interactome}
\alias{interactome}
\title{
Create an interactome from BioGrid data sets
}
\description{
This is an S4 method for objects of class 'NWA' to create an interactome
before conducting network analysis.
}
\usage{
interactome(object, interactionMatrix, species, link, reportDir = 
"HTSanalyzerReport", genetic=FALSE, verbose=TRUE)
}
\arguments{
  \item{object}{
an object of class 'NWA'
}
  \item{interactionMatrix}{
an interaction matrix including columns 'InteractionType','InteractorA'
and 'InteractorB'. If this matrix is available, the interactome can be
directly built based on it.
}
  \item{species}{
a single character value specifying the species for which the data should 
be read. The current version support one of the following species: "Dm" 
("Drosophila_melanogaster"), "Hs" ("Homo_sapiens"), "Rn" ("Rattus_norvegicus"), 
"Mm" ("Mus_musculus"), "Ce" ("Caenorhabditis_elegans").
}	
  \item{link}{
the link (url) from where the data should be downloaded (in tab2 format).
The default link is the one valid on Nov. 25 2010.
}
  \item{reportDir}{
a single character value specifying the directory to store reports. Biogrid
data set will be downloaded and stored in a subdirectory called 'Data' in
'reportDir'.
}
  \item{genetic}{
a single logical value. If TRUE, genetic interactions will be kept;
otherwise, they will be removed from the data set.
}
  \item{verbose}{
a single logical value deciding to display detailed messages (when 
verbose=TRUE) or not (when verbose=FALSE)
}
}
\details{
This function provides two options to create an interactome for network
analysis. Users can either input an interaction matrix including columns
'InteractionType', 'InteractionA' and 'InteractionB', or set 'species',
'link' and 'genetic' to download data set from Biogrid and extract
corresponding interactions to build the interactome. 

Another way to set up the interactome is to input a 'graphNEL' object
when the 'NWA' object is created (i.e. nwa=new("NWA", pvalues, phenotypes,
interactome)).    
}
\value{
In the end, this function will return an updated object with slot
'interactome' as an object of class 'graphNEL'.
}

\author{
Xin Wang \email{xw264@cam.ac.uk}
}
\seealso{
\code{biogridDataDownload}, \code{NWA}
}
\examples{
\dontrun{
##example 1, input interactome when initializing an 'NWA' object
##load p-values and phenotypes
data("pvalues","data4enrich")
##load Biogrid interactome for Drosophila Melanogaster
data("Biogrid.Dm.interactome")
##create a NWA (NetWork Analysis) object
nwa <- new("NWA",pvalues=pvalues, phenotypes=data4enrich, 
interactome=Biogrid.Dm.interactome)
##print nwa
nwa

library(BioNet)
##example 2, build an interactome from Biogrid data base
##create a NWA (NetWork Analysis) object
nwa <- new("NWA", pvalues=pvalues, phenotypes=data4enrich)
##print nwa
nwa
##download data from Biogrid and build the interactome
nwa <- interactome(nwa, species="Dm", reportDir="NWATest")
##print nwa again
nwa
}
}










