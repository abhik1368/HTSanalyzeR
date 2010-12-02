\name{biogridDataDownload}
\alias{biogridDataDownload}
\title{
Get a network interaction matrix from The Biogrid website for a given 
species
}
\description{
This functions downloads the data from the Biogrid, puts it in the user-
specified folder and extracts the interactions data for a given species. 
}
\usage{
biogridDataDownload(link="http://thebiogrid.org/downloads/archives/
Release\%20Archive/BIOGRID-3.0.64/BIOGRID-ORGANISM-3.0.64.tab2.zip", 
species = "Dm", dataDirectory=".", verbose=TRUE))
}
\arguments{
  \item{link}{
the link (url) from where the data should be downloaded (in tab2 format). 
The default link is the one valid on Nov. 25 2010.
}
  \item{species}{
a single character value specifying the species for which the data should 
be read. The current version support one of the following species: "Dm" 
("Drosophila_melanogaster"), "Hs" ("Homo_sapiens"), "Rn" ("Rattus_norvegicus"), 
"Mm" ("Mus_musculus"), "Ce" ("Caenorhabditis_elegans").
}	
  \item{dataDirectory}{
the directory to store downloaded file.
}
  \item{verbose}{
a single logical value to decide to display detailed messages (when 
verbose=TRUE) or not (when verbose=FALSE)
}
}
\details{
This function is made to work on the tab2 format from the Biogrid (i.e. 
the first line is a header containing the columns names).
}
\value{
a matrix with one row for each interaction, and three columns: InteractorA, 
InteractorB (both given by their Entrez Identifiers) and InteractionType 
(physical or genetic).
}
\references{
Stark et al. \emph{BioGRID: a general repository for interaction datasets.} 
Nucleic Acids Research 2006 34(Database Issue):D535-D539
}
\author{
Camille Terfve, Xin Wang
}
\seealso{
networkAnalysis
}
\examples{
\dontrun{
InteractionsData<-biogridDataDownload(species="Dm", dataDirectory="TestDir", 
verbose=TRUE)
}
}
