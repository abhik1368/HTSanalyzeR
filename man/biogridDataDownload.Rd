\name{biogridDataDownload}
\alias{biogridDataDownload}
\title{
Get the data from The Biogrid for a given species.
}
\description{
This functions downloads the data from the Biogrid, puts it in a new folder called "BiogridData" and extracts the interactions data for a given species. 
}
\usage{
biogridDataDownload(link , species = c("Drosophila_melanogaster", "Homo_sapiens", "Rattus_norvegicus", "Mus_musculus", "Caenorhabditis_elegans"))
}
\arguments{
  \item{link}{
The link (url) from where the data should be downloaded (in tab2 format).  The default link is the one valid on June 7th 2010.
}
  \item{species}{
The species for which the data should be read. One of the following: "Drosophila_melanogaster", "Homo_sapiens", "Rattus_norvegicus", "Mus_musculus", "Caenorhabditis_elegans"
}	
}
\details{
This function is made to work on the tab2 format from the Biogrid (i.e. the first line is a header containing the columns names).
}
\value{
A matrix with one row for each interaction, and three columns: InteractorA, InteractorB (both given by their Entrez Identifiers) and InteractionType (physical or genetic).
}
\references{
Stark et al. \emph{BioGRID: a general repository for interaction datasets.} Nucleic Acids Research 2006 34(Database Issue):D535-D539
}
\author{
Camille Terfve
}
\seealso{
enrichedSubNw, networkAnalysis
}
\examples{
\dontrun{
InteractionsData<-biogridDataDownload(species="Drosophila_melanogaster")
}
}
