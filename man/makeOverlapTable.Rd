\name{makeOverlapTable}
\alias{makeOverlapTable}
\title{
Find and save the overlapped genes between a gene set and hit list
}
\description{
This is an internal function.
}
\usage{
makeOverlapTable(geneSet, hits, mapID, filepath, filename)
}
\arguments{
\item{geneSet}{
a character vector of gene identifiers of a gene set
}
\item{hits}{
a character vector of gene identifiers of the vector of hits
}
\item{mapID}{
a list of targeting identifiers with a name for each entry of the list 
as the identifier to be mapped
}
\item{filepath}{
a single character value specifying the directory where the overlapped 
gene list will be stored
}
\item{filename}{
a single character value specifying the name of the file where the 
overlapped gene list will be stored
}
}
\details{
This is an internal function.
}
\author{
Benilton Carvalho, Xin Wang
}

\keyword{internal}
