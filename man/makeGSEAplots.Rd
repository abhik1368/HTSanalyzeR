\name{makeGSEAplots}
\alias{makeGSEAplots}
\title{
Produce GSEA plots
}
\description{
This is an internal function.
}
\usage{
makeGSEAplots(geneList, geneSet, exponent, filepath, filename, 
output='png', ...)
}
\arguments{
	\item{geneList}{
a numeric or integer vector of phenotypes in descending or ascending 
order with elements named by their Entrez IDs (no duplicates nor NA values)
}
	\item{geneSet}{
a character vector of gene identifiers of a gene set
}
	\item{exponent}{
a single integer or numeric value used in weighting phenotypes in GSEA 
}
	\item{filepath}{
a single character value specifying the directory where the GSEA plot 
will be stored
}
	\item{filename}{
a single character value specifying the name of the file where the 
GSEA plot will be stored
}
	\item{output}{
a single character value specifying the format of output image: "pdf" or
"png"
}
	\item{...}{
other arguments used by the function 'png' or 'pdf' such as width and height
}
}
\details{
This is an internal function.
}
\author{
Benilton Carvalho, Xin Wang
}

\keyword{internal}
