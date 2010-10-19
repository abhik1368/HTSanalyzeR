#This function computes the overlap between a gene list and a gene set .  
#This might be useful e.g. to check that the identifier mapping was done correctly.
overlapComputing <- function(geneList,geneSet) {
	#Check that the geneList is a named vector
	if(!is.vector(geneList) || is.null(names(geneList))) 
		stop("Please provide a geneList as a single named vector")
	#Check that the geneList does not have length zero
	if(length(geneList) == 0) 
		stop("The geneList has length zero")
	if(!is.vector(geneSet))
		stop("Please provide a geneSet as a single named vector")
	#Compute the overlap between the geneList names and the GeneSet identifiers
	if(!is.character(geneSet)) 
		geneSet<-as.character(geneSet)
	overlap<-length(intersect(names(geneList),geneSet))
	uniqueIds<-length(geneSet)
	#this prints the information from this exploration on the screen			
	cat(paste("Number of genes in your gene list: ", length(geneList), "\n", sep=""))
	cat(paste("Number of genes in your gene set(s): ", uniqueIds, "\n", sep=""))
	cat(paste("Number of genes from your gene list matched to the gene set(s): ", overlap, "(",((overlap/length(geneList))*100),"%)", "\n", sep=""))
}
