#This function generates random permutations of the identifiers in a gene list.  It is used in the GSEA process.
permutationsGeneratorGsea <-
function(nPermutations,geneList) {
	#Check that the geneList is a named vector
	if(!is.vector(geneList) || is.null(names(geneList))) 
		stop("Please provide a geneList as a single named vector")
	#Check that the geneList does not have any NA names 
	if(any(is.na(names(geneList)))) 
		stop("Some of the elements of your geneList do not have a name")
	#check that the nPermutations is a single integer
	if(!((is.integer(nPermutations) || is.numeric(nPermutations)) && length(nPermutations) == 1)) 
		stop("The nPermutations should be a single integer")
	#Generate a matrix with a random permutation of the names of the geneList in each row
	perms<-matrix(rep(names(geneList),nPermutations),nrow=nPermutations,byrow=TRUE)
	for(i in 1:nPermutations) {
		perms[i,]<-sample(names(geneList),length(geneList),replace=FALSE)
	}
	return(perms)
}

