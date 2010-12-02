##This function compute the nominal p-value associated with a GSEA for a
##collection of gene sets, from the outputs of collectionGsea

permutationPvalueCollectionGsea <- function(permScores, dataScores){
	##check 'permScores'
	if(!is.matrix(permScores)) 
		stop("The argument permScores should be a matrix")
	#check 'dataScores'
	if(!is.vector(dataScores) || is.null(names(dataScores))) 
		stop("The argument dataScores should be a named vector")
	if(!is.integer(dataScores) && !is.numeric(dataScores)) 
		stop("The argument dataScores should be a numerical vector")
	##check that the dimensions of permScores and dataScores match to 
	##what is expected
	l.dataScores <- length(dataScores)
	nPerm <- ncol(permScores)
	if(nrow(permScores) != l.dataScores) 
		warning("The number of rows of the permScores matrix is not ",
			"equal to the length of dataScores")
	##initialize a pvalues vector	
	pval <- rep(1, l.dataScores)
	##go through each element of dataScores and see how many permScores 
	##in the corresponding row are higher (in magnitude)
	##this is done separately for negative and positive scores
	##the scores of zero are simply not treated because they are left 
	##at the initial pvalue of 1	
	valid.id <- which(!is.na(dataScores))
	pval <- sapply(valid.id, function(i) {
		ifelse(
			dataScores[i] > 0,
			sum(permScores[i, ] > dataScores[i])/nPerm,
			sum(permScores[i, ] < dataScores[i])/nPerm
		)
		
	}) 
	names(pval)<-names(dataScores)		
	return(pval)
}

