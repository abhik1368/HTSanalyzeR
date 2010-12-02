##This function finds subnetworks enriched for genes with significant 
##phenotypes.

networkAnalysis <- function(pvalues, graph, fdr=0.001, verbose=TRUE) {
	##check arguments
	paraCheck("pvalues", pvalues)
	paraCheck("interactome", graph)
	paraCheck("fdr", fdr)
	paraCheck("verbose", verbose)
	cat("-Performing network analysis ... \n")
	##store the name of the nodes of the graphNEL object for which we 
	##have p-value information	
	scoredNodes<-intersect(names(pvalues),nodes(graph))
	##check that there are nodes associated with a p-value
	if(length(scoredNodes) == 0) 
		stop("The rownames of your pvalueMatrix do not match to any ",
			"name in the interactionMatrix, check that you have the ",
			"right type of identifiers.")
	if(verbose) 
		cat(paste("--Your network consists of ", length(nodes(graph)), 
			" nodes, of which ", length(scoredNodes), 
			" have an associated p-value", sep=""), "\n")
	##Get the pvalue information for the nodes of the graphNEL object 
	##only, and fit a bum model on these N.B. the fitting of the bum 
	##model will produce a diagnostic plot on the screen, to check the 
	##fitting
	dataForNw <- pvalues[scoredNodes]
	fb <- fitBumModel(dataForNw)
	##Score the nodes of the network	
	##The nodes without pvalues will get a NA value instead of a score
	scores <- scoreNodes(graph, fb = fb, fdr = fdr)
	##Compute the mean score, and set the score of all non-scored nodes 
	##(NAs) to this mean
	meanscore <- mean(scores, na.rm = TRUE)
	scoreswMean <- scores
	scoreswMean[which(is.na(scores))] <- meanscore
	##Find the optimal subnetwork	
	if(verbose) 
		cat("--Computing the optimal subnetwork", "\n")
	module <- runFastHeinz(network = graph, scores = scoreswMean)
	cat("-Network analysis complete \n")
	##Return a graphNEL object consisting of the enriched sub-network
	return(module)
}

