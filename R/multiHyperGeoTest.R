##This function performs hypergeometric tests for over-representation 
##of hits, on a list of gene sets. This function applies the 
##hyperGeoTest function to an entire list of gene sets and returns a 
##data frame.

multiHyperGeoTest <- function(collectionOfGeneSets, universe, hits, 
	minGeneSetSize = 15, pAdjustMethod = "BH", verbose = TRUE) {
	##check arguments
	paraCheck("gsc", collectionOfGeneSets)
	paraCheck("universe", universe)
	paraCheck("hits", hits)
	paraCheck("minGeneSetSize", minGeneSetSize)
	paraCheck("pAdjustMethod", pAdjustMethod)
	paraCheck("verbose", verbose)
	l.GeneSet <- length(collectionOfGeneSets)
	geneset.size <- unlist(
		lapply(
			lapply(collectionOfGeneSets, intersect, y = universe), length
		)
	)
	if(all(geneset.size < minGeneSetSize))
		stop(paste("The largest number of overlapped genes of gene ",
			"sets with universe is: ", max(geneset.size), ", which is < ", 
			minGeneSetSize, "!\n", sep = ""))
	geneset.filtered <- which(geneset.size >= minGeneSetSize)
	##if verbose, create a progress bar to monitor computation progress
	if(verbose) 
		pb <- txtProgressBar(style=3)
	results <- t(
			sapply(geneset.filtered, 
				function(i) {
					if(verbose) 
						setTxtProgressBar(pb, i/l.GeneSet)		
					hyperGeoTest(collectionOfGeneSets[i], universe, hits)
				}
			)
		)
	if(verbose) 
		close(pb)
	if(length(results) > 0) {
		##results <- t(as.data.frame(results))
		##rownames(results) <- names(collectionOfGeneSets)
		##remove gene sets with genes < minGeneSetSize in the geneList
		##Adjust pvalues
		adjPvals <- p.adjust(results[, "Pvalue"], method = pAdjustMethod)
		results <- cbind(results, adjPvals)
		colnames(results)[ncol(results)] <- "Adjusted.Pvalue"
		results <- results[order(results[, "Adjusted.Pvalue"]), , drop=FALSE]		
	} else {
		reuslts <- matrix(, nrow=0, ncol=7)
		colnames(results) <- c("Universe Size", "Gene Set Size", 
			"Total Hits", "Expected Hits", "Observed Hits", "Pvalue", 
			"Adjusted.Pvalue")
	}
	return(results)
}

