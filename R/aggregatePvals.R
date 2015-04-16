##This function takes as input a matrix of p-values for example obtained 
##from a GSEA on multiple phenotypes, with a row for each gene set and a 
##column for each phenotype and aggregates the p-values by row (i.e. one 
##aggregated p-value for each gene set) according to Fisher or Stouffer's 
##methods.

aggregatePvals <- function(pvalMatrix, method = "fishers", pAdjustMethod 
	= "BH", order=TRUE) {
	##check 'pvalMatrix'
	if(!is(pvalMatrix, "matrix")) 
		stop("The argument pvalMatrix should be a matrix")
	if(is.null(rownames(pvalMatrix))) 
		stop("The argument pvalMatrix should be a matrix with rownames")
	##check 'method'
	if(!(method %in% c("stouffers", "fishers") && length(method)==1 )) 
		stop('The argument method should be one of "stouffers" or 
			"fishers"')
	##check 'pAdjustMethod'
	paraCheck("pAdjustMethod", pAdjustMethod)
	##do aggregation
	aggr.pval <- rep(1, nrow(pvalMatrix))
	names(aggr.pval) <- rownames(pvalMatrix)
	if (method == "fishers") {
		pvalMatrixLogged <- log(pvalMatrix)
		aggr.pval <- (-2) * apply(pvalMatrixLogged, 1, sum)		
		aggr.pval <- pchisq(aggr.pval, df = 2*ncol(pvalMatrix), 
			lower.tail = FALSE)
	}	
	if (method == "stouffers") {
		aggr.pval <- apply(pvalMatrix , 1, 
				function(x) {
					newPval <- pnorm(sum(qnorm(x)) / sqrt(length(x)))
					newPval
				}
		)
	}
	adj.aggr.pval <- p.adjust(aggr.pval, method = pAdjustMethod)	
	aggr.pval.table <- cbind(aggr.pval, adj.aggr.pval)
	rownames(aggr.pval.table)<-names(aggr.pval)
	colnames(aggr.pval.table)<-c("Aggregated.p.value",
		"Adjusted.aggregated.p.value")
	if(order)
		aggr.pval.table <- aggr.pval.table[
			order(aggr.pval.table[, "Aggregated.p.value"]), ]
}	
