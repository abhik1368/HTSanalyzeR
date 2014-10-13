##This function performs a Mann-Whitney test for shift in location of 
##genes from gene sets, on a pair of phenotypes. It looks for gene sets 
##that are represented towards the 2 different ends of two ranked lists 
##of genes, i.e. whose phenotype distribution is located around two 
##different values in the two phenotypes list, rather than spread on 
##the whole list in both lists.

pairwisePhenoMannWhit <- function(gl1, gl2, gsc, minGeneSetSize = 15, 
	pAdjustMethod = "BH") {
	##check arguments
	paraCheck("genelist", gl1)
	paraCheck("genelist", gl2)
	paraCheck("gsc", gsc)
	paraCheck("minGeneSetSize", minGeneSetSize)
	paraCheck("pAdjustMethod", pAdjustMethod)

	wilcox.pval <- unlist(lapply(gsc, function(gs) {
		gl1.gs <- match(gs, names(gl1))
		gl1.gs <- gl1.gs[!is.na(gl1.gs)]
		gl2.gs <- match(gs, names(gl2))
		gl2.gs <- gl2.gs[!is.na(gl2.gs)]
		if(length(gl1.gs >= minGeneSetSize) && 
			length(gl2.gs >= minGeneSetSize)) {
			wilcox.test(x = gl1[gl1.gs], 
				y = gl2[gl2.gs], alternative = "two.sided")$p.value
		}				
	}))
	wilcox.pval <- wilcox.pval[order(wilcox.pval)]
	a.p <- p.adjust(wilcox.pval, method = pAdjustMethod)
	wilcox.pval.table <- cbind(wilcox.pval, a.p)
	rownames(wilcox.pval.table) <- names(wilcox.pval)
	colnames(wilcox.pval.table) <- c("P.value", "Adjusted.p.value")
	return(wilcox.pval.table)
}	
