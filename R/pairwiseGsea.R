##This function performs pairwise GSEA: it looks for gene sets that are 
##specifically over-represented towards the 2 different ends of two 
##ranked lists of genes, in a gene set collection.

pairwiseGsea <- function(gl1, gl2, gsc, exponent = 1, nPermutations = 1000, 
	minGeneSetSize = 15, pAdjustMethod = "BH") {
	##check arguments
	paraCheck("genelist", gl1)
	paraCheck("genelist", gl2)
	paraCheck("gsc", gsc)
	paraCheck("exponent", exponent)
	paraCheck("nPermutations", nPermutations)
	paraCheck("minGeneSetSize", minGeneSetSize)
	paraCheck("pAdjustMethod", pAdjustMethod)
	##Check that the two phenotype vectors contain data on the same genes
	if(!all(names(gl1) == names(gl2))) 
		stop("The two phenotypes must be measured on the same genes")
	##The collectionGsea function will test that the other parameters 
	##have the right format
	##Compute the individual enrichment scores and permutation based scores
	collGsea.ph1 <- collectionGsea(collectionOfGeneSets = gsc, 
		geneList = gl1, exponent = exponent, nPermutations = nPermutations, 
		minGeneSetSize = minGeneSetSize)
	collGsea.ph2 <- collectionGsea(collectionOfGeneSets = gsc, 
		geneList = gl2, exponent = exponent, nPermutations = nPermutations,
		minGeneSetSize = minGeneSetSize)
	##Compute the differences between scores
	mPh.Obs.Scores <- collGsea.ph1$Observed.scores - collGsea.ph2$Observed.scores
	names(mPh.Obs.Scores) <- names(collGsea.ph1$Observed.scores)
	mPH.Perm.Scores <- collGsea.ph1$Permutation.scores - collGsea.ph2$Permutation.scores
	rownames(mPH.Perm.Scores) <- rownames(collGsea.ph1$Permutation.scores)
	##Compute the p values
	mPh.Obs.Scores2 <- abs(mPh.Obs.Scores)
	mPh.pvalue2 <- sapply(1:length(gsc), function(i) {
		sum(abs(mPH.Perm.Scores[i,]) > mPh.Obs.Scores2[i]) / 
			length(mPH.Perm.Scores[i,])	
	})
	names(mPh.pvalue2) <- names(mPh.Obs.Scores)
	##Produce the results table
	a.pval2 <- p.adjust(mPh.pvalue2, method = pAdjustMethod)
	mPh.table2 <- cbind(collGsea.ph1$Observed.scores, 
		collGsea.ph2$Observed.scores, mPh.pvalue2, a.pval2)
	rownames(mPh.table2) <- names(mPh.pvalue2)
	colnames(mPh.table2) <- c("ES.phenotype1", "ES.phenotype2", 
		"P.value", "Adjusted.P.value")
	mPh.table2 <- mPh.table2[order(mPh.pvalue2), ]
}
