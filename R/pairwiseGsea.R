#This function performs pairwise GSEA: it looks for gene sets that are specifically 
#overrepresented towards the 2 different ends of two ranked lists of genes,
#it takes as input:
	#-gl1: a named vector where names are gene identifiers of the same type as the ones in the 
#gene set collection, and values are the measurement on phenotype 1 corresponding to those genes
#this vector MUST be ordered (decreasing of increasing)
	#-gl2: a named vector where names are gene identifiers of the same type as the ones in the 
#gene set collection, and values are the measurement on phenotype 2 corresponding to those genes
#this vector MUST be ordered
#phenotypes 1 and 2 must be measured on the same genes, i.e. the two vectors must have the same
#length and their names must match, but the two vectors must be ordered separately, i.e.
#one phenotype vector is ordered based on the values of that phenotype only
	#-gsc: a list of gene sets
	#-exponent: exponent for the GSEA, keep it at one (see the help on function collectionGsea)
	#-nPermutations: number of permutations for the GSEA (see the help on function collectionGsea)
	#-minGeneSetSize: minimum number of genes in one gene set (see the help on function collectionGsea)
	#-pAdjustMethod: a method for multiple hypothesis testing correction (one of the p.adjust methods)
#it produces as output:	
#a table containing the p value for the GSEA, and the observed scores for each of the phenotypes
#independently.  The table is ordered by the p value column.
pairwiseGsea<-function(gl1,gl2,gsc,exponent=1,nPermutations=1000,minGeneSetSize=15,pAdjustMethod="BH") {
	paraCheck("genelist",gl1)
	paraCheck("genelist",gl2)
	paraCheck("gsc",gsc)
	paraCheck("exponent",exponent)
	paraCheck("nPermutations",nPermutations)
	paraCheck("minGeneSetSize",minGeneSetSize)
	paraCheck("pAdjustMethod",pAdjustMethod)
	#Check that the two phenotype vectors contain data on the same genes
	if(!all(names(gl1)==names(gl2))) 
		stop("The two phenotypes must be measured on the same genes")
	#The collectionGsea function will test that the other parameters have the right format
	#Compute the individual enrichment scores and permutation based scores
	collGsea.ph1<-collectionGsea(collectionOfGeneSets=gsc,geneList=gl1,exponent=exponent,nPermutations=nPermutations,minGeneSetSize=minGeneSetSize)
	collGsea.ph2<-collectionGsea(collectionOfGeneSets=gsc,geneList=gl2,exponent=exponent,nPermutations=nPermutations,minGeneSetSize=minGeneSetSize)
	#Compute the differences between scores
	mPh.Obs.Scores<-collGsea.ph1$Observed.scores-collGsea.ph2$Observed.scores
	names(mPh.Obs.Scores)<-names(collGsea.ph1$Observed.scores)
	mPH.Perm.Scores<-collGsea.ph1$Permutation.scores-collGsea.ph2$Permutation.scores
	rownames(mPH.Perm.Scores)<-rownames(collGsea.ph1$Permutation.scores)
	#Compute the p values
	mPh.Obs.Scores2<-abs(mPh.Obs.Scores)
	mPh.pvalue2<-rep(1, length(mPh.Obs.Scores))
	names(mPh.pvalue2)<-names(mPh.Obs.Scores)
	for(i in 1:length(mPh.pvalue2)) {
			mPh.pvalue2[i]<-length(which(abs(mPH.Perm.Scores[i,]) > mPh.Obs.Scores2[i]))/length(mPH.Perm.Scores[i,])
	}
	#Produce the results table
	a.pval2<-p.adjust(mPh.pvalue2,method=pAdjustMethod)
	mPh.table2<-cbind(collGsea.ph1$Observed.scores,collGsea.ph2$Observed.scores,mPh.pvalue2,a.pval2)
	rownames(mPh.table2)<-names(mPh.pvalue2)
	colnames(mPh.table2)<-c("ES.phenotype1","ES.phenotype2","P.value","Adjusted.P.value")
	mPh.table2<-mPh.table2[order(mPh.pvalue2),]
}