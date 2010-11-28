#This function performs a Mann-Whitney (equivalent to a Wilcoxon rank sum tes): the null hypothesis 
#is that the distributions of x and y differ by a location shift of mu=0 
#and the alternative is that they differ by some other location shift
#It looks for gene sets that are represented towards the 2 different ends of 
#two ranked lists of genes, i.e. whose phenotype distribution is located around two different
#values in the two phenotypes list, rather than spread on the whole list in both lists.
#it takes as input:
	#-gl1: a named vector where names are gene identifiers of the same type as the ones in the 
#gene set collection, and values are the measurement on phenotype 1 corresponding to those genes
	#-gl2: a named vector where names are gene identifiers of the same type as the ones in the 
#gene set collection, and values are the measurement on phenotype 2 corresponding to those genes
	#-gsc: a gene set collection in geneSetCollection format (see the help on function collectionGsea)
	#-cutoff: minimum number of genes in one gene set (see the help on function collectionGsea)
	#-pAdjustMethod: a method for multiple hypothesis testing correction (one of the p.adjust methods)
#it produces as output:	
#a table containing the p value for the Mann-Whitney test, and the adjusted
#p value.  The table is ordered by the p value column.
pairwisePhenoMannWith<-function(gl1, gl2, gsc, minGeneSetSize=15, pAdjustMethod="BH") {
	paraCheck("genelist",gl1)
	paraCheck("genelist",gl2)
	paraCheck("gsc",gsc)
	paraCheck("minGeneSetSize",minGeneSetSize)
	paraCheck("pAdjustMethod",pAdjustMethod)
	
	wilcox.pval<-1
	names(wilcox.pval)<-"dummy"
	for(i in 1:length(gsc)){
		gl1.gs<-match(gsc[[i]],names(gl1))
		gl1.gs<-gl1.gs[!is.na(gl1.gs)]
		gl2.gs<-match(gsc[[i]],names(gl2))
		gl2.gs<-gl2.gs[!is.na(gl2.gs)]
		if(length(gl1.gs >= minGeneSetSize) && length(gl2.gs >= minGeneSetSize)) {
			wilcox.pval<-c(wilcox.pval,wilcox.test(x=gl1[gl1.gs],y=gl2[gl2.gs],alternative ="two.sided")$p.value)
			names(wilcox.pval)[length(wilcox.pval)]<-names(gsc[i])
		}
	}
	wilcox.pval<-wilcox.pval[2:length(wilcox.pval)]
	wilcox.pval<-wilcox.pval[order(wilcox.pval)]
	a.p<-p.adjust(wilcox.pval,method=pAdjustMethod)
	wilcox.pval.table<-cbind(wilcox.pval,a.p)
	rownames(wilcox.pval.table)<-names(wilcox.pval)
	colnames(wilcox.pval.table)<-c("P.value","Adjusted.p.value")
	return(wilcox.pval.table)
}	