multiHyperGeoTest <-
function(GeneSet,GeneList,hits, min.gene.set.size=15,
	pAdjustMethod=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")) {
	
	geneSets<-geneIds(GeneSet)
	results<-lapply(geneSets,hyperGeoTest,GeneList,hits)
	names(results)<-names(GeneSet)
	results<-t(as.data.frame(results))
	
	##remove gene sets with fewer genes than min.gene.set.size in the GeneList
	results<-results[which(results[,"Gene Set Size"]>min.gene.set.size),]
	
	##Adjust pvalues
	adjPvals<-p.adjust(results[,"Pvalue"],method=pAdjustMethod)
	
	results<-cbind(results,adjPvals)
	colnames(results)[7]<-"Adjusted.Pvalue"
	results<-results[order(results[,"Adjusted.Pvalue"]),]
	results
}

