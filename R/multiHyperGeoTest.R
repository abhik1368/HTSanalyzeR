multiHyperGeoTest <-
function(GeneSet,GeneList,hits, min.gene.set.size=15,
	pAdjustMethod=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"),verbose=TRUE) {

	l.GeneSet<-length(GeneSet)
	#if verbose, then create a progress bar to monitor computation progress
	if(verbose) 
		pb <- txtProgressBar(style=3)
	results<-list()
	for(i in 1:length(GeneSet)) {
		results[[i]]<-hyperGeoTest(GeneSet[i], GeneList, hits)
		if(verbose) 
			setTxtProgressBar(pb, i/l.GeneSet)
	}
	if(verbose) 
		close(pb)
	if(length(results)>0) {
		results<-t(as.data.frame(results))
		rownames(results)<-names(GeneSet)
		##remove gene sets with fewer genes than min.gene.set.size in the GeneList
		results<-results[which(results[,"Gene Set Size"]>min.gene.set.size),,drop=FALSE]
		##Adjust pvalues
		adjPvals<-p.adjust(results[,"Pvalue"],method=pAdjustMethod)
		results<-cbind(results,adjPvals)
		colnames(results)[ncol(results)]<-"Adjusted.Pvalue"
		results<-results[order(results[,"Adjusted.Pvalue"]),,drop=FALSE]
		return(results)
	} else {
		return(NULL)
	}
}

