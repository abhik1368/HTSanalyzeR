analyzeGeneSetCollections <-
function(ListOfGeneSetCollections,
			GeneList,hits,pAdjustMethod="BH",p.value.cutoff=0.05,
			npermutations=1000, min.gene.set.size=15,exponent=1,
			whichSetIsKEGGIds="none",whichSetIsGOIds="none",verbose=TRUE) {
	##Check format of inputs
	if(!is.list(ListOfGeneSetCollections))
		stop("GeneSetCollections must be supplied as elements of a list")
	######################
	###Hypergeometric test
	######################
	numGeneSetCollections<-length(ListOfGeneSetCollections)
	cat("-Performing hypergeometric analysis ...\n")
	HGTresults<-list()
	for(i in 1:length(ListOfGeneSetCollections)) {
		if(verbose) {
			cat("--For",names(ListOfGeneSetCollections)[i],"\n")
		}
		HGTresults[[i]]<-multiHyperGeoTest(ListOfGeneSetCollections[[i]], GeneList=names(GeneList), hits=hits,
			min.gene.set.size=min.gene.set.size, pAdjustMethod="none", verbose=verbose)
	}
	pvals<-NULL
	##loop combines pvalues from all gene set collections' results
	##into one vector for multiple hypothesis testing pvalue adjustement
	for(i in 1:numGeneSetCollections) {
		pv<-HGTresults[[i]][,"Pvalue"]
		pvals<-c(pvals,pv)	
	}
	##Adjustment of pvalues
	HGTpvals<-p.adjust(pvals,method=pAdjustMethod)
	##for loop adds adjusted values to data frame
	##by matching names of the adjusted pvalue vector
	##with the rownames of the results in the given list element
	for(i in 1:numGeneSetCollections) {
		ind<-match(rownames(HGTresults[[i]]),names(HGTpvals))
		Adjusted.Pvalue<-HGTpvals[ind]
		HGTresults[[i]]<-cbind(HGTresults[[i]],Adjusted.Pvalue)
	}		
	##add KEGG names to KEGG hypergeometric result dataframe
	#if(length(whichSetIsKEGGIds)>1 || whichSetIsKEGGIds!="none") {
	#	for(i in 1:length(whichSetIsKEGGIds)) {
	#		new.names<-appendKEGGNames(HGTresults[[whichSetIsKEGGIds[i]]])
	#		rownames(HGTresults[[whichSetIsKEGGIds[i]]])<-new.names
	#	}
	#}
	##add GO terms to GO hypergeometric result dataframe
	#if(length(whichSetIsGOIds)>1 || whichSetIsGOIds!="none" ) {
	#	for(i in 1:length(whichSetIsGOIds)) {
	#		new.names<-appendGOTerms(HGTresults[[whichSetIsGOIds[i]]])
	#		rownames(HGTresults[[whichSetIsGOIds[i]]])<-new.names
	#	}
	#}
	HGTAll<-NULL
	##for loop combines results from each individual gene set collection
	##into one large dataframe of all results
	for (i in 1:numGeneSetCollections) HGTAll<-rbind(HGTAll,HGTresults[[i]])
	##Results are ordered by adjusted p-values from lowest to highest
	HGTAll<-HGTAll[order(HGTAll[,"Adjusted.Pvalue"]),]
	##The dataframe of all results is added as the last element
	##in the hypergeometric results list
	HGTresults[[numGeneSetCollections+1]]<-HGTAll

	result.names<-c(names(ListOfGeneSetCollections),"All.collections")
	names(HGTresults)<-result.names
	cat("-Hypergeometric analysis complete\n")		
	######################
	###GSEA
	######################
	cat("\n")
	cat("-Performing gene set enrichment analysis ...\n")
	##Calculate enrichment scores for all gene sets in all collections
	test.collection<-list()
	for(i in 1:length(ListOfGeneSetCollections)) {
		if(verbose) {
			cat("--For",names(ListOfGeneSetCollections)[i],"\n")
		}
		test.collection[[i]]<-collectionGsea(ListOfGeneSetCollections[[i]], geneList=GeneList,exponent=exponent,
				npermutations=npermutations, cutoffGeneSetSize=min.gene.set.size,verbose=verbose)
	}
	##For loop combines observed and permutation scores for all gene set collections
	##so that fdr and pvalue functions can be performed on the entire set of results
	obs.scores<-NULL
	prm.scores<-NULL
	for(i in 1:numGeneSetCollections) {
		obs.scores<-c(obs.scores, test.collection[[i]]$Observed.scores)
		prm.scores<-rbind(prm.scores,test.collection[[i]]$Permutation.scores)
	}
	total.test.collect<-list("Observed.scores"=obs.scores,"Permutation.scores"=prm.scores)
	test.FDR.collection<-FDRcollectionGsea(permScores=total.test.collect$Permutation.scores,
			dataScores=total.test.collect$Observed.scores)
	test.pvalues.collection<-permutationPvalueCollectionGsea(
			permScores=total.test.collect$Permutation.scores,
			dataScores=total.test.collect$Observed.scores)		
	gsea.adjust.pval<-p.adjust(test.pvalues.collection,method=pAdjustMethod)
	test.GSEA.results<-cbind(total.test.collect$Observed.scores,test.pvalues.collection,
	gsea.adjust.pval,test.FDR.collection)			
	colnames(test.GSEA.results)<-c("Observed.score","Pvalue","Adjusted.Pvalue","FDR")
	GSEA.results.list<-vector("list",(numGeneSetCollections+1))
	##for loop extracts results dataframe for each gene set collection
	##and orders them by adjusted pvalue
	for(i in 1:numGeneSetCollections) {
		match.ind<-match(names(ListOfGeneSetCollections[[i]]),rownames(test.GSEA.results))
		match.ind<-match.ind[which(!is.na(match.ind))]
		GSEA.res.mat<-test.GSEA.results[match.ind,,drop=FALSE]
		GSEA.res.mat<-GSEA.res.mat[order(GSEA.res.mat[,"Adjusted.Pvalue"]),,drop=FALSE]
		GSEA.results.list[[i]]<-GSEA.res.mat
	}	
	##add KEGG names to KEGG GSEA result dataframe
	#if(length(whichSetIsKEGGIds)>1 || whichSetIsKEGGIds!="none") {
	#	for(i in 1:length(whichSetIsKEGGIds)) {
	#		new.names<-appendKEGGNames(GSEA.results.list[[whichSetIsKEGGIds[i]]])
	#		rownames(GSEA.results.list[[whichSetIsKEGGIds[i]]])<-new.names
	#	}
	#}
	
	##add GO terms to GO GSEA result dataframe
	#if(length(whichSetIsGOIds)>1 || whichSetIsGOIds!="none") {
	#	for(i in 1:length(whichSetIsGOIds)) {
	#		new.names<-appendGOTerms(GSEA.results.list[[whichSetIsGOIds[i]]])
	#		rownames(GSEA.results.list[[whichSetIsGOIds[i]]])<-new.names
	#	}
	#}
	##Produce data frame containing results for all gene set collections
	GSEA.all<-NULL
	for(i in 1:numGeneSetCollections) GSEA.all<-rbind(GSEA.all,GSEA.results.list[[i]])
	
	GSEA.all<-GSEA.all[order(GSEA.all[,"Adjusted.Pvalue"]),]
	
	GSEA.results.list[[numGeneSetCollections+1]]<-GSEA.all
	names(GSEA.results.list)<-result.names
	
	##identify gene set collections with hypergeometric test pvalues < p.value.cutoff
	sign.hgt<-lapply(HGTresults, function(x) {
		pvs<-x[,"Pvalue"]
		a<-which(pvs<p.value.cutoff)
		if (length(a)==0) {
			#print(paste("No gene sets in ",names(x),"have hypergeometric p-values <",p.value.cutoff))
			} else if(length(a)==1) {
				x<-x[a,]
			} else {
				x<-x[a,]
				x<-x[order(x[,"Pvalue"]),]
			}
		})
	##identify gene set collections with hypergeometric test adjusted pvalues < p.value.cutoff
	sign.hgt.adj<-lapply(HGTresults, function(x) {
		a<-which(x[,"Adjusted.Pvalue"]<p.value.cutoff)
		if (length(a)==0) {
			#print(paste("No gene sets in ",names(x),"have hypergeometric adjusted p-values <",p.value.cutoff))
			} else if (length(a)==1){
				x<-x[a,]
			} else {
				x<-x[a,]
				
				x<-x[order(x[,"Adjusted.Pvalue"]),]
			}
		})
	##identify gene set collections with GSEA pvalues < p.value.cutoff
	sign.gsea<-lapply(GSEA.results.list, function(x) {
		a<-which(x[,"Pvalue"]<p.value.cutoff)
		if(length(a)==0) {
			#print(paste("No gene sets in ",names(x),"have GSEA p-values <",p.value.cutoff))
			}else if (length(a)==1){
				x<-x[a,]
			}else{
				x<-x[a,]
				x<-x[order(x[,"Pvalue"]),]
				}
		})
	##identify gene set collections with adjusted GSEA pvalues < p.value.cutoff	
	sign.gsea.adj<-lapply(GSEA.results.list, function(x) {
		a<-which(x[,"Adjusted.Pvalue"]<=p.value.cutoff)
		if(length(a)==0) {
			#print(paste("No gene sets in ",names(ListOfGeneSetCollections)[nm],"have GSEA adjusted p-values <",p.value.cutoff))
			} else if(length(a)==1) {
				x<-x[a,]
			} else {
				x<-x[a,]
				x<-x[order(x[,"Adjusted.Pvalue"]),]
			}
		})
	overlap<-list()
	overlap.adj<-list()
	##identify gene set collections with significant pvalues and/or adjusted pvalues from
	##both GSEA and hypergeometric testing
	for(i in 1:(numGeneSetCollections+1)) {
		a1<-intersect(rownames(sign.gsea[[i]]),rownames(sign.hgt[[i]]))
		a2<-intersect(rownames(sign.gsea.adj[[i]]),rownames(sign.hgt.adj[[i]]))
	
		Hypergeometric.Pvalue<-HGTresults[[i]][a1,"Pvalue"]
		Hypergeometric.Adj.Pvalue<-HGTresults[[i]][a2,"Adjusted.Pvalue"]
	
		GSEA.Pvalue<-GSEA.results.list[[i]][a1,"Pvalue"]
		GSEA.Adj.Pvalue<-GSEA.results.list[[i]][a2,"Adjusted.Pvalue"]
	
		overlap[[i]]<-cbind(Hypergeometric.Pvalue,GSEA.Pvalue)
		overlap.adj[[i]]<-cbind(Hypergeometric.Adj.Pvalue,GSEA.Adj.Pvalue)
	}
	names(overlap)<-result.names	
	names(overlap.adj)<-result.names
	cat("-Gene set enrichment analysis complete \n")
	final.results<-list("HyperGeo.results"=HGTresults,"GSEA.results"=GSEA.results.list,
		"Sig.pvals.in.both"=overlap,"Sig.adj.pvals.in.both"=overlap.adj)
	return(final.results)
}

