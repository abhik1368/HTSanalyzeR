analyzeGeneSetCollections <-
function(ListOfGeneSetCollections,
					GeneList,hits,pAdjustMethod="BH",p.value.cutoff=0.05,
					npermutations=1000, min.gene.set.size=15,exponent=1,
					whichSetIsKEGGIds="none",whichSetIsGOIds="none"){
					
		numGeneSetCollections<-length(ListOfGeneSetCollections)
		
		##Check format of inputs
		if (class(ListOfGeneSetCollections)=="GeneSetCollection"){
			stop("GeneSetCollections must be supplied as elements of a list")
		}
		if (typeof(ListOfGeneSetCollections)!="list") stop(ListOfGeneSetCollections,"is not a list")
	
		l<-list()
		for (i in 1:numGeneSetCollections) {
		l[i]<-class(ListOfGeneSetCollections[[i]])
		}	
		if(any(l!="GeneSetCollection")==TRUE) {
			ind<-which (l!="GeneSetCollection")
			stop("ListOfGeneSetCollections[[",ind,"]] is not a GeneSetCollection")
		}	
		
					
			######################
			###Hypergeometric test
			######################
			print(date())
			print("Performing hypergeometric analysis...")

			HGTresults<-lapply(ListOfGeneSetCollections,multiHyperGeoTest,
				GeneList=names(GeneList), hits=hits,
				min.gene.set.size=min.gene.set.size,pAdjustMethod="none")


			
			print(paste("Number of Gene Set Collections: ",numGeneSetCollections))

			pvals<-c()
			##loop combines pvalues from all gene set collections' results
			##into one vector for multiple hypothesis testing pvalue adjustement
			for (i in 1:numGeneSetCollections){
				pv<-HGTresults[[i]][,"Pvalue"]

				pvals<-c(pvals,pv)	
				}

			##Adjustment of pvalues
			HGTpvals<-p.adjust(pvals,method=pAdjustMethod)



			##for loop adds adjusted values to data frame
			##by matching names of the adjusted pvalue vector
			##with the rownames of the results in the given list element, i
			for (i in 1:numGeneSetCollections){
				ind<-match(rownames(HGTresults[[i]]),names(HGTpvals))
				Adjusted.Pvalue<-HGTpvals[ind]
				HGTresults[[i]]<-cbind(HGTresults[[i]],Adjusted.Pvalue)

			}		
			
			##add KEGG names to KEGG hypergeometric result dataframe
			if (length(whichSetIsKEGGIds)>1 || whichSetIsKEGGIds!="none") {
				for (i in 1:length(whichSetIsKEGGIds)){
					new.names<-appendKEGGNames(HGTresults[[whichSetIsKEGGIds[i]]])
					rownames(HGTresults[[whichSetIsKEGGIds[i]]])<-new.names
				}
			}
			##add GO terms to GO hypergeometric result dataframe
			if (length(whichSetIsGOIds)>1 || whichSetIsGOIds!="none" ) {
				for (i in 1:length(whichSetIsGOIds)){
					new.names<-appendGOTerms(HGTresults[[whichSetIsGOIds[i]]])
					rownames(HGTresults[[whichSetIsGOIds[i]]])<-new.names
				}
			}
			

			HGTAll<-NULL
			
			##for loop combines results from each individual gene set collection
			##into one large dataframe of all results
			for (i in 1:numGeneSetCollections) HGTAll<-rbind(HGTAll,HGTresults[[i]])

			##Results are ordered by adjusted p-values from lowest to highest
			HGTAll<-HGTAll[order(HGTAll[,"Adjusted.Pvalue"]),]
			
			##The dataframe of all results is added as the last element
			##in the hypergeometric results list
			HGTresults[[numGeneSetCollections+1]]<-HGTAll

			##
			result.names<-c(names(ListOfGeneSetCollections),"All.collections")
			names(HGTresults)<-result.names
			print(date())
			print("Hypergeometric analysis complete.")
		print("Performing gene set enrichment analysis...")
		
		######################
		###GSEA
		######################
		
		##Generate random permutations GeneList vector for empirical significance
		##estimation
		test.perm<-permutationsGeneratorGsea(nPermutations=npermutations,geneList=GeneList)
		
		##Calculate enrichment scores for all gene sets in all collections
		test.collection<-lapply(ListOfGeneSetCollections,collectionGsea,
		 		geneList=GeneList,exponent=exponent,
				permutationsMatrix=test.perm,cutoffGeneSetSize=min.gene.set.size)
		
		
		##For loop combines observed and permutation scores for all gene set collections
		##so that fdr and pvalue functions can be performed on the entire set of results
		
		obs.scores<-c()
		prm.scores<-NULL
		for (i in 1:numGeneSetCollections){
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
		for (i in 1:numGeneSetCollections){
			a<-match(names(ListOfGeneSetCollections[[i]]),rownames(test.GSEA.results))
			b<-which(a!="NA")
			a<-a[b]
			GSEA.res.mat<-test.GSEA.results[a,]
			GSEA.res.mat<-GSEA.res.mat[order(GSEA.res.mat[,"Adjusted.Pvalue"]),]
			GSEA.results.list[[i]]<-GSEA.res.mat
		}
		
		
		##add KEGG names to KEGG GSEA result dataframe
		if (length(whichSetIsKEGGIds)>1 || whichSetIsKEGGIds!="none") {
			for (i in 1:length(whichSetIsKEGGIds)){
				new.names<-appendKEGGNames(GSEA.results.list[[whichSetIsKEGGIds[i]]])
				rownames(GSEA.results.list[[whichSetIsKEGGIds[i]]])<-new.names
			}
		}
		
		##add GO terms to GO GSEA result dataframe
		if (length(whichSetIsGOIds)>1 || whichSetIsGOIds!="none") {
			for (i in 1:length(whichSetIsGOIds)){
				new.names<-appendGOTerms(GSEA.results.list[[whichSetIsGOIds[i]]])
				rownames(GSEA.results.list[[whichSetIsGOIds[i]]])<-new.names
			}
		}
		
		
		#test.GSEA.results<-test.GSEA.results[order(test.GSEA.results[,"Adjusted.Pvalue"]),]
		
		
		##Produce data frame containing results for all gene set collections
		GSEA.all<-NULL
		for (i in 1:numGeneSetCollections) GSEA.all<-rbind(GSEA.all,GSEA.results.list[[i]])
		
		GSEA.all<-GSEA.all[order(GSEA.all[,"Adjusted.Pvalue"]),]
		
		
		
		#GSEA.results.list[[numGeneSetCollections+1]]<-test.GSEA.results
		
		GSEA.results.list[[numGeneSetCollections+1]]<-GSEA.all
		names(GSEA.results.list)<-result.names
		
		##identify gene set collections with hypergeometric test pvalues < p.value.cutoff
		sign.hgt<-lapply(HGTresults, function(x) {
			pvs<-x[,"Pvalue"]
			a<-which(pvs<p.value.cutoff)
			if (length(a)==0) {
				#print(paste("No gene sets in ",names(x),"have hypergeometric p-values <",p.value.cutoff))
				}else if (length(a)==1){
					x<-x[a,]
				}else{
					x<-x[a,]
					x<-x[order(x[,"Pvalue"]),]
					}
			})
		##identify gene set collections with hypergeometric test adjusted pvalues < p.value.cutoff
		sign.hgt.adj<-lapply(HGTresults, function(x) {
			a<-which(x[,"Adjusted.Pvalue"]<p.value.cutoff)
			if (length(a)==0) {
				#print(paste("No gene sets in ",names(x),"have hypergeometric adjusted p-values <",p.value.cutoff))
				}else if (length(a)==1){
					x<-x[a,]
				}else{
					x<-x[a,]
					
					x<-x[order(x[,"Adjusted.Pvalue"]),]
					}
			})
		
		##identify gene set collections with GSEA pvalues < p.value.cutoff
		sign.gsea<-lapply(GSEA.results.list, function(x) {
			a<-which(x[,"Pvalue"]<p.value.cutoff)
			if (length(a)==0) {
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
			if (length(a)==0) {
				#print(paste("No gene sets in ",names(ListOfGeneSetCollections)[nm],"have GSEA adjusted p-values <",p.value.cutoff))
				}else if (length(a)==1){
					x<-x[a,]
				}else{
					x<-x[a,]
					x<-x[order(x[,"Adjusted.Pvalue"]),]
					}
			})
			
		overlap<-list()
		overlap.adj<-list()
		
		##identify gene set collections with significant pvalues and/or adjusted pvalues from
		##both GSEA and hypergeometric testing
		for (i in 1:(numGeneSetCollections+1)){
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
		
		
		final.results<-list("HyperGeo.results"=HGTresults,"GSEA.results"=GSEA.results.list,
			"Sig.pvals.in.both"=overlap,"Sig.adj.pvals.in.both"=overlap.adj)

		return(final.results)
		
			
}

