collectionGsea <-
function(collectionOfGeneSets,geneList,exponent=1,npermutations,cutoffGeneSetSize,verbose=TRUE) {
	geneList.names <- names(geneList)
	#Check that the collectionOfGeneSets is in GeneSetCollection format
	if(!is.list(collectionOfGeneSets)) 
		stop("Please provide a collection of gene sets as a list")
	#Check that the geneList is a named vector
	if(!is.vector(geneList) || is.null(names(geneList))) 
		stop("Please provide a geneList as a single named and ordered vector")	
	#Check that the geneList does not have any NA names 
	if(any(is.na(geneList.names))) 
		stop("Some of the elements of your geneList do not have a name")
	#check that the exponent is a single integer
	if(!((is.integer(exponent) || is.numeric(exponent)) && length(exponent) == 1))  
		stop("The exponent should be a single integer (default (advised) = 1, see help(gseaScores))")
	#check that the cutoffGeneSetSize is a single integer
	if(!((is.integer(cutoffGeneSetSize) || is.numeric(cutoffGeneSetSize)) && length(cutoffGeneSetSize) == 1))  
		stop("The 'cutoffGeneSetSize' parameter should be a single integer")		
	#check that the permutationsMatrix is a matrix
	if(!((is.integer(npermutations) || is.numeric(npermutations)) && length(npermutations) == 1 && npermutations>0))
		stop("The npermutations should be an integer greater than 0")
	#check that the dimensions of the permutationsMatrix are coherent (i.e. that we have a permutations of the names in the geneList on each row)
	#tag the gene sets that can be used in the analysis, i.e. those that are smaller than the size of the gene list
	#and that have more than 'cutoffGeneSetSize' elements that can be found in the geneList	
	nGeneSets<-length(collectionOfGeneSets)
	tagGeneSets<-rep(FALSE,nGeneSets)
	tagGeneSets[which(unlist(lapply(collectionOfGeneSets, length))<length(geneList))] <- TRUE
	tagGeneSets[which(unlist(lapply(lapply(collectionOfGeneSets, intersect, y=geneList.names),length))<cutoffGeneSetSize)] <- FALSE
	#check that there are actually some gene sets that pass the max and min cutoffs
	n.tagGeneSets <- sum(tagGeneSets)
	if(n.tagGeneSets == 0) 
		stop("There are no gene sets in your collection that pass the cutoffs on size")
	#Generate a matrix to store the permutation-based scores, with one row for each gene set (that has been tagged) and one column for each permutation	
	scoresperm<-matrix(rep(0,(npermutations*n.tagGeneSets)),nrow=n.tagGeneSets)
	rownames(scoresperm)<-names(collectionOfGeneSets)[which(tagGeneSets)]
	#Generate a vector to store the experimental scores	one entry for each gene set (that has been tagged)
	scoresObserved<-rep(0,n.tagGeneSets)
	names(scoresObserved)<-names(collectionOfGeneSets)[which(tagGeneSets)]
	#Compute the scores	
	#create permutation gene list
	perm.gL<-sapply(1:npermutations, function(n) as.integer(names(geneList))[sample(1:length(geneList), length(geneList),replace=FALSE)])
	perm.gL<-cbind(as.integer(names(geneList)),perm.gL)
	#check if package snow has been loaded and a cluster object has been created for HTSanalyzeR	
	if(is(getOption("htsa.cluster"),"cluster") && "package:snow" %in% search()) {
		scores<-gseaScoresBatchParallel(geneList,geneNames.perm=perm.gL,collectionOfGeneSets=collectionOfGeneSets[which(tagGeneSets)],
				exponent=exponent,npermutations=npermutations, mode="score")
		junk<-sapply(1:n.tagGeneSets, function(i) {scoresperm[i,]<<-unlist(scores["scoresperm",i]); scoresObserved[i]<<-unlist(scores["scoresObserved",i])})
	} else {
		if(verbose) 
			pb <- txtProgressBar(style=3)
		for(i in 1:n.tagGeneSets) {
			scores<-gseaScoresBatch(geneList,geneNames.perm=perm.gL,geneSet=as.integer(collectionOfGeneSets[[which(tagGeneSets)[i]]]),
					exponent=exponent,npermutations=npermutations, mode="score")
			scoresObserved[i]<-scores$scoresObserved
			scoresperm[i,]<-scores$scoresperm
			if(verbose) 
				setTxtProgressBar(pb, i/n.tagGeneSets)
		}	
		if(verbose) 
			close(pb)
	}
	return(list("Observed.scores"=scoresObserved,"Permutation.scores"=scoresperm))	
}

