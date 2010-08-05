collectionGsea <-
function(collectionOfGeneSets,geneList,exponent=1,permutationsMatrix,cutoffGeneSetSize){
#Check that the collectionOfGeneSets is in GeneSetCollection format
	if(class(collectionOfGeneSets) != "GeneSetCollection") stop("Please provide a collection of gene sets in 'GeneSetCollection' format (see package GSEABase)")
#Check that the geneList is a named vector
	if(is.null(dim(geneList)) != TRUE | is.null(names(geneList))) stop("Please provide a geneList as a single named and ordered vector")	
#Check that the geneList does not have any NA names 
	if(length(which(is.na(names(geneList)))) != 0) stop("Some of the elements of your geneList do not have a name")
#check that the exponent is a single integer
	if(length(exponent) != 1 | (class(exponent) != "integer" && class(exponent) != "numeric")) stop("The exponent should be a single integer (default (advised) = 1, see help(gseaScores))")
#check that the cutoffGeneSetSize is a single integer
	if(length(cutoffGeneSetSize) != 1 | (class(cutoffGeneSetSize) != "integer" && class(cutoffGeneSetSize) != "numeric")) stop("The 'cutoffGeneSetSize' parameter should be a single integer")		
	print(date())
#check that the permutationsMatrix is a matrix
	if(is.matrix(permutationsMatrix) == FALSE) stop("The permutationsMatrix should be a matrix (see help(permutationsGenerator))")
#check that the dimensions of the permutationsMatrix are coherent (i.e. that we have a permutations of the names in the geneList on each row)
	if(dim(permutationsMatrix)[2] != length(geneList)) warning("The permutationsMatrix does not have the right dimensions, it should have a row for each permutation of the names of your geneList (see help(permutationsGenerator))")
#tag the gene sets that can be used in the analysis, i.e. those that are smaller than the size of the gene list
#and that have more than 'cutoffGeneSetSize' elements that can be found in the geneList	
	nGeneSets<-length(collectionOfGeneSets)
	tagGeneSets<-rep(0,nGeneSets)
	for(i in 1:nGeneSets){
		if(length(geneList) > length(geneIds(collectionOfGeneSets[[i]]))) tagGeneSets[i]<-1
		if(length(intersect(names(geneList),geneIds(collectionOfGeneSets[[i]]))) 
			< cutoffGeneSetSize) tagGeneSets[i]<-0
		}
#check that there are actually some gene sets that pass the max and min cutoffs
	if(sum(tagGeneSets) == 0) stop("There are no gene sets in your collection that pass the cutoffs on size")
#Generate a matrix to store the permutation-based scores, with one row for each gene set (that has been tagged) and one column for each permutation	
	scoresperm<-matrix(rep(0,(dim(permutationsMatrix)[1]*sum(tagGeneSets))),nrow=sum(tagGeneSets))
	rownames(scoresperm)<-names(collectionOfGeneSets)[which(tagGeneSets != 0)]
#Generate a vector to store the experimental scores	one entry for each gene set (that has been tagged)
	scoresObserved<-rep(0,sum(tagGeneSets))
	names(scoresObserved)<-names(collectionOfGeneSets)[which(tagGeneSets != 0)]
#Compute the scores	
	for(i in 1:length(scoresObserved)){
		scoresObserved[i]<-gseaScores(geneList=geneList,
			geneSet=geneIds(collectionOfGeneSets[[which(tagGeneSets == 1)[i]]]),
			exponent=exponent,mode="score")$Enrichment.Score
		scoresperm[i,]<-permutationScoresGsea(permutationsMatrix=permutationsMatrix,geneList=geneList,
			geneSet=geneIds(collectionOfGeneSets[[which(tagGeneSets == 1)[i]]]),
			exponent=exponent)
		}
	print(date())
	print(paste("Number of gene sets used: ",sum(tagGeneSets)))
	return(list("Observed.scores"=scoresObserved,"Permutation.scores"=scoresperm))	
	}

