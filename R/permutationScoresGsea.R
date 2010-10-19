#!mod: Abandoned
permutationScoresGsea <-
function(permutationsMatrix,geneList,geneSet,exponent){
	#check that the permutationsMatrix argument is a matrix
	if(!is.matrix(permutationsMatrix)) 
		stop("The argument 'permutationsMatrix' should be a matrix, see help(permutationsGenerator)")
	#Initialize the permutations scores vector
	permutationScores<-rep(0,dim(permutationsMatrix)[1]);
	#compute all the permutation-based scores: go through each random permutation of the gene names and	
	#1.create a fake geneList based on the real data but with names randomly shuffled
	#2.compute the ES for that fake geneList
	#The gseaScores function will check:
		#that the geneList has the right format
		#that the fake list has been named properly
		#that the exponent is a single integer
		#that the geneSet is a single vector	
	for(i in 1:length(permutationScores)) {
		permgeneList<-geneList;
		names(permgeneList)<-permutationsMatrix[i,];
		permutationScores[i]<-gseaScores(geneList=permgeneList,geneSet=geneSet,exponent=exponent,mode="score")
	}
	return(permutationScores);	
}

