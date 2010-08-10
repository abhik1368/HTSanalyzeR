#Compute the nominal p-value associated with a GSEA for one gene set
permutationPvalueGsea <-
function(permutationsMatrix,geneList,geneSet,exponent,observedES){
#check that the permutationsMatrix argument is a matrix
	if(is.matrix(permutationsMatrix) == FALSE) stop("The argument 'permutationsMatrix' should be a matrix, see help(permutationsGenerator)")
#check that the observedES is a single numerical value
	if(length(observedES) != 1 | class(observedES) != "numeric") stop("The observedES should be a single numerical value")
#initialize the observed pvalue	
	pval<-1
#initialize the vector of permutation-based scores	
	permutationScores<-rep(0,dim(permutationsMatrix)[1])
#compute all the permutation-based scores: go thorugh each random permutation of the gene names and	
#1.create a fake geneList based on the real data but with names randomly shuffled
#2.compute the ES for that fake geneList
#The gseaScores function will check:
	#that the geneList has the right format
	#that the fake list has been named properly
	#that the exponent is a single integer
	#that the geneSet is a single vector
	for(i in 1:length(permutationScores)){
		permgeneList<-geneList
		names(permgeneList)<-permutationsMatrix[i,]
		permutationScores[i]<-gseaScores(geneList=permgeneList,geneSet=geneSet,exponent=exponent,
			mode="score")$Enrichment.Score
		}
#Compute the pvalue, separatley for negative observed scores and positive observed scores		
	if(observedES < 0){pval=length(which(permutationScores < observedES))/dim(permutationsMatrix)[1]
		}else{
			pval=length(which(permutationScores > observedES))/dim(permutationsMatrix)[1]
			}
	return(pval)		
	}

