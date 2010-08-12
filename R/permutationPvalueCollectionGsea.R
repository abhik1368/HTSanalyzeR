#Compute the nominal p-value associated with a GSEA for a collection of gene sets, from the outputs of collectionGsea
permutationPvalueCollectionGsea <-
function(permScores,dataScores){
	print(date())
#Check that the permScores is a numerical matrix
	if(is.matrix(permScores) == FALSE) stop("The argument permScores should be a matrix")
#Check that dataScores is a numerical vector, named
	if(is.null(dim(dataScores)) != TRUE | is.null(names(dataScores))) stop("The argument dataScores should be a named vector")
	if(class(dataScores) != "integer" && class(dataScores)!= "numeric") stop("The argument dataScores should be a numerical vector")
#Check that the dimensions of permScores and dataScores match to what is expected
	l.dataScores<-length(dataScores)
	if(dim(permScores)[1] != l.dataScores) warning("The number of rows of the permScores matrix is not equal to the length of dataScores")
#Initialize a pvalues vector	
	pval<-rep(1,l.dataScores)
#Go through each element of dataScores and see how many permScores in the corresponding row are higher (in magnitude)
#this is done separately for negative and positive scores
#the scores of zero are simply not treated because they are left at the initial pvalue of 1	
	for(i in 1:l.dataScores){
		if(is.na(dataScores[i])){
			pval[i]=1
			}else{
				if(dataScores[i] < 0){
					pval[i]=length(which(permScores[i,] < dataScores[i]))/dim(permScores)[2]
					}else{
						pval[i]=length(which(permScores[i,] > dataScores[i]))/dim(permScores)[2]
						}
					}	
		}			
	names(pval)<-names(dataScores)		
	print(date())
	return(pval)
	}

