##This function computes the FDR associated with a permutation-based
##p-value from the GSEA on a list of gene sets

FDRcollectionGsea <- function(permScores, dataScores){
	##check arguments
	if(!is.matrix(permScores)) 
		stop("'permScores' should be a matrix!\n")
	if(!is.numeric(dataScores) && !is.integer(dataScores))
		stop("'dataScores' should be an integer or numeric vector!\n")
	if(is.null(names(dataScores))) 
		stop("'dataScores' should be named (by gene set identifier)")
	if(nrow(permScores) != length(dataScores)) 
		stop(paste("The number of rows of the 'permScores' matrix ",
			"should be the same as the length of the 'dataScores' vector",
			sep=""))
	##create a vector to store the FDRs	
	ldataScores<-length(dataScores)
	FDRgeneset=rep(0,ldataScores)
	##Compute the normalized enrichment score (i.e. divide each ES 
	##(experimental or permutation-based) by the mean of all negative/positive 
	##(depending on the sign of the ES) permutation based scores for that gene set)	
	##This is done gene set by gene set
	sapply(1:ldataScores, function(i) {
		##Get the indices of all negative permutation-based score
		neg<-which(permScores[i,] <= 0)
		##Get the indices of all positive permutation-based score			
		pos<-which(permScores[i,] >= 0)
		##Average the values, separately for positive and negative scores		
		NegAvg<-abs(mean(permScores[i,neg]))
		PosAvg<-abs(mean(permScores[i,pos]))
		##Normalize the permutation-based scores, 
		##separately for negative and positive scores		
		permScores[i,neg]<<-permScores[i,neg]/NegAvg
		permScores[i,pos]<<-permScores[i,pos]/PosAvg
		##Normalize the observed scores, separately for negative and positive scores		
		dataScores[i] <<- ifelse((dataScores[i] < 0), 
			(dataScores[i]/NegAvg), (dataScores[i]/PosAvg))	
	})
	##Compute the total number of negative/positive scores across all 
	##permutations and all gene sets
	negtot <- length(which(permScores <= 0))
	postot <- length(which(permScores >= 0))
	##Compute the FDR by comparing: 
	##for negative observed ES:
		##the number of permutation-based scores under the observed ES 
		##divided by the total number of negative permutation based scores
		##to the number of observed scores under the observed ES divided 
		##by the total number of negative observed scores
	##for positive observed ES:
		##the number of permutation-based scores over the observed ES 
		##divided by the total number of positive permutation based scores
		##to the number of observed scores over the observed ES divided 
		##by the total number of positive observed scores	
	sapply(1:ldataScores, function(i) {
		if(is.na(dataScores[i])) {
			FDRgeneset[i] <<- 1
		} else if(dataScores[i] < 0) {
				FDRgeneset[i] <<- (sum(permScores <= dataScores[i])/negtot)/
					(sum(dataScores <= dataScores[i])/sum(dataScores <= 0))
		} else {
				FDRgeneset[i] <<- (sum(permScores >= dataScores[i])/postot)/
					(sum(dataScores >= dataScores[i])/sum(dataScores >= 0))
		}
		FDRgeneset[i] <<- ifelse(FDRgeneset[i]>1, 1, FDRgeneset[i])
	})
	#name the FDRs (by gene set name) and return the vector		
	names(FDRgeneset) <- names(dataScores)
	return(FDRgeneset)	
}

