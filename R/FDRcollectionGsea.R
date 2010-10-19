#This function computes the FDR associated with a permutation-based p-value from the GSEA on a collection of gene sets
FDRcollectionGsea <-
function(permScores,dataScores){
	#Check that permscores and dataScores have the right format
	if(!is.matrix(permScores)) 
		stop("The argument 'permScores' should be a matrix")
	if(is.matrix(dataScores)) 
		stop("The argument 'dataScores' should be a vector")
	if(nrow(permScores) != length(dataScores)) 
		stop("The number of rows of the permScores matrix should be the same as the length of the dataScores vector")
	if(is.null(names(dataScores))) 
		stop("The dataScores vector should be named (by gene set identifier)")
	#Create a vector to store the FDRs	
	ldataScores<-length(dataScores)
	FDRgeneset=rep(0,ldataScores)
	#Compute the normalized enrichment score (i.e. divide each ES (experimental or permutation-based) by the mean of all negative/positive 
	#(depending on the sign of the ES) permutation based scores for that gene set)	
	#This is done gene set by gene set
	for(i in 1:ldataScores) {
		#Get the indices of all negative permutation-based score for that gene set	
		neg<-which(permScores[i,] <= 0)
		#Get the indices of all positive permutation-based score for that gene set			
		pos<-which(permScores[i,] >= 0)
		#Average the values, separately for positive and negative scores		
		NegAvg<-abs(mean(permScores[i,neg]))
		PosAvg<-abs(mean(permScores[i,pos]))
		#If the score is 0, just do not normalize anything	
		if(is.na(dataScores[i]) || dataScores[i]==0) {
			NegAvg<-1
			PosAvg<-1
		}
		#Normalize the permutation-based scores, separately for negative and positive scores		
		permScores[i,neg]=permScores[i,neg]/NegAvg
		permScores[i,pos]=permScores[i,pos]/PosAvg
		#Normalize the observed scores, separately for negative and positive scores		
		dataScores[i]<-ifelse((dataScores[i] < 0), (dataScores[i]/NegAvg), (dataScores[i]/PosAvg))
	}
	#Compute the total number of negative/positive scores across all permutations and all gene sets
	negtot<-length(which(permScores <= 0))
	postot<-length(which(permScores >= 0))
	#Compute the FDR by comparing: 
	#for negative observed ES:
		#the number of permutation-based scores under the observed ES divided by the total number of negative permutation based scores
		#to the number of observed scores under the observed ES divided by the total number of negative observed scores
	#for positive observed ES:
		#the number of permutation-based scores over the observed ES divided by the total number of positive permutation based scores
		#to the number of observed scores over the observed ES divided by the total number of positive observed scores	
	for(i in 1:ldataScores) {
		if(is.na(dataScores[i])) {
			FDRgeneset[i]=1
		} else {
			if(dataScores[i]==0) {
				FDRgeneset[i]=1
			} else {
				if(dataScores[i] < 0) {
					FDRgeneset[i]=(length(which(permScores <= dataScores[i]))/negtot)/
						(length(which(dataScores <= dataScores[i]))/length(which(dataScores <= 0)))
				} else {
					FDRgeneset[i]=(length(which(permScores >= dataScores[i]))/postot)/
						(length(which(dataScores >= dataScores[i]))/length(which(dataScores >= 0)))
				}
			}
		}			
	}	
	#name the FDRs (by gene set name) and return the vector		
	names(FDRgeneset)<-names(dataScores)
	return(FDRgeneset)	
}

