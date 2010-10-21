###############################################################################
# TODO: 
# Xin Wang <xw264@cam.ac.uk>
# Advisor: Florian Markowetz <florian.markowetz@cancer.org.uk> 
# University of Cambridge Deparment of Oncology
# Cancer Research UK - Cambridge Research Institute
# At 17:08:31, on 7 Oct 2010
###############################################################################
#Function to compute enrichment scores for GSEA, running score and position of hits for a gene set.
gseaScoresBatch <-
function(geneList,geneNames.perm,geneSet,exponent,npermutations, mode=c("graph","score")) {

	#Check that the geneList is a named vector
	if(!is.null(dim(geneList)) || is.null(names(geneList))) 
		stop("Please provide a geneList as a single named and ordered vector")	
	#Check that the geneList does not have any NA names 
	if(any(is.na(names(geneList)))) 
		stop("Some of the elements of your geneList do not have a name")
	#check that the exponent is a single integer
	if(!((is.integer(exponent) || is.numeric(exponent)) && length(exponent) == 1))  
		stop("The exponent should be a single integer (default (advised) = 1, see help(gseaScores))")
	#check that the geneSet is a single vector (I am not precising the class of that vector because although characters are expected
	#most of the time, it is actually possible to have numerical values in the case of Entrez IDs for example)
	if(!is.vector(geneSet)) 
		stop("The geneSet should be a vector")
	#check that the mode argument is correctly specified
	if(!(mode %in% c("graph", "score"))) 
		stop("The mode argument is not correctly specified, please provide one of the following character strings: 'graph', 'score'")

	geneList.names<-names(geneList)
	#Compute the size of the gene set and of the genelist	
	nh<-length(geneSet)
	N<-length(geneList)
	#The geneSet should be a subset of the gene universe, i.e. we keep only those element of the gene set that appear in the geneList		
	geneSet<-intersect(geneList.names,geneSet)
	ES<-rep(0,npermutations+1)
	Phit<-matrix(0,nrow=N,ncol=npermutations+1)
	Pmiss<-Phit
	runningES<-NULL

	if(nh>N) {
		stop("Gene Set is larger than Gene List")
	} else {
		hits<-matrix(FALSE,nrow=N,ncol=npermutations+1) 	
		hits[which(!is.na(match(geneNames.perm, geneSet)))]<-TRUE	
		hits<-matrix(hits,ncol=npermutations+1,byrow=FALSE)		
		if(sum(hits[,1])>0) {
			junk<-sapply(1:(npermutations+1),function(i) Phit[which(hits[,i]),i]<<-abs(geneList[which(hits[,i])])^exponent)	
			NR<-colSums(Phit)		
			Pmiss[which(!hits)]<-1/(N-nh)		
			Pmiss<-sapply(1:(npermutations+1),function(i) cumsum(Pmiss[,i]))
			Phit<-sapply(1:(npermutations+1), function(i) cumsum(Phit[,i])/NR[i])		
			runningES<-Phit-Pmiss		
			ESrange<-sapply(1:(npermutations+1),function(i) range(runningES[,i]))
			ES<-sapply(1:(npermutations+1), function(i) ESrange[which.max(abs(ESrange[,i])),i])	
			if(is.list(ES)) ES<-unlist(ES)
		}
	}
	#Return the relevant information according to mode		
	if(mode=="score") {
		ES<-list(scoresObserved=ES[1],scoresperm=ES[2:(npermutations+1)])
		return(ES)	
	} else if(mode=="graph")
		return(list("Enrichment.Score"=ES, "Running.Score"=runningES, "Positions"=lapply(hits, as.integer)))
}


