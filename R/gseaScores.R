#Function to compute enrichment scores for GSEA, running score and position of hits for a gene set.
gseaScores <-
function(geneList,geneSet,exponent,mode=c("graph","score")) {
	#Check that the geneList is a named vector
	if(!is.vector(geneList) || is.null(names(geneList))) 
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
	if(!(mode %in% c("graph", "score"))) {
		stop("The mode argument is not correctly specified, please provide one of the following character strings: 'graph', 'score'")
	}
	#The geneSet should be a subset of the gene universe, i.e. we keep only those element of the gene set that appear in the geneList		
	geneSet<-intersect(names(geneList),geneSet)
	#Compute the size of the gene set and of the genelist	
	nh<-length(geneSet)
	N<-length(geneList)
	#Initialize the ES, runningES and the Phit and Pmiss by position (the actual values of Phit and Pmiss are actually cumulative sums of these 'by position' values)	
	ES<-0
	Phit<-rep(0,N)
	Pmiss<-rep(0,N)
	runningES<-rep(0,N)
	#Stop if the geneSet is larger than the gene universe	
	if(nh > N) {
		stop("Gene Set is larger than Gene List")
	} else {
		#Compute the positions of the hits in the geneList (0 if there is no match, 1 if there is a match)	
		hits<-rep(FALSE,N)
		hits[which(!is.na(match(names(geneList),geneSet)))]<-TRUE
		#If sum(hits)=0 then there is no match between geneList and geneSet, and all scores stay at 0.		
		if(sum(hits)!=0) {
			#Fill the Phit by position		
			Phit[which(hits)]<-abs(geneList[which(hits)])^exponent
			NR=sum(Phit)
			#Fill the Pmiss by positions			
			Pmiss[which(!hits)]<-1/(N-nh)
			#Do the cumulative sums	and compute the runningES		
			Phit=cumsum(Phit/NR)
			Pmiss=cumsum(Pmiss)
			runningES<-Phit-Pmiss
			#Compute the maximal (positive) and minimal (or maximal negative) values of the ES, and choose which one is kept			
			ESmax<-max(runningES)
			ESmin<-min(runningES)
			ES<-ifelse(abs(ESmin)>abs(ESmax), ESmin, ESmax)
		}
	}
	#Return the relevant information according to mode		
	if(mode=="score")
		return(ES)	
	if(mode=="graph")	
		return(list("Enrichment.Score"=ES, "Running.Score"=runningES, "Positions"=as.integer(hits)))
}

