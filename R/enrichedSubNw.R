#This function takes a matrix of p-values generated for a set of genes and searches for enriched sub-networks in an interaction data set from The Biogrid.
enrichedSubNw <-
function(pvaluesMatrix,columns=c("t.test.pvalues.two.samples","t.test.pvalues.one.sample"),species=c("Dm"),fdr=0.001,genetic=FALSE,biogridObject=NA,
		link="http://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-3.0.64/BIOGRID-ORGANISM-3.0.64.tab2.zip",order=TRUE,
		biogridDataDir="Data",verbose=TRUE) {
	#Check that the pvaluesMatrix is a matrix
	if(!is.matrix(pvaluesMatrix)) 
		stop("The pvaluesMatrix should be a matrix")
	#Check that the columns argument matches to column names in the pvaluesMatrix
	#if(length(match(columns,colnames(pvaluesMatrix))) != length(columns)) stop("The names specified in the 'columns' argument do not match column names of your pvaluesMatrix")
	if(!all(columns %in% colnames(pvaluesMatrix))) 
		stop("The names specified in the 'columns' argument do not match column names of your pvaluesMatrix")
	#Check that the pvaluesMatrix has row names
	if(is.null(rownames(pvaluesMatrix))) 
		stop("The pvaluesMatrix should have row names (which should be Entrez Gene identifiers unless you provide a biogridObject with a different type of identifiers)")
	#Create a vector of pvalues or aggregated pvalues if mutliple columns are specified
	if(length(columns) == 1) {
		pvalues<-pvaluesMatrix[,columns]
	} else {
		pvalues<-pvaluesMatrix[,columns]
		pvalues<-aggregatePvals(pvaluesMatrix,order=order)
	}
	names(pvalues)<-rownames(pvaluesMatrix)
	#Download the data from the BioGRID, if no data matrix is specified by the argument 'biogridObject'	
	if(is.null(dim(biogridObject))) {
		#This checks that the species argument corresponds to a single character string that can be matched to the names of
		#one of the BioGRID data files	
		#if(is.character(species) != TRUE | length(species) != 1)
		if(!is.character(species) || length(species) != 1) 
			stop("The species argument does not have the right format")
		#if(is.na(match(species,c("Drosophila_melanogaster","Homo_sapiens","Rattus_norvegicus","Mus_musculus","Caenorhabditis_elegans")))) {
		if(!(species %in% c("Dm","Hs","Rn","Mm","Ce"))) {
			stop("The species argument does not match any of the names recognized by this function, please provide one of the following character strings: 
				'Dm' ('Drosophila_melanogaster'), 'Hs' ('Homo_sapiens'), 'Rn' ('Rattus_norvegicus'), 'Mm' ('Mus_musculus'), 'Ce' ('Caenorhabditis_elegans')")
		}	
		InteractionsData<-biogridDataDownload(link=link,species=species,dataDirectory=biogridDataDir,verbose=verbose)
	} else {
		#If a data matrix is specified, check that it contains the right columns
		if(!is.matrix(biogridObject)) 
			stop("The biogridObject should be a matrix")
		#if(length(match(c("InteractionType","InteractorA","InteractorB"),colnames(biogridObject))) != 3)
		if(!all(c("InteractionType","InteractorA","InteractorB") %in% colnames(biogridObject))) 
			stop("The biogridObject should contain the following named columns: 'InteractionType','InteractorA','InteractorB'")
		InteractionsData=biogridObject
	}
	#If it is specified that genetic interactions should be discarded, remove those rows			
	if(!genetic) 
		InteractionsData<-InteractionsData[-which(InteractionsData[,"InteractionType"]=="genetic"),]
	#Make a graphNEL object from the data 	
	graph<-makeNetwork(source=InteractionsData[,"InteractorA"],target=InteractionsData[,"InteractorB"],edgemode="undirected")
	#Store the name of the nodes of the graphNEL object for which we have pvalue information	
	scoredNodes<-intersect(names(pvalues),nodes(graph))
	#Check that there are nodes associated with a p-value
	if(length(scoredNodes) == 0) 
		stop("The rownames of your pvalueMatrix do not match to any name in the biogridObject, check that you have the right type of identifiers.")
	if(verbose) 
		cat(paste("--Your network consists of ",length(nodes(graph))," nodes, of which ",length(scoredNodes)," have an associated p-value",sep=""),"\n")
	#Get the pvalue information for the nodes of the graphNEL object only, and fit a bum model on these	
	#N.B. the fitting of the bum model will produce a diagnostic plot on the screen, to check the fitting
	dataForNw<-pvalues[scoredNodes]
	fb<-fitBumModel(dataForNw)
	#Score the nodes of the network	
	#The nodes without pvalues will get a NA value instead of a score
	scores<-scoreNodes(graph,fb=fb,fdr=fdr)
	#Compute the mean score, and set the score of all non-scored nodes (NAs) to this mean
	meanscore<-mean(scores,na.rm=TRUE)
	scoreswMean<-scores
	scoreswMean[which(is.na(scores))]<-meanscore
	#Find the optimal subnetwork	
	if(verbose) 
		cat("--Computing the optimal subnetwork","\n")
	module<-runFastHeinz(network=graph,scores=scoreswMean)
	#Return a graphNEL object consisting of the enriched sub-network	
	return(module)
}

