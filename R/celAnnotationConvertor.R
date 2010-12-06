##This function converts an initial named data vector to the same
##vector but with a different identifier category. This function can
##also take a matrix, with rows=gene id's. This function removes the
##genes for which no mapping were found.

celAnnotationConvertor <- function(geneList, initialIDs="Entrez.gene", 
	finalIDs="Entrez.gene", keepMultipleMappings=TRUE, verbose=TRUE) {
	##check arguments
	paraCheck("genelist.general", geneList)
	paraCheck("cel.initialIDs", initialIDs)
	paraCheck("cel.finalIDs", finalIDs)
	paraCheck("keepMultipleMappings", keepMultipleMappings)
	paraCheck("verbose", verbose)
	fromto<-"dummystring"
	#check the environment to be used for the mapping
	#If the type of initial identifiers is not "Entrez.gene", then the 
	#mapping will automatically be from one of the following to Entrez 
	#Gene identifiers	
	if(initialIDs == "Ensembl.transcript") 
		fromto <- org.Ce.egENSEMBLTRANS2EG
	else if(initialIDs == "Ensembl.prot") 
		fromto <- org.Ce.egENSEMBLPROT2EG
	else if(initialIDs == "Ensembl.gene") 
		fromto <- org.Ce.egENSEMBL2EG
	else if(initialIDs == "RefSeq") 
		fromto <- org.Ce.egREFSEQ2EG
	else if(initialIDs == "Symbol") 
		fromto <- org.Ce.egSYMBOL2EG
	else if(initialIDs == "GenBank") 
		fromto <- org.Ce.egACCNUM2EG
	#If the initial identifiers is 	"Entrez.gene", then the mapping will 
	#automatically be from Entrez Gene identifiers	to one of the following
	if(initialIDs == "Entrez.gene") {
		if(finalIDs == "Ensembl.gene") 
			fromto <- org.Ce.egENSEMBL
		if(finalIDs == "Ensembl.transcript") 
			fromto <- org.Ce.egENSEMBLTRANS			
		if(finalIDs == "Ensembl.prot") 
			fromto <- org.Ce.egENSEMBLPROT
		if(finalIDs == "RefSeq") 
			fromto <- org.Ce.egREFSEQ
		if(finalIDs == "Symbol") 
			fromto <- org.Ce.egSYMBOL
		if(finalIDs == "GenBank") 
			fromto <- org.Ce.egACCNUM
		if(finalIDs == "wormbase") 
			fromto <- org.Ce.egWORMBASE
	}
	#Check that the environment has been correctly determined
	if(class(fromto) != "AnnDbBimap") 
		stop(paste("Please provide a valid type of identifiers for the",
			" 'initialIDs' and 'finalIDs' parameters ",
			"(see help(celAnnotationConvertor))", sep=""))		
	#for a named vector		
	if(!is.matrix(geneList)) {
		#Create a list with an element for each name in the geneList, 
		#containing a vector of identifiers of the type finalIDs mapped 
		#to that name in the geneList		
		list.new.names <- mget(names(geneList), fromto, ifnotfound = NA)
		#Create a vector that will hold the new names, and a vector that 
		#will tag the names that were mapped to multiple identifiers		
		n.new.names <- length(list.new.names)
		new.names <- rep(0, n.new.names)
		tag.multiples <- rep(FALSE, n.new.names)
	
		sapply(1:n.new.names, function(i) {
			new.names[i] <<- list.new.names[[i]][1]			
			if(length(list.new.names[[i]]) > 1) {
				if(keepMultipleMappings){
					if(verbose) {
						cat("--The following identifier was mapped to ",
							"more than one value (only the first value",
							" is kept): \n");
						cat("--", list.new.names[i], "\n")
					}
				} else {
					if(verbose) {
						cat("--The following identifier was mapped to ",
							"more than one value (this entry will be ",
							"discarded): \n");
						cat("--", list.new.names[i], "\n")
					}
					tag.multiples[i] <<- TRUE
				}
			}
			NULL
		})
		#If multiple mappings should be kept
		if(keepMultipleMappings) {	
			newdata <- geneList	
			names(newdata) <- new.names
			newdata <- newdata[!is.na(names(newdata))]
			if(verbose) 
				cat("--", paste((length(geneList)-length(newdata)),
					" genes (out of ", length(geneList),
					") could not be mapped to any identifier, ",
					"and were removed from the data. \n"))
		}
		#If multiple mappings should be discarded
		else {	
			newdata <- geneList[which(!tag.multiples)]	
			names(newdata) <- new.names[which(!tag.multiples)]
			newdata <- newdata[!is.na(names(newdata))]
			if(verbose)
				cat("--", paste((length(geneList) - length(newdata)), 
					" genes (out of ", length(geneList),
					") could not be mapped to any identifier ",
					"(or were mapped to multiple identifiers), ",
					"and were removed from the data. \n"))
		}
	} 
	#if a matrix 
	else {
		list.new.names <- mget(rownames(geneList), fromto, ifnotfound = NA)	
		n.new.names <- length(list.new.names)
		new.names <- rep(0, n.new.names)
		tag.multiples <- rep(0, n.new.names)
		
		sapply(1:n.new.names, function(i) {
			new.names[i]<<-list.new.names[[i]][1]
			if(length(list.new.names[[i]]) > 1) {
				if(keepMultipleMappings) {
					if(verbose) {
						cat("--The following identifier was mapped to ",
							"more than one value (only the first value ",
							"is kept): \n");
						cat("--", list.new.names[i], "\n")
					}
				} else {
					if(verbose) {
						cat("--The following identifier was mapped to ",
							"more than one value (this entry will be ",
							"discarded): \n")
						cat("--", list.new.names[i], "\n")
					}
					tag.multiples[i]<<-TRUE
				}
			}
			NULL
		})		
		if(keepMultipleMappings) {	
			newdata<-geneList	
			rownames(newdata)<-new.names
			newdata<-newdata[!is.na(rownames(newdata)),]
			if(verbose)
				cat("--",paste(((dim(geneList)[1])-(dim(newdata)[1])),
					" genes (out of ",dim(geneList)[1] ,
					") could not be mapped to any identifier, and were ",
					"removed from the data. \n"))
		} else {	
			newdata<-geneList[which(!tag.multiples)]	
			rownames(newdata)<-new.names[which(!tag.multiples)]
			newdata<-newdata[!is.na(rownames(newdata)),]
			if(verbose)
				cat("--",paste(((dim(geneList)[1])-(dim(newdata)[1])),
					" genes (out of ",dim(geneList)[1] ,
					") could not be mapped to any identifier ",
					"(or were mapped to multiple identifiers), 
					and were removed from the data. \n"))
			}
	}
	return(newdata);	
}

