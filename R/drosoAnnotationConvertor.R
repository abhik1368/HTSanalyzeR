##This function converts an initial named data vector to the same vector 
##but with a different identifier category for species Drosophila 
##Melanogaster. This function can also take a matrix, with rows=gene id's. 
##This function removes the genes for which no mapping were found.

drosoAnnotationConvertor <- function(geneList, initialIDs = "Entrez.gene", 
	finalIDs = "Entrez.gene", keepMultipleMappings = TRUE, verbose = TRUE) {
	##check arguments
	paraCheck("genelist.general",geneList)
	paraCheck("dro.initialIDs",initialIDs)
	paraCheck("dro.finalIDs",finalIDs)
	paraCheck("keepMultipleMappings",keepMultipleMappings)
	paraCheck("verbose",verbose)
	##Determine the environment to be used for the mapping
	##If the type of initial identifiers is not "Entrez.gene", then the mapping will 
	##automatically be from one of the following to Entrez Gene identifiers		
#	if(initialIDs == "Ensembl.transcript") 
#		fromto<-org.Dm.egENSEMBLTRANS2EG
#	else if(initialIDs == "Ensembl.prot") 
#		fromto<-org.Dm.egENSEMBLPROT2EG
#	else if(initialIDs == "Ensembl.gene") 
#		fromto<-org.Dm.egENSEMBL2EG
#	else if(initialIDs == "RefSeq") 
#		fromto<-org.Dm.egREFSEQ2EG
#	else if(initialIDs == "Symbol") 
#		fromto<-org.Dm.egSYMBOL2EG
#	else if(initialIDs == "GenBank") 
#		fromto<-org.Dm.egACCNUM2EG
#	else if(initialIDs == "Flybase") 
#		fromto<-org.Dm.egFLYBASE2EG
#	else if(initialIDs == "FlybaseCG") 
#		fromto<-org.Dm.egFLYBASECG2EG
#	else if(initialIDs == "FlybaseProt") 
#		fromto<-org.Dm.egFLYBASEPROT2EG
		
	if(initialIDs == "Ensembl.transcript") 
		fromto<-tryCatch(get("org.Dm.egENSEMBLTRANS2EG"), error=function(e) NULL)
	else if(initialIDs == "Ensembl.prot") 
		fromto<-tryCatch(get("org.Dm.egENSEMBLPROT2EG"), error=function(e) NULL)
	else if(initialIDs == "Ensembl.gene") 
		fromto<-tryCatch(get("org.Dm.egENSEMBL2EG"), error=function(e) NULL)
	else if(initialIDs == "RefSeq") 
		fromto<-tryCatch(get("org.Dm.egREFSEQ2EG"), error=function(e) NULL)
	else if(initialIDs == "Symbol") 
		fromto<-tryCatch(get("org.Dm.egSYMBOL2EG"), error=function(e) NULL)
	else if(initialIDs == "GenBank") 
		fromto<-tryCatch(get("org.Dm.egACCNUM2EG"), error=function(e) NULL)
	else if(initialIDs == "Flybase") 
		fromto<-tryCatch(get("org.Dm.egFLYBASE2EG"), error=function(e) NULL)
	else if(initialIDs == "FlybaseCG") 
		fromto<-tryCatch(get("org.Dm.egFLYBASECG2EG"), error=function(e) NULL)
	else if(initialIDs == "FlybaseProt") 
		fromto<-tryCatch(get("org.Dm.egFLYBASEPROT2EG"), error=function(e) NULL)
	##If the initial identifiers is 	"Entrez.gene", then the mapping will 
	##automatically be from Entrez Gene identifiers	to one of the following	
	if(initialIDs== "Entrez.gene") {
		if(finalIDs == "Ensembl.gene") 
			fromto<-tryCatch(get("org.Dm.egENSEMBL"), error=function(e) NULL)
		else if(finalIDs == "Ensembl.transcript") 
			fromto<-tryCatch(get("org.Dm.egENSEMBLTRANS"), error=function(e) NULL)			
		else if(finalIDs == "Ensembl.prot") 
			fromto<-tryCatch(get("org.Dm.egENSEMBLPROT"), error=function(e) NULL)
		else if(finalIDs == "RefSeq") 
			fromto<-tryCatch(get("org.Dm.egREFSEQ"), error=function(e) NULL)
		else if(finalIDs == "Symbol") 
			fromto<-tryCatch(get("org.Dm.egSYMBOL"), error=function(e) NULL)
		else if(finalIDs == "GenBank") 
			fromto<-tryCatch(get("org.Dm.egACCNUM"), error=function(e) NULL)
		else if(finalIDs == "Flybase") 
			fromto<-tryCatch(get("org.Dm.egFLYBASE"), error=function(e) NULL)
		else if(finalIDs == "FlybaseCG") 
			fromto<-tryCatch(get("org.Dm.egFLYBASECG"), error=function(e) NULL)
		else if(finalIDs == "FlybaseProt") 
			fromto<-tryCatch(get("org.Dm.egFLYBASEPROT"), error=function(e) NULL)
	}	
#	if(initialIDs== "Entrez.gene") {
#		if(finalIDs == "Ensembl.gene") 
#			fromto<-org.Dm.egENSEMBL
#		else if(finalIDs == "Ensembl.transcript") 
#			fromto<-org.Dm.egENSEMBLTRANS			
#		else if(finalIDs == "Ensembl.prot") 
#			fromto<-org.Dm.egENSEMBLPROT
#		else if(finalIDs == "RefSeq") 
#			fromto<-org.Dm.egREFSEQ
#		else if(finalIDs == "Symbol") 
#			fromto<-org.Dm.egSYMBOL
#		else if(finalIDs == "GenBank") 
#			fromto<-org.Dm.egACCNUM
#		else if(finalIDs == "Flybase") 
#			fromto<-org.Dm.egFLYBASE
#		else if(finalIDs == "FlybaseCG") 
#			fromto<-org.Dm.egFLYBASECG
#		else if(finalIDs == "FlybaseProt") 
#			fromto<-org.Dm.egFLYBASEPROT
#	}	
	##Check that the environment has been correctly determined
	annopc<-paste("org", "Dm", "eg", "db", sep=".")
	if(is.null(fromto))
		stop(paste('Please load library ', annopc, 
			' before running this function!', sep=""))
	if(!is(fromto,"AnnDbBimap")) 
		stop(paste("Please provide a valid type of identifiers for the ",
			"'initialIDs' and 'finalIDs' parameters ",
			"(see help(celAnnotationConvertor))",sep=""))		
	##if a named vector				
	if(!is.matrix(geneList)) {
		##Create a list with an element for each name in the geneList, 
		##containing a vector of identifiers of the type finalIDs mapped 
		##to that name in the geneList			
		list.new.names <- mget(names(geneList), fromto, ifnotfound=NA)
		##Create a vector that will hold the new names, and a vector 
		##that will tag the names that were mapped to multiple identifiers		
		n.new.names<-length(list.new.names)
		new.names<-rep(0, n.new.names)
		tag.multiples<-rep(FALSE, n.new.names)
		#Go through the list of names and:
		#1. assign the first result in each element to the corresponding 
		##position in the new names vector
		#2. check if the element of the list contained more than one result
		#3. if the user asked to keep multiple mappings, just inform the 
		##user that this entry was mapped multiple times
		#4. if the user asked to discard multiple mappings, tag this entry 
		##and inform the user that this entry was mapped multiple times		
		sapply(1:n.new.names, function(i) {
			new.names[i]<<-list.new.names[[i]][1]
			if(length(list.new.names[[i]]) > 1) {
				if(keepMultipleMappings){
					if(verbose) {
						cat("--The following identifier was mapped to ",
							"more than one value (only the first value ",
							"is kept): \n")
						print(list.new.names[i])	
					}
				} else {
					if(verbose) {
						cat("--The following identifier was mapped to ",
							"more than one value (this entry will be ",
							"discarded): \n")
						print(list.new.names[i])
					}
					tag.multiples[i] <<- TRUE
				}
			}			
		})		
		##If the user asked to keep multiple mappings, the vector of new 
		##names is used to set the names of the vector as it is and the 
		##entries that could not be mapped to any new names are removed 
		##from the data the user is informed of how many entries were 
		##removed because they were not mapped						
		if(keepMultipleMappings) {	
			newdata <- geneList	
			names(newdata) <- new.names
			newdata <- newdata[!is.na(names(newdata))]
			if(verbose)
				cat("--", paste((length(geneList) - length(newdata)), 
					" genes (out of ", length(geneList) ,
					") could not be mapped to any identifier, ",
					"and were removed from the data. \n"))
		}
		##If the user asked to discard multiple mappings, the data is 
		##trimmed of all entries that were mapped multiple times
		##and so is the vector of new names.  This is done using the 
		##information in the tag.multiples vector					
		else {	
			newdata <- geneList[which(!tag.multiples)]	
			names(newdata) <- new.names[which(!tag.multiples)]
			newdata <- newdata[!is.na(names(newdata))]
			if(verbose)
				cat("--", paste((length(geneList)-length(newdata)), 
					" genes (out of ",length(geneList) ,
					") could not be mapped to any identifier ",
					"(or were mapped to multiple identifiers), 
					and were removed from the data. \n"))
		}
	} else {
		##This is identical to what is done for vectors, except that we 
		##work on row names			
		list.new.names <- mget(rownames(geneList), fromto, ifnotfound = NA)
		n.new.names <- length(list.new.names)
		new.names <- rep(0,n.new.names)
		tag.multiples <- rep(0,n.new.names)
		sapply(1:n.new.names, function(i) {
			new.names[i] <<- list.new.names[[i]][1]
			if(length(list.new.names[[i]]) > 1) {
				if(keepMultipleMappings) {
					if(verbose) {
						cat("--The following identifier was mapped to ",
						"more than one value (only the first value is ",
						"kept): \n")
						print(list.new.names[i])
					}
				} else {
					if(verbose) {
						cat("--The following identifier was mapped to ",
						"more than one value (this entry will be ",
						"discarded): \n")
						print(list.new.names[i])
					}
					tag.multiples[i]<<-TRUE
				}
			}		
		})
		##This is identical to what is done for vectors, except that we 
		##work on row names and that we discard the whole row of data when 
		##there is no mapping						
		if(keepMultipleMappings) {	
			newdata<-geneList	
			rownames(newdata)<-new.names
			newdata<-newdata[!is.na(rownames(newdata)),]
			if(verbose) {
				cat("--",paste(((dim(geneList)[1])-(dim(newdata)[1])),
					" genes (out of ",dim(geneList)[1] ,
					") could not be mapped to any identifier, ",
					"and were removed from the data. \n"))
			}
		} else {	
			newdata<-geneList[which(!tag.multiples)]	
			rownames(newdata)<-new.names[which(!tag.multiples)]
			newdata<-newdata[!is.na(rownames(newdata)),]
			if(verbose)
				cat("--",paste(((dim(geneList)[1])-(dim(newdata)[1])),
					" genes (out of ", dim(geneList)[1] ,
					") could not be mapped to any identifier ",
					"(or were mapped to multiple identifiers), 
					and were removed from the data. \n"))
		}
	}
	return(newdata)	
}

