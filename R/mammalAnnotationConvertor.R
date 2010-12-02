##This function converts an initial named data vector to the same vector 
##but with a different identifier category for Mammalian species (this 
##function also works on matrices with rows=genes (named)).

mammalAnnotationConvertor <- function(geneList, initialIDs="Entrez.gene", 
	finalIDs="Entrez.gene", species="Hs", keepMultipleMappings=TRUE, 
	verbose=TRUE) {
	##check arguments
	paraCheck("genelist.general", geneList)
	paraCheck("mam.initialIDs", initialIDs)
	paraCheck("mam.finalIDs", finalIDs)
	paraCheck("mam.species", species)
	paraCheck("keepMultipleMappings", keepMultipleMappings)
	paraCheck("verbose", verbose)
	if(species == "Hs") {
		##Determine the environment to be used for the mapping
		##If the type of initial identifiers is not "Entrez.gene", then 
		##the mapping will automatically be from one of the following to 
		##Entrez Gene identifiers		
		if(initialIDs == "Ensembl.transcript") 
			fromto <- org.Hs.egENSEMBLTRANS2EG
		else if(initialIDs == "Ensembl.prot") 
			fromto <- org.Hs.egENSEMBLPROT2EG
		else if(initialIDs == "Ensembl.gene") 
			fromto <- org.Hs.egENSEMBL2EG
		else if(initialIDs == "RefSeq") 
			fromto <- org.Hs.egREFSEQ2EG
		else if(initialIDs == "Symbol") 
			fromto <- org.Hs.egSYMBOL2EG
		else if(initialIDs == "GenBank") fromto<-org.Hs.egACCNUM2EG
		else if(initialIDs== "Entrez.gene") {
			##If the initial identifiers is "Entrez.gene", then the 
			##mapping will automatically be from Entrez Gene identifiers
			##to one of the following		
			if(finalIDs == "Ensembl.gene") 
				fromto<-org.Hs.egENSEMBL
			else if(finalIDs == "Ensembl.transcript") 
				fromto<-org.Hs.egENSEMBLTRANS
			else if(finalIDs == "Ensembl.prot") 
				fromto<-org.Hs.egENSEMBLPROT
			else if(finalIDs == "RefSeq") 
				fromto<-org.Hs.egREFSEQ
			else if(finalIDs == "Symbol") 
				fromto<-org.Hs.egSYMBOL
			else if(finalIDs == "GenBank") 
				fromto<-org.Hs.egACCNUM
		}	
		##Check that the environment has been correctly determined
		if(!is(fromto, "AnnDbBimap"))
			stop("Please provide a valid type of identifiers for the ",
				"'initialIDs' and 'finalIDs' parameters ",
				"(see help(mammalAnnotationConvertor))")	
	} else if(species == "Mm") {
		##Determine the environment to be used for the mapping
		##If the type of initial identifiers is not "Entrez.gene", then 
		##the mapping will automatically be from one of the following 
		##to Entrez Gene identifiers			
		if(initialIDs == "Ensembl.transcript") 
			fromto <- org.Mm.egENSEMBLTRANS2EG
		else if(initialIDs == "Ensembl.prot") 
			fromto <- org.Mm.egENSEMBLPROT2EG
		else if(initialIDs == "Ensembl.gene") 
			fromto <- org.Mm.egENSEMBL2EG
		else if(initialIDs == "RefSeq") 
			fromto <- org.Mm.egREFSEQ2EG
		else if(initialIDs == "Symbol") 
			fromto <- org.Mm.egSYMBOL2EG
		else if(initialIDs == "GenBank") 
			fromto <- org.Mm.egACCNUM2EG
		else if(initialIDs == "Entrez.gene") {
			##If the initial identifiers is "Entrez.gene", then the 
			##mapping will automatically be from Entrez Gene identifiers
			##to one of the following				
			if(finalIDs == "Ensembl.gene") 
				fromto <- org.Mm.egENSEMBL
			else if(finalIDs == "Ensembl.transcript") 
				fromto <- org.Mm.egENSEMBLTRANS
			else if(finalIDs == "Ensembl.prot") 
				fromto <- org.Mm.egENSEMBLPROT
			else if(finalIDs == "RefSeq") 
				fromto <- org.Mm.egREFSEQ
			else if(finalIDs == "Symbol") 
				fromto <- org.Mm.egSYMBOL
			else if(finalIDs == "GenBank") 
				fromto <- org.Mm.egACCNUM
		}	
		##Check that the environment has been correctly determined
		if(!is(fromto, "AnnDbBimap"))
			stop("Please provide a valid type of identifiers for the ",
				"'initialIDs' and 'finalIDs' parameters ",
				"(see help(mammalAnnotationConvertor))")					
	} else if(species == "Rn") {
		##Determine the environment to be used for the mapping
		##If the type of initial identifiers is not "Entrez.gene", then 
		##the mapping will automatically be from one of the following to 
		##Entrez Gene identifiers			
		if(initialIDs == "Ensembl.transcript") 
			fromto <- org.Rn.egENSEMBLTRANS2EG
		else if(initialIDs == "Ensembl.prot") 
			fromto <- org.Rn.egENSEMBLPROT2EG
		else if(initialIDs == "Ensembl.gene") 
			fromto <- org.Rn.egENSEMBL2EG
		else if(initialIDs == "RefSeq") 
			fromto <- org.Rn.egREFSEQ2EG
		else if(initialIDs == "Symbol") 
			fromto <- org.Rn.egSYMBOL2EG
		else if(initialIDs == "GenBank") 
			fromto <- org.Rn.egACCNUM2EG
		else if(initialIDs == "Entrez.gene") {
			##If the initial identifiers is "Entrez.gene", then the 
			##mapping will automatically be from Entrez Gene identifiers	
			##to one of the following				
			if(finalIDs == "Ensembl.gene") 
				fromto <- org.Rn.egENSEMBL
			else if(finalIDs == "Ensembl.transcript") 
				fromto <- org.Rn.egENSEMBLTRANS
			else if(finalIDs == "Ensembl.prot") 
				fromto <- org.Rn.egENSEMBLPROT
			else if(finalIDs == "RefSeq") 
				fromto <- org.Rn.egREFSEQ
			else if(finalIDs == "Symbol") 
				fromto <- org.Rn.egSYMBOL
			else if(finalIDs == "GenBank") 
				fromto <- org.Rn.egACCNUM
		}	
		##Check that the environment has been correctly determined
		if(!is(fromto, "AnnDbBimap")) 
			stop("Please provide a valid type of identifiers for the ",
				"'initialIDs' and 'finalIDs' parameters ",
				"(see help(mammalAnnotationConvertor))")				
	}
	##if a named vector				
	if(!is.matrix(geneList)) {
		##Create a list with an element for each name in the geneList, 
		##containing a vector of identifiers of the type finalIDs mapped 
		##to that name in the geneList			
		list.new.names <- mget(names(geneList), fromto, ifnotfound = NA)
		##Create a vector that will hold the new names, and a vector 
		##that will tag the names that were mapped to multiple identifiers		
		l.new.names<-length(list.new.names)
		new.names<-rep(0,l.new.names)
		tag.multiples<-rep(FALSE,l.new.names)
		##Go through the list of names and:
		##1. assign the first result in each element to the corresponding 
		##position in the new names vector
		##2. check if the element of the list contained more than one result
		##3. if the user asked to keep multiple mappings, just inform the 
		##user that this entry was mapped multiple times
		##4. if the user asked to discard multiple mappings, tag this 
		##entry and inform the user that this entry was mapped multiple times
		sapply(1:l.new.names, function(i){
			new.names[i]<<-list.new.names[[i]][1]
			if(length(list.new.names[[i]]) > 1) {
				if(keepMultipleMappings) {
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
					tag.multiples[i]<<-TRUE
				}
			}		
		})				

		##If the user asked to keep multiple mappings, the vector of new 
		##names is used to set the names of the vector as it is and the 
		##entries that could not be mapped to any new names are removed 
		##from the data the user is informed of how many entries were 
		##removed because they were not mapped						
		if(keepMultipleMappings) {	
			newdata<-geneList	
			names(newdata)<-new.names
			newdata<-newdata[!is.na(names(newdata))]
			if(verbose)
				cat("--",paste((length(geneList)-length(newdata)), 
					" genes (out of ",length(geneList) ,
					") could not be mapped to any identifier, ",
					"and were removed from the data. \n"))
		}
		##If the user asked to discard multiple mappings, the data is 
		##trimmed of all entries that were mapped multiple times
		##and so is the vector of new names.  This is done using the 
		##information in the tag.multiples vector					
		else {	
			newdata<-geneList[which(!tag.multiples)]	
			names(newdata)<-new.names[which(!tag.multiples)]
			newdata<-newdata[!is.na(names(newdata))]
			if(verbose)
				cat("--",paste((length(geneList)-length(newdata)),
					" genes (out of ",length(geneList) ,
					") could not be mapped to any identifier ",
					"(or were mapped to multiple identifiers), ",  
					"and were removed from the data. \n"))
		}
	} else {	
		list.new.names<-mget(rownames(geneList), fromto, ifnotfound=NA)
		##This is identical to what is done for vectors, except that we 
		##work on row names		
		l.new.names<-length(list.new.names)
		new.names<-rep(0,l.new.names)
		tag.multiples<-rep(FALSE,l.new.names)
		sapply(1:l.new.names, function(i) {
			new.names[i]<<-list.new.names[[i]][1]
			if(length(list.new.names[[i]]) > 1) {
				if(keepMultipleMappings) {
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
			if(verbose) 
				cat("--",paste(((dim(geneList)[1])-(dim(newdata)[1])),
					" genes (out of ",dim(geneList)[1] ,
					") could not be mapped to any identifier, ",
					"and were removed from the data. \n"))
		} else {	
			newdata<-geneList[which(!tag.multiples)]	
			rownames(newdata)<-new.names[which(!tag.multiples)]
			newdata<-newdata[!is.na(rownames(newdata)),]
			if(verbose)
				cat("--",paste(((dim(geneList)[1])-(dim(newdata)[1])),
					" genes (out of ",dim(geneList)[1] ,
					") could not be mapped to any identifier ",
					"(or were mapped to multiple identifiers), ",
					"and were removed from the data. \n"))
		}
	}		
	return(newdata)	
}

