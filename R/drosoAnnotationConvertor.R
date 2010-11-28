#This function takes a named vector or matrix of data and converts the names into another type of identifiers
drosoAnnotationConvertor <-
	function(
			geneList,initialIDs="Entrez.gene", finalIDs="Entrez.gene", keepMultipleMappings=TRUE,verbose=TRUE) {
	#check input arguments
	paraCheck("genelist.general",geneList)
	paraCheck("initialIDs",initialIDs)
	paraCheck("finalIDs",finalIDs)
	paraCheck("keepMultipleMappings",keepMultipleMappings)
	paraCheck("verbose",verbose)
	#Determine the environment to be used for the mapping
	#If the type of initial identifiers is not "Entrez.gene", then the mapping will 
	#automatically be from one of the following to Entrez Gene identifiers		
	if(initialIDs == "Ensembl.transcript") fromto<-org.Dm.egENSEMBLTRANS2EG
	if(initialIDs == "Ensembl.prot") fromto<-org.Dm.egENSEMBLPROT2EG
	if(initialIDs == "Ensembl.gene") fromto<-org.Dm.egENSEMBL2EG
	if(initialIDs == "RefSeq") fromto<-org.Dm.egREFSEQ2EG
	if(initialIDs == "Symbol") fromto<-org.Dm.egSYMBOL2EG
	if(initialIDs == "GenBank") fromto<-org.Dm.egACCNUM2EG
	if(initialIDs == "Flybase") fromto<-org.Dm.egFLYBASE2EG
	if(initialIDs == "FlybaseCG") fromto<-org.Dm.egFLYBASECG2EG
	if(initialIDs == "FlybaseProt") fromto<-org.Dm.egFLYBASEPROT2EG
	#If the initial identifiers is 	"Entrez.gene", then the mapping will 
	#automatically be from Entrez Gene identifiers	to one of the following	
	if(initialIDs== "Entrez.gene") {
		if(finalIDs == "Ensembl.gene") fromto<-org.Dm.egENSEMBL
		if(finalIDs == "Ensembl.transcript") fromto<-org.Dm.egENSEMBLTRANS			
		if(finalIDs == "Ensembl.prot") fromto<-org.Dm.egENSEMBLPROT
		if(finalIDs == "RefSeq") fromto<-org.Dm.egREFSEQ
		if(finalIDs == "Symbol") fromto<-org.Dm.egSYMBOL
		if(finalIDs == "GenBank") fromto<-org.Dm.egACCNUM
		if(finalIDs == "Flybase") fromto<-org.Dm.egFLYBASE
		if(finalIDs == "FlybaseCG") fromto<-org.Dm.egFLYBASECG
		if(finalIDs == "FlybaseProt") fromto<-org.Dm.egFLYBASEPROT
	}	
	#Check that the environment has been correctly determined
	if(!is(fromto,"AnnDbBimap")) 
		stop("Please provide a valid type of identifiers for the 'initialIDs' and 'finalIDs' parameters (see help(celAnnotationConvertor))")		
	#Determine if the geneList is a matrix (named rows) or a named vector 
	#the computation will be done separately depending on this (i.e. based on names or rownames)					
	if(!is.matrix(geneList)) {
		#Create a list with an element for each name in the geneList, containing a vector 
		#of identifiers of the type finalIDs mapped to that name in the geneList			
		list.new.names <- mget(names(geneList), fromto, ifnotfound=NA)
		#Create a vector that will hold the new names, and a vector that will tag the names that were mapped to multiple identifiers		
		n.new.names<-length(list.new.names)
		new.names<-rep(0,n.new.names)
		tag.multiples<-rep(0,n.new.names)
		#Go through the list of names and:
		#1. assign the first result in each element to the corresponding position in the new names vector
		#2. check if the element of the list contained more than one result
		#3. if the user asked to keep multiple mappings, just inform the user that this entry was mapped multiple times
		#4. if the user asked to discard multiple mappings, tag this entry and inform the user that this entry was mapped multiple times				
		for(i in 1:n.new.names) {
			new.names[i]<-list.new.names[[i]][1]
			if(length(list.new.names[[i]]) > 1) {
				if(keepMultipleMappings){
					if(verbose) {
						cat("--The following identifier was mapped to more than one value (only the first value is kept): \n")
						print(list.new.names[i])	
					}
				} else {
					if(verbose) {
						cat("--The following identifier was mapped to more than one value (this entry will be discarded): \n")
						print(list.new.names[i])
					}
					tag.multiples[i]<-1
				}
			}
		}
		#If the user asked to keep multiple mappings, the vector of new names is used to set the names of the vector as it is
		#and the entries that could not be mapped to any new names are removed from the data
		#the user is informed of how many entries were removed because they were not mapped						
		if(keepMultipleMappings) {	
			newdata<-geneList	
			names(newdata)<-new.names
			newdata<-newdata[!is.na(names(newdata))]
			if(verbose)
				cat("--",paste((length(geneList)-length(newdata))," genes (out of ",length(geneList) ,
					") could not be mapped to any identifier, and were removed from the data. \n"))
		}
		#If the user asked to discard multiple mappings, the data is trimmed of all entries that were mapped multiple times
		#and so is the vector of new names.  This is done using the information in the tag.multiples vector					
		else {	
			newdata<-geneList[which(tag.multiples == 0)]	
			names(newdata)<-new.names[which(tag.multiples == 0)]
			newdata<-newdata[!is.na(names(newdata))]
			if(verbose)
				cat("--",paste((length(geneList)-length(newdata))," genes (out of ",length(geneList) ,
					") could not be mapped to any identifier (or were mapped to multiple identifiers), 
					and were removed from the data. \n"))
		}
	} else {
		#This is identical to what is done for vectors, except that we work on row names			
		list.new.names<-mget(rownames(geneList), fromto, ifnotfound=NA)
		n.new.names<-length(list.new.names)
		new.names<-rep(0,n.new.names)
		tag.multiples<-rep(0,n.new.names)
		for(i in 1:n.new.names) {
			new.names[i]<-list.new.names[[i]][1]
			if(length(list.new.names[[i]]) > 1) {
				if(keepMultipleMappings) {
					if(verbose) {
						cat("--The following identifier was mapped to more than one value (only the first value is kept): \n")
						print(list.new.names[i])
					}
				} else {
					if(verbose) {
						cat("--The following identifier was mapped to more than one value (this entry will be discarded): \n")
						print(list.new.names[i])
					}
					tag.multiples[i]<-1
				}
			}
		}
		#This is identical to what is done for vectors, except that we work on row names and that we discard the whole row of data when 
		#there is no mapping						
		if(keepMultipleMappings) {	
			newdata<-geneList	
			rownames(newdata)<-new.names
			newdata<-newdata[!is.na(rownames(newdata)),]
			if(verbose) {
				cat("--",paste(((dim(geneList)[1])-(dim(newdata)[1]))," genes (out of ",dim(geneList)[1] ,
								") could not be mapped to any identifier, and were removed from the data. \n"))
			}
		} else {	
			newdata<-geneList[which(tag.multiples == 0)]	
			rownames(newdata)<-new.names[which(tag.multiples == 0)]
			newdata<-newdata[!is.na(rownames(newdata)),]
			if(verbose)
				cat("--",paste(((dim(geneList)[1])-(dim(newdata)[1]))," genes (out of ",dim(geneList)[1] ,
					") could not be mapped to any identifier (or were mapped to multiple identifiers), 
					and were removed from the data. \n"))
		}
	}
	return(newdata)	
}

