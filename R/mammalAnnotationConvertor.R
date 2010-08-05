#This function takes a named vector or matrix of data and converts the names into another type of identifiers
mammalAnnotationConvertor <-
function(geneList,initialIDs=c(
	"Ensembl.transcript","Ensembl.prot","Ensembl.gene","Entrez.gene","RefSeq","Symbol","GenBank"),
	finalIDs=c(
	"Ensembl.transcript","Ensembl.prot","Ensembl.gene","Entrez.gene","RefSeq","Symbol","GenBank"),
	species=c("Homo_sapiens","Rattus_norvegicus","Mus_musculus"),keepMultipleMappings=TRUE){
#This checks that the species argument corresponds to a single character string that can be matched to one of our species
if(is.character(species) != TRUE | length(species) != 1) stop("The species argument does not have the right format")
if(is.na(match(species,c("Homo_sapiens","Rattus_norvegicus","Mus_musculus")))) {
	stop("The species argument does not match any of the names recognized by this function, please provide one of the following character strings: 'Homo_sapiens','Rattus_norvegicus','Mus_musculus'")
	}	
	if(species == "Homo_sapiens"){
#Determine the environment to be used for the mapping
#If the type of initial identifiers is not "Entrez.gene", then the mapping will 
#automatically be from one of the following to Entrez Gene identifiers		
		if(initialIDs == "Ensembl.transcript") fromto<-org.Hs.egENSEMBLTRANS2EG
		if(initialIDs == "Ensembl.prot") fromto<-org.Hs.egENSEMBLPROT2EG
		if(initialIDs == "Ensembl.gene") fromto<-org.Hs.egENSEMBL2EG
		if(initialIDs == "RefSeq") fromto<-org.Hs.egREFSEQ2EG
		if(initialIDs == "Symbol") fromto<-org.Hs.egSYMBOL2EG
		if(initialIDs == "GenBank") fromto<-org.Hs.egACCNUM2EG
		if(initialIDs== "Entrez.gene") {
#If the initial identifiers is 	"Entrez.gene", then the mapping will 
#automatically be from Entrez Gene identifiers	to one of the following		
			if(finalIDs == "Ensembl.gene") fromto<-org.Hs.egENSEMBL
			if(finalIDs == "Ensembl.transcript") fromto<-org.Hs.egENSEMBLTRANS
			if(finalIDs == "Ensembl.prot") fromto<-org.Hs.egENSEMBLPROT
			if(finalIDs == "RefSeq") fromto<-org.Hs.egREFSEQ
			if(finalIDs == "Symbol") fromto<-org.Hs.egSYMBOL
			if(finalIDs == "GenBank") fromto<-org.Hs.egACCNUM
			}	
#Check that the environment has been correctly determined
		if(class(fromto) != "AnnDbBimap") {
			stop("Please provide a valid type of identifiers for the 'initialIDs' and 'finalIDs' parameters (see help(mammalAnnotationConvertor))")	
			}	
		}	
	if(species == "Mus_musculus"){
#Determine the environment to be used for the mapping
#If the type of initial identifiers is not "Entrez.gene", then the mapping will 
#automatically be from one of the following to Entrez Gene identifiers			
		if(initialIDs == "Ensembl.transcript") fromto<-org.Mm.egENSEMBLTRANS2EG
		if(initialIDs == "Ensembl.prot") fromto<-org.Mm.egENSEMBLPROT2EG
		if(initialIDs == "Ensembl.gene") fromto<-org.Mm.egENSEMBL2EG
		if(initialIDs == "RefSeq") fromto<-org.Mm.egREFSEQ2EG
		if(initialIDs == "Symbol") fromto<-org.Mm.egSYMBOL2EG
		if(initialIDs == "GenBank") fromto<-org.Mm.egACCNUM2EG
		if(initialIDs== "Entrez.gene") {
#If the initial identifiers is 	"Entrez.gene", then the mapping will 
#automatically be from Entrez Gene identifiers	to one of the following				
			if(finalIDs == "Ensembl.gene") fromto<-org.Mm.egENSEMBL
			if(finalIDs == "Ensembl.transcript") fromto<-org.Mm.egENSEMBLTRANS
			if(finalIDs == "Ensembl.prot") fromto<-org.Mm.egENSEMBLPROT
			if(finalIDs == "RefSeq") fromto<-org.Mm.egREFSEQ
			if(finalIDs == "Symbol") fromto<-org.Mm.egSYMBOL
			if(finalIDs == "GenBank") fromto<-org.Mm.egACCNUM
			}	
#Check that the environment has been correctly determined
		if(class(fromto) != "AnnDbBimap") {
			stop("Please provide a valid type of identifiers for the 'initialIDs' and 'finalIDs' parameters (see help(mammalAnnotationConvertor))")	
			}				
		}	
	if(species == "Rattus_norvegicus"){
#Determine the environment to be used for the mapping
#If the type of initial identifiers is not "Entrez.gene", then the mapping will 
#automatically be from one of the following to Entrez Gene identifiers			
		if(initialIDs == "Ensembl.transcript") fromto<-org.Rn.egENSEMBLTRANS2EG
		if(initialIDs == "Ensembl.prot") fromto<-org.Rn.egENSEMBLPROT2EG
		if(initialIDs == "Ensembl.gene") fromto<-org.Rn.egENSEMBL2EG
		if(initialIDs == "RefSeq") fromto<-org.Rn.egREFSEQ2EG
		if(initialIDs == "Symbol") fromto<-org.Rn.egSYMBOL2EG
		if(initialIDs == "GenBank") fromto<-org.Rn.egACCNUM2EG
		if(initialIDs== "Entrez.gene") {
#If the initial identifiers is 	"Entrez.gene", then the mapping will 
#automatically be from Entrez Gene identifiers	to one of the following				
			if(finalIDs == "Ensembl.gene") fromto<-org.Rn.egENSEMBL
			if(finalIDs == "Ensembl.transcript") fromto<-org.Rn.egENSEMBLTRANS
			if(finalIDs == "Ensembl.prot") fromto<-org.Rn.egENSEMBLPROT
			if(finalIDs == "RefSeq") fromto<-org.Rn.egREFSEQ
			if(finalIDs == "Symbol") fromto<-org.Rn.egSYMBOL
			if(finalIDs == "GenBank") fromto<-org.Rn.egACCNUM
			}	
#Check that the environment has been correctly determined
		if(class(fromto) != "AnnDbBimap") {
			stop("Please provide a valid type of identifiers for the 'initialIDs' and 'finalIDs' parameters (see help(mammalAnnotationConvertor))")	
			}				
		}
#Determine if the geneList is a matrix (named rows) or a named vector 
#the computation will be done separately depending on this (i.e. based on names or rownames)					
	if(is.matrix(geneList)==FALSE){
#Check that the data vector is named 	
		if(is.null(names(geneList))) stop("Your data vector is not named")
#Create a list with an element for each name in the geneList, containing a vector 
#of identifiers of the type finalIDs mapped to that name in the geneList			
		list.new.names<-mget(names(geneList), fromto, ifnotfound=NA)
#Create a vector that will hold the new names, and a vector that will tag the names that were mapped to multiple identifiers		
		l.new.names<-length(list.new.names)
		new.names<-rep(0,l.new.names)
		tag.multiples<-rep(0,l.new.names)
#Go through the list of names and:
	#1. assign the first result in each element to the corresponding position in the new names vector
	#2. check if the element of the list contained more than one result
	#3. if the user asked to keep multiple mappings, just inform the user that this entry was mapped multiple times
	#4. if the user asked to discard multiple mappings, tag this entry and inform the user that this entry was mapped multiple times				
		for(i in 1:l.new.names){
			new.names[i]<-list.new.names[[i]][1]
			if(length(list.new.names[[i]]) > 1) {
				if(keepMultipleMappings == TRUE){
					print("The following identifier was mapped to more than one value (only the first value is kept): ")
					print(list.new.names[i])
					}
				if(keepMultipleMappings == FALSE){
					print("The following identifier was mapped to more than one value (this entry will be discarded): ")
					print(list.new.names[i])
					tag.multiples[i]<-1
					}
				}
			}
#If the user asked to keep multiple mappings, the vector of new names is used to set the names of the vector as it is
#and the entries that could not be mapped to any new names are removed from the data
#the user is informed of how many entries were removed because they were not mapped						
		if(keepMultipleMappings == TRUE){	
			newdata<-geneList	
			names(newdata)<-new.names
			newdata<-newdata[!is.na(names(newdata))]
			print(paste((length(geneList)-length(newdata))," genes (out of ",length(geneList) ,
				") could not be mapped to any identifier, and were removed from the data."))
			}
#If the user asked to discard multiple mappings, the data is trimmed of all entries that were mapped multiple times
#and so is the vector of new names.  This is done using the information in the tag.multiples vector					
		if(keepMultipleMappings == FALSE){	
			newdata<-geneList[which(tag.multiples == 0)]	
			names(newdata)<-new.names[which(tag.multiples == 0)]
			newdata<-newdata[!is.na(names(newdata))]
			print(paste((length(geneList)-length(newdata))," genes (out of ",length(geneList) ,
				") could not be mapped to any identifier (or were mapped to multiple identifiers), 
				and were removed from the data."))
			}
	}else{
#Check that the data matrix has row names	
		if(is.null(rownames(geneList))) stop("Your data matrix does not have row names")		
		list.new.names<-mget(rownames(geneList), fromto, ifnotfound=NA)
#This is identical to what is done for vectors, except that we work on row names		
		l.new.names<-length(list.new.names)
		new.names<-rep(0,l.new.names)
		tag.multiples<-rep(0,l.new.names)
		for(i in 1:l.new.names){
			new.names[i]<-list.new.names[[i]][1]
			if(length(list.new.names[[i]]) > 1) {
				if(keepMultipleMappings == TRUE){
					print("The following identifier was mapped to more than one value (only the first value is kept): ")
					print(list.new.names[i])
					}
				if(keepMultipleMappings == FALSE){
					print("The following identifier was mapped to more than one value (this entry will be discarded): ")
					print(list.new.names[i])
					tag.multiples[i]<-1
					}
				}
			}
#This is identical to what is done for vectors, except that we work on row names and that we discard the whole row of data when 
#there is no mapping				
		if(keepMultipleMappings == TRUE){	
			newdata<-geneList	
			rownames(newdata)<-new.names
			newdata<-newdata[!is.na(rownames(newdata)),]
			print(paste(((dim(geneList)[1])-(dim(newdata)[1]))," genes (out of ",dim(geneList)[1] ,
				") could not be mapped to any identifier, and were removed from the data."))
			}
		if(keepMultipleMappings == FALSE){	
			newdata<-geneList[which(tag.multiples == 0)]	
			rownames(newdata)<-new.names[which(tag.multiples == 0)]
			newdata<-newdata[!is.na(rownames(newdata)),]
			print(paste(((dim(geneList)[1])-(dim(newdata)[1]))," genes (out of ",dim(geneList)[1] ,
				") could not be mapped to any identifier (or were mapped to multiple identifiers), 
				and were removed from the data."))
			}
		}		
	return(newdata)	
	}

