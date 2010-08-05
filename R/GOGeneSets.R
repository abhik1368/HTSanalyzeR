#This function creates gene set collections based on GO terms. 
#It is species-specific, and returns a GeneSetCollection objects 
#with the elements of the gene sets represented by Entrez Gene IDs.
GOGeneSets <-
function(species=c("Drosophila_melanogaster", "Homo_sapiens", "Rattus_norvegicus", "Mus_musculus", "Caenorhabditis_elegans"),
	ontologies=c("BP","MF","CC")){
#This checks that the species argument corresponds to a single character string that can be matched
	if(is.character(species) != TRUE | length(species) != 1) stop("The species argument does not have the right format")
	if(is.na(match(species,c("Drosophila_melanogaster","Homo_sapiens","Rattus_norvegicus","Mus_musculus","Caenorhabditis_elegans")))) {
		stop("The species argument does not match any of the names recognized by this function, please provide one of the following character strings: 'Drosophila_melanogaster','Homo_sapiens','Rattus_norvegicus','Mus_musculus','Caenorhabditis_elegans'")
		}
#This checks that the ontologies argument is appropriate
	if(length(intersect(ontologies,c("BP","MF","CC"))) != length(ontologies)) stop('The "ontologies" argument is not corectly specified, please provide a vector containing any non redundant combination of "BP","MF","CC" ')
#This creates a list with an element for each GO term, containing the Entrez Gene identifiers and corresponding evidence codes
	if(species == "Homo_sapiens") lstGO<-as.list(org.Hs.egGO2EG)
	if(species == "Mus_musculus") lstGO<-as.list(org.Mm.egGO2EG)
	if(species == "Rattus_norvegicus") lstGO<-as.list(org.Rn.egGO2EG)
	if(species == "Drosophila_melanogaster") lstGO<-as.list(org.Dm.egGO2EG)
	if(species == "Caenorhabditis_elegans") lstGO<-as.list(org.Ce.egGO2EG)
#This creates a list of all GO terms, and a vector that will be used to tag the terms from the lstGO list that belong to the right ontology	
	goList<-as.list(GOTERM)
	l.lstGO<-length(lstGO)
	tag.ontology<-rep(0,l.lstGO)
#Go through each ontology specified in the "ontologies" parameter	
	for(n in 1:length(ontologies)){
#Go through each element of the lstGO list and keep only one record of each Entrez identifier (there can be multiples when the evidence codes are different)
#Remove also the names of each entry (evidence codes) because otherwise the GeneSetCollection function does not work
#Tag the ontology if it is in the ontologies parameter	
		for(i in 1:l.lstGO){
			lstGO[[i]]<-unique(lstGO[[i]])
			names(lstGO[[i]])<-c()
			if(Ontology(goList[[names(lstGO)[i]]]) == ontologies[n]){tag.ontology[i]=1}
			}	
		}
#keep only those elements of the lstGO that have been tagged as belonging to the right ontology		
	lstGO<-lstGO[which(tag.ontology == 1)]	
#create and return a GeneSetCollection object	
	gsc <- GeneSetCollection(mapply(function(geneIds, GOId) {
		GeneSet(geneIds,geneIdType=EntrezIdentifier(),setName=GOId)
		}, lstGO, names(lstGO)))	
	return(gsc)
	}

