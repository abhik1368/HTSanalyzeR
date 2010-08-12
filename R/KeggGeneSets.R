#This function creates gene set collections based on Kegg pathways terms. 
#It is species-specific, and returns a GeneSetCollection objects with 
#the elements of the gene sets represented by Entrez Gene IDs.
KeggGeneSets <-
function(species=c("Drosophila_melanogaster", "Homo_sapiens", "Rattus_norvegicus", "Mus_musculus", "Caenorhabditis_elegans")){
#This checks that the species argument corresponds to a single character string that can be matched to one of our Kegg species codes
if(is.character(species) != TRUE | length(species) != 1) stop("The species argument does not have the right format")
if(is.na(match(species,c("Drosophila_melanogaster","Homo_sapiens","Rattus_norvegicus","Mus_musculus","Caenorhabditis_elegans")))) {
	stop("The species argument does not match any of the names recognized by this function, please provide one of the following character strings: 'Drosophila_melanogaster','Homo_sapiens','Rattus_norvegicus','Mus_musculus','Caenorhabditis_elegans'")
	}
#Create a list with an element for each pathway, each element containing a vector of gene identifiers	
	Kegg <- as.list(KEGGPATHID2EXTID)
#Keep only those elements of the list that correspond to pathways from the gievn species	
	if(species=="Homo_sapiens") Kegg<-Kegg[grep("hsa",names(Kegg))]
	if(species=="Mus_musculus") Kegg<-Kegg[grep("mmu",names(Kegg))]
	if(species=="Rattus_norvegicus") Kegg<-Kegg[grep("rno",names(Kegg))]
	if(species=="Caenorhabditis_elegans") Kegg<-Kegg[grep("cel",names(Kegg))]
#The drosophila pathways are slightly more complicated because the external IDs are flybase gene IDs instead of Entrez Gene IDs
#So we first extract those pathways, 
#then create a list that maps from flybase to EG, 
#then go through each element of the vector in each element of the list, and keep the first EG ID that maps to 
#that particular flybase gene ID
	if(species=="Drosophila_melanogaster") {
		Kegg<-Kegg[grep("dme",names(Kegg))]
		fbgn2eg <- as.list(org.Dm.egFLYBASE2EG)
		for(i in 1:length(Kegg)){
			for(n in 1:length(Kegg[[i]])){
				if(length(fbgn2eg[[Kegg[[i]][n]]]) != 0){Kegg[[i]][n]<-fbgn2eg[[Kegg[[i]][n]]][1]}
				}
			}
		}
#Create the gene set collection		
	gsc.kegg <- GeneSetCollection(mapply(function(geneIds, keggId) {GeneSet(unique(geneIds),
		geneIdType=EntrezIdentifier(),collectionType=KEGGCollection(keggId),setName=keggId)}, Kegg, names(Kegg)))
	return(gsc.kegg)
	}

