#This function downloads the data from The BioGRID (all organisms) and extracts the data for a chosen species (specified by the argument "species") into a matrix
#with one lign for each interaction, containing the Entrez identifiers for the interaction partners and the type of interaction (physical or genetic)
biogridDataDownload <-
function(link="http://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-3.0.64/BIOGRID-ORGANISM-3.0.64.tab2.zip",
	species=c("Dm","Hs","Rn","Mm","Ce"),dataDirectory,verbose=TRUE) {
	#This checks that the species argument corresponds to a single character string that can be matched to the names of
	#one of the BioGRID data files	
	if(!(is.character(species) && length(species) == 1)) 
		stop("The species argument does not have the right format")
	if(!(species %in% c("Dm","Hs","Rn","Mm","Ce"))) {
		stop("The species argument does not match any of the names recognized by this function, please provide one of the following 
			character strings: 'Dm' ('Drosophila_melanogaster'), 'Hs' ('Homo_sapiens'), 'Rn' ('Rattus_norvegicus'), 'Mm' ('Mus_musculus'), 'Ce' ('Caenorhabditis_elegans')")
	}
	bionet.species<-c("Drosophila_melanogaster", "Homo_sapiens","Rattus_norvegicus","Mus_musculus","Caenorhabditis_elegans")
	names(bionet.species)<-c("Dm","Hs","Rn","Mm","Ce")
	#download the data from the BioGRID and put the compressed file into the above directory, then unzip the file into the directory created above	
	if(verbose) 
		cat("--Downloading BioGRID interactome dataset \n")
	if(!file.exists(dataDirectory)) 
		dir.create(dataDirectory)
	download.file(url=link,destfile=file.path(dataDirectory,"Biogrid-all-org"),quiet=TRUE)
	unzip(zipfile=file.path(dataDirectory,"Biogrid-all-org"),exdir=dataDirectory)
	#the data directory now contains one tab delimited file for each species in the biogrid data
	#the next few lines find and read the file that contains the data for the species that we want to use
	listfiles<-list.files(dataDirectory)
	biogridSpecies<-read.table(file=file.path(dataDirectory,grep(bionet.species[species],listfiles,value=TRUE)),header=TRUE,sep="\t",fill=TRUE,comment.char="")
	#the following code extracts the relevant columns from the tab-delimited file that was read		
	source<-as.character(biogridSpecies[,"Entrez.Gene.Interactor.A"])
	target<-as.character(biogridSpecies[,"Entrez.Gene.Interactor.B"])
	interac.type<-as.character(biogridSpecies[,"Experimental.System.Type"])
	#this creates a matrix from the data, and names its columns accordingly	
	interactions<-cbind(source,target,interac.type)
	colnames(interactions)<-c("InteractorA","InteractorB","InteractionType")
	return(interactions)
}

