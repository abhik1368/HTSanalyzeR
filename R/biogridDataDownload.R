##This functions downloads the data from the Biogrid, puts it in the
##user-specified folder and extracts the interactions data for a given 
##species.

biogridDataDownload <- function(link, species = "Dm", dataDirectory = ".", 
		verbose = TRUE) {
	#check arguments
	if(!missing(link) && !is.null(link))
		paraCheck("link", link)
	else 
		link <- "http://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-3.1.71/BIOGRID-ORGANISM-3.1.71.tab2.zip"
	paraCheck("species", species)
	paraCheck("dataDirectory", dataDirectory)
	paraCheck("verbose",verbose)
	bionet.species <- c("Drosophila_melanogaster", "Homo_sapiens", 
		"Rattus_norvegicus", "Mus_musculus", "Caenorhabditis_elegans")
	names(bionet.species) <- c("Dm", "Hs", "Rn", "Mm", "Ce")
	#download the data from the BioGRID and put the compressed file 
	##into directory, then unzip the file
	if(verbose) 
		cat("--Downloading BioGRID interactome dataset ...\n")
	if(!file.exists(dataDirectory)) 
		dir.create(dataDirectory)
	download.file(url = link, destfile = file.path(dataDirectory, 
		"Biogrid-all-org"), quiet=TRUE)
	unzip(zipfile = file.path(dataDirectory, "Biogrid-all-org"), 
		exdir = dataDirectory)
	#the data directory now contains one tab delimited file for each 
	#species in the biogrid data the next few lines find and read the 
	#file that contains the data for the species that we want to use
	listfiles <- list.files(dataDirectory)
	biogridSpecies <- read.table(file = file.path(dataDirectory, 
		grep(bionet.species[species], listfiles, value = TRUE)), 
		header = TRUE, sep = "\t", fill = TRUE, comment.char = "")
	#Extract the relevant columns from the tab-delimited file that was read		
	source <- as.character(biogridSpecies[, "Entrez.Gene.Interactor.A"])
	target <- as.character(biogridSpecies[, "Entrez.Gene.Interactor.B"])
	interac.type <- as.character(biogridSpecies[, "Experimental.System.Type"])
	#Create a matrix from the data, and names its columns accordingly	
	interactions <- cbind(source, target, interac.type)
	colnames(interactions) <- c("InteractorA", "InteractorB", "InteractionType")
	return(interactions)
}

