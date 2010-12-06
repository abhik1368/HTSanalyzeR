##This function converts an initial named data vector to the same
##vector but with a different identifier category. This function can
##also take a matrix, with rows=gene id's. This function removes the
##genes for which no mapping were found.

annotationConvertor<-function(geneList, species="Dm", 
	initialIDs="Entrez.gene", finalIDs="Entrez.gene", 
	keepMultipleMappings=TRUE, verbose=TRUE) {
	
	##check arguments
	paraCheck("genelist.general", geneList)

	paraCheck("keepMultipleMappings", keepMultipleMappings)
	paraCheck("verbose", verbose)
	
	##convert annotations if initialIDs are not Entrez 
	if(initialIDs != "Entrez.gene") {
		##check species
		paraCheck(name = "species", para = species)	
		if(species == "Dm") {
			paraCheck("dro.initialIDs", initialIDs)
			paraCheck("dro.finalIDs", finalIDs)
			if(!("package:org.Dm.eg.db" %in% search())) 
				library(org.Dm.eg.db)
			geneListEntrez <- drosoAnnotationConvertor(
				geneList = geneList, initialIDs = initialIDs, 
				finalIDs = finalIDs, verbose = verbose)
		} else if(species == "Hs") {
			paraCheck("mam.initialIDs", initialIDs)
			paraCheck("mam.finalIDs", finalIDs)
			if(!("package:org.Hs.eg.db" %in% search())) 
				library(org.Hs.eg.db)
			geneListEntrez <- mammalAnnotationConvertor( 
				geneList = geneList, initialIDs = initialIDs, 
				finalIDs = finalIDs, species = species, verbose = verbose)
		} else if(species == "Rn") {
			paraCheck("mam.initialIDs", initialIDs)
			paraCheck("mam.finalIDs", finalIDs)
			if(!("package:org.Rn.eg.db" %in% search())) 
				library(org.Rn.eg.db)
			geneListEntrez <- mammalAnnotationConvertor(
				geneList = geneList, initialIDs = initialIDs, 
				finalIDs = finalIDs, species = species, verbose = verbose)
		} else if(species == "Mm") {
			paraCheck("mam.initialIDs", initialIDs)
			paraCheck("mam.finalIDs", finalIDs)
			if(!("package:org.Mm.eg.db" %in% search())) 
				library(org.Mm.eg.db)
			geneListEntrez <- mammalAnnotationConvertor(
				geneList = geneList, initialIDs = initialIDs, 
				finalIDs = finalIDs, species = species, verbose = verbose)
		} else if(species == "Ce") {
			paraCheck("cel.initialIDs", initialIDs)
			paraCheck("cel.finalIDs", finalIDs)
			if(!("package:org.Ce.eg.db" %in% search())) 
				library(org.Ce.eg.db)
			geneListEntrez <- celAnnotationConvertor(
				geneList = geneList, initialIDs = initialIDs , 
				finalIDs = finalIDs, verbose = verbose)
		}
	} else {
		geneListEntrez<-geneList
		names(geneListEntrez)<-names(geneList)
	}
	return(geneListEntrez)
}

