###############################################################################
# convert gene annotations between species
# Xin Wang <xw264@cam.ac.uk>
# Advisor: Florian Markowetz <florian.markowetz@cancer.org.uk> 
# University of Cambridge Deparment of Oncology
# Cancer Research UK - Cambridge Research Institute
# At 12:12:31, on 17 Nov 2010
###############################################################################
annotationConvertor<-function(geneList, species="Dm", initialIDs="Entrez.gene", finalIDs="Entrez.gene", keepMultipleMappings=TRUE, verbose=TRUE) {
	paraCheck("genelist.general",geneList)
	paraCheck("initialIDs",initialIDs)
	paraCheck("finalIDs",finalIDs)
	paraCheck("keepMultipleMappings",keepMultipleMappings)
	paraCheck("verbose",verbose)
	if(initialIDs != "Entrez.gene") {
		#This checks that the species argument corresponds to a single character string 
		#that can be matched to one of the annotationConvertor functions
		paraCheck(name="species",para=species)
		#The annotationConvertor functions will check that the initialIDs argument has the right format			
		if(species == "Dm") {
			if(!("package:org.Dm.eg.db" %in% search())) library(org.Dm.eg.db)
			geneListEntrez<-drosoAnnotationConvertor(geneList=geneList,initialIDs=initialIDs,finalIDs=finalIDs, verbose=verbose)
		}
		if(species == "Hs") {
			if(!("package:org.Hs.eg.db" %in% search())) library(org.Hs.eg.db)
			geneListEntrez<-mammalAnnotationConvertor(geneList=geneList,initialIDs=initialIDs,finalIDs=finalIDs,species=species, verbose=verbose)
		}
		if(species == "Rn") {
			if(!("package:org.Rn.eg.db" %in% search())) library(org.Rn.eg.db)
			geneListEntrez<-mammalAnnotationConvertor(geneList=geneList,initialIDs=initialIDs,finalIDs=finalIDs,species=species, verbose=verbose)
		}
		if(species == "Mm") {
			if(!("package:org.Mm.eg.db" %in% search())) library(org.Mm.eg.db)
			geneListEntrez<-mammalAnnotationConvertor(geneList=geneList,initialIDs=initialIDs,finalIDs=finalIDs,species=species, verbose=verbose)
		}
		if(species == "Ce") {
			if(!("package:org.Ce.eg.db" %in% search())) library(org.Ce.eg.db)
			geneListEntrez<-celAnnotationConvertor(geneList=geneList,initialIDs=initialIDs,finalIDs=finalIDs, verbose=verbose)
		}
	} else {
		geneListEntrez<-geneList
		names(geneListEntrez)<-names(geneList)
	}
	
	return(geneListEntrez)
}

