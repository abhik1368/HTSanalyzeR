#This function creates gene set collections based on Kegg pathways terms. 
#It is species-specific, and returns a GeneSetCollection objects with 
#the elements of the gene sets represented by Entrez Gene IDs.
KeggGeneSets <- function(species = "Dm") {
		paraCheck("species", species)
		##Create a list with an element for each pathway, each element 
		##containing a vector of gene identifiers	
		Kegg <- as.list(KEGGPATHID2EXTID)
		##Keep only those elements of the list that correspond to 
		##pathways from the gievn species	
		if(species == "Hs") 
			this.pw.kegg <- Kegg[grep("hsa", names(Kegg))]
		else if(species=="Mm") 
			this.pw.kegg <- Kegg[grep("mmu", names(Kegg))]
		else if(species=="Rn") 
			this.pw.kegg <- Kegg[grep("rno", names(Kegg))]
		else if(species=="Ce") 
			this.pw.kegg <- Kegg[grep("cel", names(Kegg))]
		##The drosophila pathways are slightly more complicated because 
		##the external IDs are flybase gene IDs instead of Entrez Gene IDs
		##So we first extract those pathways, 
		##then create a list that maps from flybase to EG, 
		##then go through each element of the vector in each element of 
		##the list, and keep the first EG ID that maps to that particular 
		##flybase gene ID
		if(species == "Dm") {
			this.pw.kegg <- Kegg[grep("dme", names(Kegg))]
			fbgn2eg <- as.list(org.Dm.egFLYBASE2EG)
			junk <- sapply(1:length(this.pw.kegg), 
				function(l) {
					this.pw.kegg[[l]] <<- 
						unlist(sapply(this.pw.kegg[[l]], function(i) {
							if(length(fbgn2eg[[i]])>0) 
								fbgn2eg[[i]][1]
					}))
				}
			)
		}
		return(this.pw.kegg)
}

