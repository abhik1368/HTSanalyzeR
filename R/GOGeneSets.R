#This function creates gene set collections based on GO terms. 
#It is species-specific, and returns a GeneSetCollection objects 
#with the elements of the gene sets represented by Entrez Gene IDs.
GOGeneSets <-
	function(species="Dm", ontologies="MF") {
	#check input arguments
	paraCheck("species",species)
	paraCheck("ontologies",ontologies)
	#This creates a list with an element for each GO term, containing the Entrez Gene identifiers and corresponding evidence codes
	if(species == "Hs") ThisGO<-org.Hs.egGO2EG
	if(species == "Mm") ThisGO<-org.Mm.egGO2EG
	if(species == "Rn") ThisGO<-org.Rn.egGo2EG
	if(species == "Dm") ThisGO<-org.Dm.egGO2EG
	if(species == "Ce") ThisGO<-org.Ce.egGO2EG

	all.go.id<-GOID(GOTERM)
	names(all.go.id)<-NULL
	all.go.ontology<-Ontology(GOTERM)
	this.go.list<-as.list(ThisGO)
	this.go.list<-lapply(this.go.list,unique)
	this.go.id<-names(this.go.list)
	this.go.len<-length(this.go.id)
	junk<-sapply(1:this.go.len,function(s) names(this.go.list[[s]])<<-NULL)
	this.go.ontology.tag<-rep(FALSE,this.go.len)
	this.go.ontology.id<-unlist(sapply(1:length(ontologies), 
					function(n) {		
						match.id<-match(this.go.id,all.go.id)
						which(all.go.ontology[match.id]==ontologies[n])
					}))
	this.go.ontology.tag[this.go.ontology.id]<-TRUE
	this.go.list<-this.go.list[which(this.go.ontology.tag)]	
	return(this.go.list)	
}

