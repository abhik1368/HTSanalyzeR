getGeneSetInfo <-
	function(GeneList,GeneSet,exponent=1,species=c("Dm","Hs","Rn","Mm","Ce")) {
	
	##Returns the EntrezID, symbol, numeric phenotype, the phenotype rank,
	## and whenther the gene is in the leading edge subset.
	
	##calculate running enrichment scores
	gseaResults<-gseaScores(geneList=GeneList,geneSet=GeneSet, exponent=exponent,mode="graph")
	ind<-which(gseaResults$Positions==1) 	#which positions in ranked GeneList are in GeneSet
	Genes.In.Universe<-names(GeneList[ind])
	
	##Get gene symbols
	if(species == "Hs") Symbol<-mget(Genes.In.Universe, env=org.Hs.egSYMBOL)
	if(species == "Dm") Symbol<-mget(Genes.In.Universe, env=org.Dm.egSYMBOL)
	if(species == "Rn") Symbol<-mget(Genes.In.Universe, env=org.Rn.egSYMBOL)
	if(species == "Mm") Symbol<-mget(Genes.In.Universe, env=org.Mm.egSYMBOL)
	if(species == "Ce") Symbol<-mget(Genes.In.Universe, env=org.Ce.egSYMBOL)
	
	
	Phenotype<-GeneList[ind]
	Phenotype.Rank<-ind
	highscore<-which(gseaResults$Running.Score==gseaResults$Enrichment.Score)
	if(gseaResults$Enrichment.Score > 0) {
		peak<-which(Genes.In.Universe==names(GeneList[highscore]))
	} else {
		peak<-which(Genes.In.Universe==names(GeneList[(highscore+1)]))
	}
	if(gseaResults$Enrichment.Score > 0) {
		Leading.Edge<-c(rep(1,peak),rep(0,length(Genes.In.Universe)-peak))
	} else {
		Leading.Edge<-c(rep(0,(peak-1)),rep(1,length(Genes.In.Universe)-peak+1))
	}
	EntrezID<-Genes.In.Universe
	results<-cbind(EntrezID,Symbol, Phenotype,Phenotype.Rank,Leading.Edge)
	results<-as.data.frame(results)
	results
}

