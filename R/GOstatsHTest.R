GOstatsHTest <-
function(hits,GeneList,annotation,ontology=c("BP","CC","MF"),pvalueCutoff,conditional=FALSE,testDirection="over") {
	##Implements hyperGTest from GOstats package
	##creates GOHyperGOParams object within function
	
	params<-new("GOHyperGParams",geneIds=hits,universeGeneIds=GeneList,
		annotation=annotation,ontology=ontology,pvalueCutoff=pvalueCutoff,
		conditional=conditional,testDirection=testDirection)
		
		hgtResults<-hyperGTest(params)
		hgtResults
	}

