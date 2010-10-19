appendKEGGNames <-
function(resultsDataFrame){
	#This function appends the corresponding KEGG term to KEGG Path ID 
	#in the analyzeGeneSetCollection results data frames.
	#Gene set names must have been supplied in xxx00000 or 00000 format,
	#where xxx is the organism code.
	KEGGIds<-rownames(resultsDataFrame)
	KEGGNumIds<-sapply(KEGGIds, 
			function(x) {
				if (nchar(x)==8) substr(x,4,nchar(x)) ##removes organism code
			}
	)
	KEGGNames<-mget(KEGGNumIds, env=KEGGPATHID2NAME)
	KEGGNamesVec<-unlist(KEGGNames)
	GeneSetNames<-cbind(KEGGIds,KEGGNamesVec)
	GeneSetNames<-paste(GeneSetNames[,1],"::",GeneSetNames[,2])
	GeneSetNames
}

