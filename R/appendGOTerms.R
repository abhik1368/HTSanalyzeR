appendGOTerms <-
function(resultsDataFrame){
	#This function appends the corresponding GO term to GO ID 
	#in the analyzeGeneSetCollection results data frames.
	
	#Gene set names must have been supplied in format "GO:0000000"

	GOIds<-rownames(resultsDataFrame)
	if (nchar(GOIds[1]==10)) substr(GOIds,1,3)<-"GO:" ## ensures proper prefix (hyperGeoTest 
													  ##outputs them as "GO." for unknown reason)
	GONames<-Term(GOIds)
	GOontols<-Ontology(GOIds)
	new.names<-cbind(GOIds,GONames,"(",GOontols,")")
	new.names<-paste(new.names[,1],"::",new.names[,2],
		"(",GOontols,")")
	new.names
}

