HTSanalysis4CellHTS2 <- function(
		normCellHTSobject,
		scoreSign="-",
		scoreMethod="zscore",
		summarizeMethod="mean",
		annotationColumn="GeneID",
		species=c("Dm"),
		initialIDs="FlybaseCG",
		duplicateRemoverMethod="max",
		orderAbsValue=FALSE,
		listOfGeneSetCollections,
		cutoffHitsEnrichment=2,
		pAdjMethod="BH",
		nPermutations=1000,
		minGeneSetSize=15,
		exponent=1,
		whichSetIsKEGGIds="none",
		whichSetIsGOIds="none",
		nwStatsControls="neg",
		nwStatsAlternative="two.sided",
		nwStatsTests=c("T-test"),
		nwStatsColumns=c("t.test.pvalues.two.samples","t.test.pvalues.one.sample"),
		nwAnalysisFdr=0.001,
		nwAnalysisGenetic=FALSE,
		networkObject=NA,
		nwAnalysisOrder=2,
		nGseaPlots=10,
		reportdir="HTSanalyzerReport",
		verbose=TRUE
	) {
	#check that 'normCellHTSobject' is a cellHTS cellHTSobject
	if(class(normCellHTSobject) != "cellHTS") 
		stop("The argument 'normCellHTSobject' is not a cellHTS object")
	#check that the cellHTS2 cellHTSobject is in the right format 
	#(configured, so that we know which rows are samples and which are controls
	#annotated, so that we can match names of constructs
	#normalized, but not scored)
	if(!state(normCellHTSobject)["configured"]) 
		stop("The cellHTS object should be configured to perform the statistical tests")
	if(!state(normCellHTSobject)["normalized"]) 
		warning("Your cellHTS object has not been normalized, this could impact the results of these tests",immediate.=TRUE)
	if(state(normCellHTSobject)["scored"]) 
		stop("This cellHTS object has been scored; the statistical analysis should be performed on the normalized signal intensities",immediate.=TRUE)
	if(!state(normCellHTSobject)["annotated"]) 
		stop("This cellHTS object has not been annotated",immediate.=TRUE)
	#check that the annotationColumn has been specified as a single character string
	if(!is.character(annotationColumn) || length(annotationColumn) !=1 ) 
		stop("The 'annotationColumn' parameter does not have the right format")
	#check that the annotationColumn matches a column in the fData(cellHTSobject) dataframe	
	if(!(annotationColumn %in% colnames(fData(normCellHTSobject)))) 
		stop("The 'annotationColumn' parameter does not match to any column in your cellHTS object")
	#Prepare the data for the enrichment analysis: 
	#1.we need a scored cellHTS2 object
	scoredCellHTSobject<-scoreReplicates(normCellHTSobject,sign=scoreSign,method=scoreMethod)
	scoredCellHTSobject<-summarizeReplicates(scoredCellHTSobject,summary=summarizeMethod)
	#2.we need a named vector, with no NA names
	data4enrich<-as.vector(Data(scoredCellHTSobject));
	names(data4enrich)<-fData(scoredCellHTSobject)[,annotationColumn] 
	data4enrich<-data4enrich[-which(is.na(names(data4enrich)))]
	#3.we need Entrez.gene identifiers
	if(initialIDs != "Entrez.gene") {
	#This checks that the species argument corresponds to a single character string 
	#that can be matched to one of the annotationConvertor functions
		if(!is.character(species) || length(species) != 1) 
			stop("The species argument does not have the right format")
		if(!(species %in% c("Dm","Hs","Rn","Mm","Ce"))) {
			stop("The species argument does not match any of the names recognized by this function, please provide one of the following character strings: 
				'Dm' ('Drosophila_melanogaster'), 'Hs' ('Homo_sapiens'), 'Rn' ('Rattus_norvegicus'), 'Mm' ('Mus_musculus'), 'Ce' ('Caenorhabditis_elegans')")
			}
		#The AnnotationConvertor functions will check that the initialIDs argument has the right format			
		if(species == "Dm") {
			library(org.Dm.eg.db)
			data4enrichentrez<-drosoAnnotationConvertor(geneList=data4enrich,initialIDs=initialIDs)
		}
		if(species == "Hs") {
			library(org.Hs.eg.db)
			data4enrichentrez<-mammalAnnotationConvertor(geneList=data4enrich,initialIDs=initialIDs,species=species)
		}
		if(species == "Rn") {
			library(org.Rn.eg.db)
			data4enrichentrez<-mammalAnnotationConvertor(geneList=data4enrich,initialIDs=initialIDs,species=species)
		}
		if(species == "Mm") {
			library(org.Mm.eg.db)
			data4enrichentrez<-mammalAnnotationConvertor(geneList=data4enrich,initialIDs=initialIDs,species=species)
		}
		if(species == "Ce") {
			library(org.Ce.eg.db)
			data4enrichentrez<-celAnnotationConvertor(geneList=data4enrich,initialIDs=initialIDs)
		}
	} else {
		data4enrichentrez<-data4enrich
		names(data4enrichentrez)<-names(data4enrich)
	}
	data4enrichentrez<-data4enrichentrez[-which(is.na(data4enrichentrez))]
	#4.we need 	to remove the duplicates, and order the gene list
	#This function will check that the method is correctly specified
	data4enrichentrez<-duplicateRemover(geneList=data4enrichentrez,method=duplicateRemoverMethod,absValue=orderAbsValue)
	#Now we can perform the enrichment analysis
	#check that cutoffHitsEnrichment is a single integer, 
	#and is comprised in the range of the data
	#this cutoff is taken in terms of absolute value
	if(length(cutoffHitsEnrichment) != 1 || (!is.integer(cutoffHitsEnrichment) && !is.numeric(cutoffHitsEnrichment))) 
		stop("The cutoffHitsEnrichment should be a single integer ")
	if(cutoffHitsEnrichment >= max(data4enrichentrez) && cutoffHitsEnrichment >= abs(min(data4enrichentrez))) {
		stop(paste("The cutoffHitsEnrichment argument should be comprised between 0 and ",max(max(data4enrichentrez),abs(min(data4enrichentrez)))))
	}
	if(cutoffHitsEnrichment <= 0) {
		stop(paste("The cutoffHitsEnrichment argument should be comprised between 0 and ",max(max(data4enrichentrez),abs(min(data4enrichentrez)))))
	}	
	if(whichSetIsKEGGIds == "none") {
		if(length(grep("kegg",ignore.case=FALSE,names(listOfGeneSetCollections))) != 0) 
			whichSetIsKEGGIds<-grep("kegg",ignore.case=FALSE,names(listOfGeneSetCollections))
	}
	if(whichSetIsGOIds == "none") {
		if(length(grep("GO",ignore.case=FALSE,names(listOfGeneSetCollections))) != 0)
			whichSetIsGOIds<-grep("GO",ignore.case=FALSE,names(listOfGeneSetCollections))
	}		
	enrichment.analysis<-analyzeGeneSetCollections(
		ListOfGeneSetCollections=listOfGeneSetCollections,
		GeneList=data4enrichentrez,
		hits=names(data4enrichentrez)[which(abs(data4enrichentrez)>cutoffHitsEnrichment)],
		pAdjustMethod=pAdjMethod,
		npermutations=nPermutations,
		min.gene.set.size=minGeneSetSize,
		exponent=exponent,
		whichSetIsKEGGIds=whichSetIsKEGGIds,
		whichSetIsGOIds=whichSetIsGOIds,
		verbose=verbose
	)
	#Now we can perform the network analysis
	#This function will test that the controls argument is correctly set
	#that the alternative, columns, tests arguments are correctly set
	test.module<-networkAnalysis(
		cellHTSobject=normCellHTSobject,
		annotationColumn=annotationColumn,
		controls=nwStatsControls,
		alternative=nwStatsAlternative,
		logged=FALSE,
		tests=nwStatsTests,
		columns=nwStatsColumns,
		species=species,
		initialIDs=initialIDs,
		fdr=nwAnalysisFdr,
		genetic=nwAnalysisGenetic,
		biogridObject=networkObject,
		#phenotypeVector=data4enrichentrez,
		order=nwAnalysisOrder,
		reportdir=reportdir,
		verbose=verbose
	)			
	#Now we can write the report
	writeReportHTSA(
		experimentName=deparse(substitute(normCellHTSobject)),
		enrichmentAnalysis=enrichment.analysis,
		cutoffHits=cutoffHitsEnrichment,
		hits=names(data4enrichentrez)[which(abs(data4enrichentrez)>cutoffHitsEnrichment)],
		listOfGeneSetCollections=listOfGeneSetCollections,
		geneList=data4enrichentrez,
		p.adj.method=pAdjMethod,
		nPermutations=nPermutations,
		min.gene.set.size=minGeneSetSize,
		exponent=exponent,
		nwAnalysisOutput=test.module,
		nwAnalysisGraphFile="EnrichedSubNw.png",
		controls=nwStatsControls,
		alternative=nwStatsAlternative,
		tests=nwStatsTests,
		columns=nwStatsColumns,
		species=species,
		fdr=nwAnalysisFdr,
		genetic=nwAnalysisGenetic,
		networkObject=networkObject,
		nGseaPlots=nGseaPlots,
		geneListName=paste("scored (",scoreMethod,")",deparse(substitute(normCellHTSobject))),
		whichSetIsKEGGIds=whichSetIsKEGGIds,
		whichSetIsGOIds=whichSetIsGOIds,
		reportdir=reportdir
	)
}