##HTSanalyzeR pipeline for cellHTS2 objects
HTSanalyzeR4cellHTS2 <- function(
		normCellHTSobject,
		scoreSign = "-",
		scoreMethod = "zscore",
		summarizeMethod = "mean",
		annotationColumn = "GeneID",
		species = c("Dm"),
		initialIDs = "FlybaseCG",
		duplicateRemoverMethod = "max",
		orderAbsValue = FALSE,
		listOfGeneSetCollections,
		cutoffHitsEnrichment = 2,
		pValueCutoff = 0.05,
		pAdjustMethod = "BH",
		nPermutations = 1000,
		minGeneSetSize = 15,
		exponent = 1,
		keggGSCs,
		goGSCs,
		nwStatsControls = "neg",
		nwStatsAlternative = "two.sided",
		nwStatsTests = c("T-test"),
		nwStatsColumns = c("t.test.pvalues.two.samples", "t.test.pvalues.one.sample"),
		nwAnalysisFdr = 0.001,
		nwAnalysisGenetic = FALSE,
		interactionMatrix = NULL,
		nwAnalysisOrder = 2,
		ntop = NULL,
		allSig = TRUE,
		reportDir = "HTSanalyzerReport",
		verbose = TRUE
) {
	###############################
	##0. 	CellHTS2 part
	###############################
	##check that 'normCellHTSobject' is a cellHTS cellHTSobject
	if(class(normCellHTSobject) != "cellHTS") 
		stop("The argument 'normCellHTSobject' is not a cellHTS object")
	##check that the cellHTS2 cellHTSobject is in the right format 
	##(configured, so that we know which rows are samples and which are controls
	##annotated, so that we can match names of constructs normalized, but not scored)
	if(!state(normCellHTSobject)["configured"]) 
		stop("The cellHTS object should be configured to perform the statistical tests")
	if(!state(normCellHTSobject)["normalized"]) 
		warning("Your cellHTS object has not been normalized, this could impact the results of these tests", immediate. = TRUE)
	if(state(normCellHTSobject)["scored"]) 
		stop("This cellHTS object has been scored; the statistical analysis should be performed on the normalized signal intensities", immediate. = TRUE)
	if(!state(normCellHTSobject)["annotated"]) 
		stop("This cellHTS object has not been annotated", immediate. = TRUE)
	##check that the annotationColumn has been specified as a single character string
	if(!is.character(annotationColumn) || length(annotationColumn) !=1 ) 
		stop("The 'annotationColumn' parameter does not have the right format")
	##check that the annotationColumn matches a column in the fData(cellHTSobject) dataframe	
	if(!(annotationColumn %in% colnames(fData(normCellHTSobject)))) 
		stop("The 'annotationColumn' parameter does not match to any column in your cellHTS object")
	##Prepare the data for the enrichment analysis: 
	##1.we need a scored cellHTS2 object
	scoredCellHTSobject <- scoreReplicates(normCellHTSobject, 
		sign=scoreSign, method=scoreMethod)
	scoredCellHTSobject <- summarizeReplicates(scoredCellHTSobject, 
		summary=summarizeMethod)
	##2.we need a named vector, with no NA names
	data4enrich <- as.vector(Data(scoredCellHTSobject))
	names(data4enrich) <- fData(scoredCellHTSobject)[, annotationColumn] 
	data4enrich <- data4enrich[which(!is.na(names(data4enrich)) & 
		!is.na(data4enrich) )]
	paraCheck(cutoffHitsEnrichment)
	if(cutoffHitsEnrichment >= max(data4enrich) && 
		cutoffHitsEnrichment >= abs(min(data4enrich)))
		stop(paste("The cutoffHitsEnrichment argument should be comprised between 0 and ", 
			max(max(data4enrich),abs(min(data4enrich)))))
	if(cutoffHitsEnrichment <= 0)
		stop(paste("The cutoffHitsEnrichment argument should be comprised between 0 and ", 
			max(max(data4enrich),abs(min(data4enrich)))))
	################################
	##1. Gene set enrichment analysis
	################################
	##create a GSCA object 
	gsca <- new("GSCA", listOfGeneSetCollections = listOfGeneSetCollections,
		geneList = data4enrich,
		hits = names(data4enrich)[which(abs(data4enrich) > cutoffHitsEnrichment)])
	##preprocessing of input gene list and hit list * remove NA; 
	##duplicate operations; annotation conversions; order phenotypes
	gsca <- preprocess(gsca, species = species, initialIDs = initialIDs, 
		keepMultipleMappings = TRUE, duplicateRemoverMethod = "max", 
		orderAbsValue = FALSE, verbose = verbose)
	##do analysis
	gsca <- analyze(gsca, para = list(pValueCutoff = pValueCutoff, 
		pAdjustMethod = pAdjustMethod, nPermutations = nPermutations, 
		minGeneSetSize = minGeneSetSize,exponent = exponent), 
		verbose = verbose)
	################################
	##2. 	Network analysis
	################################	
	test.stats <- cellHTS2OutputStatTests(cellHTSobject = normCellHTSobject, 
		annotationColumn = annotationColumn, alternative = nwStatsAlternative, 
		controls = nwStatsControls, tests = nwStatsTests)
	##Check that the columns argument matches to column names in the pvaluesMatrix
	if(!all(nwStatsColumns %in% colnames(test.stats))) 
		stop("The names specified in the 'columns' argument do not match column names of your pvaluesMatrix")
	##Check that the pvaluesMatrix has row names
	if(is.null(rownames(test.stats))) 
		stop("The pvaluesMatrix should have row names (which should be Entrez Gene identifiers unless you provide a biogridObject with a different type of identifiers)")
	##Create a vector of pvalues or aggregated pvalues if mutliple columns are specified
	if(length(nwStatsColumns) == 1) {
		pvalues <- test.stats[, nwStatsColumns]
	} else {
		pvalues <- test.stats[, nwStatsColumns]
		pvalues <- aggrPvals(test.stats, order = nwAnalysisOrder, plot = FALSE)
	}
	##create a NWA (NetWork Analysis) object
	nwa <- new("NWA", pvalues = pvalues, phenotypes = data4enrich)
	##preprocessing
	nwa <- preprocess(nwa, species = species, initialIDs = initialIDs, 
		keepMultipleMappings = TRUE, duplicateRemoverMethod = duplicateRemoverMethod, 
		verbose = verbose)
	##create an interactome
	nwa <- interactome(nwa, interactionMatrix = interactionMatrix, 
		species = species, reportDir = reportDir, genetic = nwAnalysisGenetic, 
		verbose = verbose)
	##do analysis
	nwa <- analyze(nwa, fdr = nwAnalysisFdr, species = species, verbose = verbose)	
	################################
	##3. 	write report
	################################	
	reportAll(gsca = gsca, nwa = nwa, experimentName = 
		deparse(substitute(normCellHTSobject)), species = species, 
		ntop = ntop, allSig = allSig, keggGSCs = keggGSCs, 
		goGSCs = goGSCs, reportDir = reportDir)	
	##save GSCA and NWA objects to RData
	save(gsca, nwa, file = file.path(reportDir, 
		paste(deparse(substitute(normCellHTSobject)), ".RData", sep="")))
}
