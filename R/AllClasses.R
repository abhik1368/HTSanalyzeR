##Class GSCA (Gene Set Collection Analyses)
##definition of class gene set collection analysis
setClass(
		"GSCA",
		representation(
				listOfGeneSetCollections="list",
				geneList="numeric_Or_integer",
				hits="character",
				para="list",
				result="list",
				summary="list",
				preprocessed="logical"
		),
		prototype=list(
				listOfGeneSetCollections=list(),
				geneList=numeric(),
				hits=character(),
				para=list(pValueCutoff = 0.05,pAdjustMethod = "BH", nPermutations = 1000, minGeneSetSize = 15,exponent = 1),
				result=list(),
				summary=list(),
				preprocessed=FALSE
		)
)

##class NWA (NetWork Analyses)
##definition of class NWA
setClass(
		"NWA",
		representation(
				pvalues = "numeric",
				phenotypes = "numeric_Or_integer_Or_NULL",
				interactome = "graphNEL_Or_NULL",
				fdr = "numeric",
				result = "list",
				summary = "list",
				preprocessed = "logical"
		),
		prototype = list(
				pvalues = numeric(),
				phenotypes = NULL,
				interactome = NULL,
				fdr = 0.001,
				result = list(),
				summary = list(),
				preprocessed = FALSE
		)
)