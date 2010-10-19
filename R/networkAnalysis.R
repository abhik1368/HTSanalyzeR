#This function finds subnetworks enriched for genes with significant phenotypes.  
#It is a wrapper around the whole network analysis process in this package.
networkAnalysis <-
function(
	cellHTSobject,
	annotationColumn="GeneID",
	controls="neg",
	alternative=c("two.sided","less","greater"),
	logged=FALSE,
	tests=c("T-test"),
	columns=c("t.test.pvalues.two.samples","t.test.pvalues.one.sample"),
	species=c("Dm"),
	initialIDs="FlybaseCG",
	fdr=0.001,
	genetic=FALSE,
	biogridObject=NA,
	order=2,
	link="http://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-3.0.64/BIOGRID-ORGANISM-3.0.64.tab2.zip",
	reportdir="HTSanalyzerReport",
	verbose=TRUE
) {
	#Do the statistical tests: this function will test:
	#-that the cellHTSobject is a cellHTS object
	#that the cellHTS2 cellHTSobject is in the right format (configured, so that we know which rows are samples and which are controls)
	#that the annotationColumn has been specified as a single character string
	#that the annotationColumn matches a column in the fData(cellHTSobject) dataframe
	#that the 'controls' parameter has been specified as a single character string
	#that the  'controls' parameter matches a status in the 'controlStatus' column of the cellHTS cellHTSobject
	#that the 'alternative' parameter has been correctly set 	
	cat("-Performing network analysis ... \n")
	test.stats<-cellHTS2OutputStatTests(cellHTSobject=cellHTSobject,alternative=alternative,tests=tests)
	#Check that the species is correctly specified (corresponds to a single character string 
	#that can be matched to an input for the AnnotationConvertor functions)
	if(!is.character(species) || length(species) != 1) 
		stop("The species argument does not have the right format")
	if(!(species %in% c("Dm","Hs","Rn","Mm","Ce")))
		stop("The species argument does not match any of the names recognized by this function, please provide one of the following character strings: 
			\"Dm\" for Drosophila_melanogaster,\"Hs\" for Homo_sapiens,\"Rn\" for Rattus_norvegicus, \"Mm\" for Mus_musculus, \"Ce\" for Caenorhabditis_elegans")
	#The biogridObject contains Entrez.gene IDs, so if this is not the initial type of identifiers,
	#they need to be mapped to Entrez.gene IDs		
	if(initialIDs != "Entrez.gene") {
		#Check that the initialIDs argument is correctly specified
		if(!(initialIDs %in% c("Ensembl.transcript", "Ensembl.prot", "Ensembl.gene", "Entrez.gene", "RefSeq", "Symbol", "GenBank", "Flybase", "FlybaseCG", "FlybaseProt"))) 
			stop('The initialIDs argument is not correctly specified, please provide one of the following character strings: "Ensembl.transcript", "Ensembl.prot", "Ensembl.gene", "Entrez.gene", "RefSeq", "Symbol", "GenBank", "Flybase", "FlybaseCG", "FlybaseProt"')
		#map the identifiers to Entrez.gene identifiers			
		if(species == "Dm") test.stats.entrez<-drosoAnnotationConvertor(geneList=test.stats,initialIDs=initialIDs,finalIDs="Entrez.gene",verbose=verbose)
		if(species == "Ce") test.stats.entrez<-celAnnotationConvertor(geneList=test.stats,initialIDs=initialIDs,finalIDs="Entrez.gene",verbose=verbose)
		if(species == "Hs") test.stats.entrez<-mammalAnnotationConvertor(geneList=test.stats,species=species,initialIDs=initialIDs,finalIDs="Entrez.gene",verbose=verbose)
		if(species == "Rn") test.stats.entrez<-mammalAnnotationConvertor(geneList=test.stats,species=species,initialIDs=initialIDs,finalIDs="Entrez.gene",verbose=verbose)
		if(species == "Mm") test.stats.entrez<-mammalAnnotationConvertor(geneList=test.stats,species=species,initialIDs=initialIDs,finalIDs="Entrez.gene",verbose=verbose)
		}
	#create folders for biogrid date downloading
	biogridDataDir=file.path(reportdir,"Data")
	if(!file.exists(reportdir)) 
		dir.create(reportdir)
	subnw<-enrichedSubNw(
		pvaluesMatrix=test.stats.entrez,
		columns=columns,
		species=species,
		fdr=fdr,genetic=genetic,
		biogridObject=biogridObject,order=order,biogridDataDir=biogridDataDir,verbose=verbose
	)
	nodesModule<-nodes(subnw)
	#To represent the network in a more convenient format, the symbol identifiers will be mapped and given to the user (more readable than Entrez.gene IDs)	
	if(species == "Dm") map <- org.Dm.egSYMBOL
	if(species == "Ce") map <- org.Ce.egSYMBOL
	if(species == "Hs") map <- org.Hs.egSYMBOL
	if(species == "Rn") map <- org.Rn.egSYMBOL
	if(species == "Mm") map <- org.Mm.egSYMBOL
	par(mfrow=c(1,1))
	map<-as.list(map)
	cat("-Network analysis complete \n")
	#This return a list with a graphNEL object (the module), and a mapping from the nodes of this module to their symbol identifiers	
	return(list(subnw=subnw,labels=map[nodesModule]))
}

