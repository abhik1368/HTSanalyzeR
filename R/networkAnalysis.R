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
	species=c("Drosophila_melanogaster"),
	initialIDs="FlybaseCG",
	fdr=0.001,
	genetic=FALSE,
	phenotypeVector=NA,
	biogridObject=NA,
	order=2,
	link="http://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-3.0.64/BIOGRID-ORGANISM-3.0.64.tab2.zip"){
#Do the statistical tests: this function will test:
#-that the cellHTSobject is a cellHTS object
#that the cellHTS2 cellHTSobject is in the right format (configured, so that we know which rows are samples and which are controls)
#that the annotationColumn has been specified as a single character string
#that the annotationColumn matches a column in the fData(cellHTSobject) dataframe
#that the 'controls' parameter has been specified as a single character string
#that the  'controls' parameter matches a status in the 'controlStatus' column of the cellHTS cellHTSobject
#that the 'alternative' parameter has been correctly set 	
	test.stats<-cellHTS2OutputStatTests(cellHTSobject=cellHTSobject,alternative=alternative,tests=tests)
#Check that the species is correctly specified (corresponds to a single character string 
#that can be matched to an input for the AnnotationConvertor functions)
	if(is.character(species) != TRUE | length(species) != 1) stop("The species argument does not have the right format")
	if(is.na(match(species,c("Drosophila_melanogaster","Homo_sapiens","Rattus_norvegicus","Mus_musculus","Caenorhabditis_elegans")))) {
		stop("The species argument does not match any of the names recognized by this function, please provide one of the following character strings: 'Drosophila_melanogaster','Homo_sapiens','Rattus_norvegicus','Mus_musculus','Caenorhabditis_elegans'")
		}	
#The biogridObject contains Entrez.gene IDs, so if this is not the initial type of identifiers,
#they need to be mapped to Entrez.gene IDs		
	if(initialIDs != "Entrez.gene"){
#Check that the initialIDs argument is correctly specified
		if(length(intersect(initialIDs, c("Ensembl.transcript", "Ensembl.prot", "Ensembl.gene", "Entrez.gene", "RefSeq", "Symbol", "GenBank", "Flybase", "FlybaseCG", "FlybaseProt"))) != 1) {
			stop('The initialIDs argument is not correctly specified, please provide one of the following character strings: "Ensembl.transcript", "Ensembl.prot", "Ensembl.gene", "Entrez.gene", "RefSeq", "Symbol", "GenBank", "Flybase", "FlybaseCG", "FlybaseProt"')
			}
#map the identifiers to Entrez.gene identifiers			
		if(species == "Drosophila_melanogaster") test.stats.entrez<-drosoAnnotationConvertor(geneList=test.stats,initialIDs=initialIDs,finalIDs="Entrez.gene")
		if(species == "Caenorhabditis_elegans") test.stats.entrez<-celAnnotationConvertor(geneList=test.stats,initialIDs=initialIDs,finalIDs="Entrez.gene")
		if(species == "Homo_sapiens") test.stats.entrez<-mammalAnnotationConvertor(geneList=test.stats,species=species,initialIDs=initialIDs,finalIDs="Entrez.gene")
		if(species == "Rattus_norvegicus") test.stats.entrez<-mammalAnnotationConvertor(geneList=test.stats,species=species,initialIDs=initialIDs,finalIDs="Entrez.gene")
		if(species == "Mus_musculus") test.stats.entrez<-mammalAnnotationConvertor(geneList=test.stats,species=species,initialIDs=initialIDs,finalIDs="Entrez.gene")
		}
	subnw<-enrichedSubNw(
	pvaluesMatrix=test.stats.entrez,
	columns=columns,
	species=species,
	fdr=fdr,genetic=genetic,
	biogridObject=biogridObject,order=order)
	nodesModule<-nodes(subnw)
#To represent the network in a more convenient format, the symbol identifiers will be mapped and given to the user (more readable than Entrez.gene IDs)	
	if(species == "Drosophila_melanogaster") map <- org.Dm.egSYMBOL
	if(species == "Caenorhabditis_elegans") map <- org.Ce.egSYMBOL
	if(species == "Homo_sapiens") map <- org.Hs.egSYMBOL
	if(species == "Rattus_norvegicus") map <- org.Rn.egSYMBOL
	if(species == "Mus_musculus") map <- org.Mm.egSYMBOL
	par(mfrow=c(1,1))
	map<-as.list(map)
#If no phenotype vector is specified, then we can just plot the module	
	if(length(phenotypeVector) == 1){
		png("EnrichedSubNw.png",width = 900, height = 900)
		plotModule(subnw,labels=map[nodesModule])
		dev.off()
		}else{
#If a phenotype vector is specified, it will be used to color the nodes
#check that the phenotypeVector has the right format (single numerical vector)
			if(class(phenotypeVector) != "numeric" && class(phenotypeVector) != "integer") stop("The phenotypeVector should be a single vector of numerical phenotypes")
#Check that the phenotypeVector is a named vector
			if(is.null(dim(phenotypeVector)) != TRUE | is.null(names(phenotypeVector))) stop("Please provide a phenotypeVector as a single named vector")	
#this vector holds the phenotype for the nodes of the sub-network	
			diff.expr<-phenotypeVector[nodes(subnw)]
			names(diff.expr)<-nodes(subnw)
#this vector contains the information of wether a node has an associated phenotype (1) or not (-1)
#this information will be used to give a different shape to the nodes of the network
			present<-rep(1,length(nodes(subnw)))
			present[which(is.na(diff.expr))]<--1
#this replaces all phenotypes of non-phenotyped nodes by a zero		
			diff.expr[which(is.na(diff.expr))]<-0
			names(present)<-nodes(subnw)
#Plot the module	
			if(length(nodes(subnw)) == 1){
				png("EnrichedSubNw.png",width = 900, height = 900)
				plotModule(subnw,labels=map[nodesModule],scores=present,diff.expr=diff.expr)
				dev.off()
				}else{
					Tcolors<-diff.expr
					Tcolors2<-diff.expr
					if(max(abs(Tcolors))<5)  {Tcolors <- Tcolors*5}
					# set red colors
					if(any(Tcolors>0)){
						maxRed <- max(ceiling(abs(Tcolors[which(Tcolors>0)])))
						redCols <- colorRampPalette(colors=c("white", "red"))
						redVec <- redCols(maxRed)
						Tcolors2[which(Tcolors>0)] <- redVec[ceiling(abs(Tcolors[which(Tcolors>0)]))]
						}
					#set the greens
					if(any(Tcolors<0)){
						maxGreen <- max(ceiling(abs(Tcolors[which(Tcolors<0)])))
						greenCols <- colorRampPalette(colors=c("white", "green"))
						greenVec <- greenCols(maxGreen)
						Tcolors2[which(Tcolors<0)] <- greenVec[ceiling(abs(Tcolors[which(Tcolors<0)]))]
						}	  	
					colScale<-unique(Tcolors2)
					colboundary<-rep(0,length(colScale))
					for(i in 1:length(colScale)){
						values<-diff.expr[which(Tcolors2 == colScale[i])]
						colboundary[i]<-values[which(abs(values) == max(abs(values)))[1]]
						}
					colMatrix<-cbind(colboundary[order(colboundary)],colScale[order(colboundary)])
					png("EnrichedSubNw.png", width = 900, height = 900)
					plotModule(subnw,labels=map[nodesModule],scores=present,diff.expr=diff.expr)
					points(x=rep(-1.2,length(unique(Tcolors2))),y=seq(1.2,(1.2-(0.05*length(colMatrix[,2]))),length.out=length(colMatrix[,2])),pch=15,col=colMatrix[,2])
					text(x=rep(-1.1,length(unique(Tcolors2))),y=seq(1.2,(1.2-(0.05*length(colMatrix[,2]))),length.out=length(colMatrix[,2])),labels=signif(as.numeric(colMatrix[,1]),digits=2),cex=0.8)
					dev.off()
					}
			}
#This plots the module on the screen, in 2D and in 3D (the latter is a dynamic view)	
	plotModule(subnw,labels=map[nodesModule])
	if(length(nodes(subnw)) > 1) {plot3dModule(subnw,labels=map[nodesModule])}
#This saves the module in a sif and in a tab delimited format	
	if(length(nodes(subnw)) > 1) {
		saveNetwork(subnw,name="EnrichedSubNw",file="EnrichedSubNw.sif",type="sif")
		saveNetwork(subnw,name="EnrichedSubNw",file="EnrichedSubNw.txt",type="table")
		}
#This saves the symbol identifiers of the nodes of the network, to be used as attributes in cytoscape for example	
	write.table(x=cbind(nodes(subnw),map[nodesModule]),file="EnrichedSubNwnodeattr.txt",row.names=FALSE,col.names=FALSE,sep="\t")	
#This return a list with a graphNEL object (the module), and a mapping from the nodes of this module to their symbol identifiers	
	return(list(subnw=subnw,labels=map[nodesModule]))
	}

