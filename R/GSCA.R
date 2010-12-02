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
##initialization method
setMethod("initialize",
		"GSCA",
		function(.Object, listOfGeneSetCollections,geneList,hits) {
			#######################
			##check input arguments
			#######################
			paraCheck(name="gscs",listOfGeneSetCollections)
			paraCheck(name="genelist",geneList)
			paraCheck(name="hits",hits)
			#[para<-paraCheck(name="gsca.para",para)
			#######################
			##	initialization
			#######################
			.Object@listOfGeneSetCollections<-listOfGeneSetCollections
			.Object@geneList<-geneList
			.Object@hits<-hits
			##[.Object@para<-para]
			##check result summary and preprocessed to make sure they are not specified by users
			.Object@result<-list()
			.Object@preprocessed<-FALSE
			
			##summary info--framework
			###gene set collections
			sum.info.gsc<-matrix(,length(listOfGeneSetCollections),2)
			rownames(sum.info.gsc)<-names(listOfGeneSetCollections)
			colnames(sum.info.gsc)<-c("input",paste("above min size",sep=""))
			###gene list
			sum.info.gl<-matrix(,1,4)
			colnames(sum.info.gl)<-c("input","valid","duplicate removed","converted to entrez")
			rownames(sum.info.gl)<-"Gene List"
			###hits
			sum.info.hits<-matrix(,1,2)
			colnames(sum.info.hits)<-c("input","preprocessed")
			rownames(sum.info.hits)<-"Hits"
			###parameters
			sum.info.para<-list()
			sum.info.para$hypergeo<-matrix(,1,3)
			colnames(sum.info.para$hypergeo)<-c("minGeneSetSize","pValueCutoff","pAdjustMethod")
			rownames(sum.info.para$hypergeo)<-"HyperGeo Test"
			
			sum.info.para$gsea<-matrix(,1,5)
			colnames(sum.info.para$gsea)<-c("minGeneSetSize","pValueCutoff","pAdjustMethod","nPermutations","exponent")
			rownames(sum.info.para$gsea)<-"GSEA"			
			###results
			sum.info.results<-matrix(,3,length(listOfGeneSetCollections))
			colnames(sum.info.results)<-names(listOfGeneSetCollections)
			rownames(sum.info.results)<-c("HyperGeo","GSEA","Both")
			
			##summary info--initialize
			sum.info.gsc[,"input"]<-unlist(lapply(listOfGeneSetCollections,length))
			sum.info.gl[,"input"]<-length(geneList)
			sum.info.hits[,"input"]<-length(hits)
			##[sum.info.para$hypergeo[1,]<-c(para$minGeneSetSize,para$pValueCutoff,para$pAdjustMethod)
			##[sum.info.para$gsea[1,]<-c(para$minGeneSetSize,para$pValueCutoff,para$pAdjustMethod,para$nPermutations,para$exponent)
			##[sum.info.para$hypergeo<-data.frame(sum.info.para$hypergeo)
			##[sum.info.para$gsea<-data.frame(sum.info.para$gsea)
			
			
			.Object@summary<-list(gsc=sum.info.gsc,gl=sum.info.gl,hits=sum.info.hits,para=sum.info.para,results=sum.info.results)			
			.Object
		}
)
##pre-processing
setMethod(
		"preprocess",
		"GSCA",
		function(object, species="Dm", initialIDs="FlybaseCG", keepMultipleMappings=TRUE, duplicateRemoverMethod="max", orderAbsValue=FALSE, verbose=TRUE) {
			#######################
			##check input arguments
			#######################
			paraCheck(name="species",para=species)
			paraCheck(name="initialIDs",para=initialIDs)
			paraCheck(name="keepMultipleMappings",para=keepMultipleMappings)
			paraCheck(name="duplicateRemoverMethod",para=duplicateRemoverMethod)
			paraCheck(name="orderAbsValue",para=orderAbsValue)
			paraCheck(name="verbose",para=verbose)
			#######################
			##	 preprocessing
			#######################
			cat("-Preprocessing for input gene list and hit list ...\n")
			
			genelist<-object@geneList
			##remove NA in geneList
			if(verbose) cat("--Removing genes without values in geneList ...\n")
			genelist<-genelist[which((!is.na(genelist)) & (names(genelist)!="") & (!is.na(names(genelist))))]
			object@summary$gl[,"valid"]<-length(genelist)				#genes with valid values
			if(length(genelist)==0)
				stop("Input 'geneList' contains no useful data!\n")
			##duplicate remover
			if(verbose) cat("--Removing duplicated genes ...\n")
			genelist<-duplicateRemover(geneList=genelist,method=duplicateRemoverMethod)
			object@summary$gl[,"duplicate removed"]<-length(genelist)	#genes after removing duplicates 

			hits<-object@hits[object@hits!="" & !is.na(object@hits)]
			if(length(hits)==0)
				stop("Input 'hits' contains no useful data!\n")
			hits<-unique(hits)
			match.ind<-match(hits,names(genelist))
			
			hits.vec<-genelist[match.ind[!is.na(match.ind)]]
			if(length(hits.vec)==0) 
				stop("Hits and geneList have no overlaps!\n")
				
			##annotation convertor
			if(initialIDs!="Entrez.gene") {
				if(verbose) cat("--Converting annotations ...\n")
				genelist<-annotationConvertor(
						geneList=genelist,
						species=species,
						initialIDs=initialIDs,
						finalIDs="Entrez.gene",
						keepMultipleMappings=keepMultipleMappings,
						verbose=verbose
				)
				hits.vec<-annotationConvertor(
						geneList=hits.vec,
						species=species,
						initialIDs=initialIDs,
						finalIDs="Entrez.gene",
						keepMultipleMappings=keepMultipleMappings,
						verbose=verbose
				)
			}
			
			object@summary$gl[,"converted to entrez"]<-length(genelist)	#genes after annotation conversion
			
			if(verbose) cat("--Ordering Gene List decreasingly ...\n")
			if(!orderAbsValue)
				genelist<-genelist[order(genelist,decreasing=TRUE)]
			else
				genelist<-abs(genelist)[order(abs(genelist),decreasing=TRUE)]	
			hits<-names(hits.vec)
			object@summary$hits[,"preprocessed"]<-length(hits)			#hits after preprocessed
			##update genelist and hits, and return object
			object@geneList<-genelist
			object@hits<-hits
			object@preprocessed<-TRUE
			
			cat("-Preprocessing complete!\n\n")
			
			object
		}
)
##Gene set enrichment analyses
setMethod(
		"analyze",
		"GSCA",
		function(object, para=list(pValueCutoff = 0.05,pAdjustMethod = "BH", nPermutations = 1000, minGeneSetSize = 15,exponent = 1), verbose=TRUE) {
			paraCheck(name="verbose",para=verbose)
			object@para<-paraCheck(name="gsca.para",para=para)
			object@summary$para$hypergeo[1,]<-c(object@para$minGeneSetSize,object@para$pValueCutoff,object@para$pAdjustMethod)
			object@summary$para$gsea[1,]<-c(object@para$minGeneSetSize,object@para$pValueCutoff,object@para$pAdjustMethod,object@para$nPermutations,object@para$exponent)
			##if(!is.data.frame(object@summary$para$hypergeo))
			##	object@summary$para$hypergeo<-data.frame(object@summary$para$hypergeo)
			##if(!is.data.frame(object@summary$para$gsea))
			##	object@summary$para$gsea<-data.frame(object@summary$para$gsea)
			#######################
			##	  do analysis
			#######################
			object@result<-
					analyzeGeneSetCollections(
					listOfGeneSetCollections=object@listOfGeneSetCollections,
					geneList=object@geneList,
					hits=object@hits,
					pAdjustMethod=object@para$pAdjustMethod,
					pValueCutoff=object@para$pValueCutoff,
					nPermutations=object@para$nPermutations, 
					minGeneSetSize=object@para$minGeneSetSize,
					exponent=object@para$exponent,
					verbose=verbose
			)
			##update summary information
			cols<-colnames(object@summary$results)
			object@summary$gsc[,2]<-unlist(lapply(object@result$HyperGeo.results, nrow))[cols]

			object@summary$results["HyperGeo",]<-unlist(lapply(object@result$HyperGeo.results,function(df) {sum(df[,"Adjusted.Pvalue"]<object@para$pValueCutoff)}))[cols]
			object@summary$results["GSEA",]<-unlist(lapply(object@result$GSEA.results,function(df) {sum(df[,"Adjusted.Pvalue"]<object@para$pValueCutoff)}))[cols]
			object@summary$results["Both",]<-unlist(lapply(object@result$Sig.adj.pvals.in.both,nrow))[cols]
			
			object
		}
)
##show summary information on screen
setMethod(
		"show",
		"GSCA",
		function(object) {
			cat("A GSCA (Gene Set Collection Analysis) object:\n")
			summary(object, what=c("GSC","GeneList","Hits","Para"))
		}
)
##print summary information on screen
setMethod(
		"summary",
		"GSCA",
		function(object, what="ALL") {
			paraCheck("what.gsca",what)
			##what can be "GSC" (gene set collection), "GeneList", "Hits", "Para", "Result"
			if(any(c("ALL","GSC") %in% what)) {
				cat("\n")
				cat("-No of genes in Gene set collections: \n")
				print(object@summary$gsc,quote=FALSE)
				cat("\n")
			}
			if(any(c("ALL","GeneList") %in% what)) {
				cat("\n")
				cat("-No of genes in Gene List: \n")
				print(object@summary$gl,quote=FALSE)
				cat("\n")
			}
			if(any(c("ALL","Hits") %in% what)) {
				cat("\n")
				cat("-No of hits: \n")
				print(object@summary$hits,quote=FALSE)
				cat("\n")
			}
			if(any(c("ALL","Para") %in% what)) {
				cat("\n")
				cat("-Parameters for analysis: \n")
				print(object@summary$para$hypergeo,quote=FALSE)
				cat("\n")
				print(object@summary$para$gsea,quote=FALSE)
				cat("\n")
			}
			if(any(c("ALL","Result") %in% what)) {
				cat("\n")
				cat("-Significant gene sets (adjusted p-value<",object@para$pValueCutoff,"): \n")
				print(object@summary$results,quote=FALSE)
				cat("\n")
			}	
		}
)
##select top significant gene sets 
setMethod(
		"getTopGeneSets",
		"GSCA",
		function(object, resultName, gscs, ntop=NULL, allSig=FALSE) {
			##check arguments			
			if(missing(gscs))
				stop("Please specify the name(s) of Gene Set Collections in 'gscs'! \n")
			paraCheck(name="gscs.names",para=gscs)
			gsc.names<-names(object@result[[resultName]])
			if(!all(gscs %in% gsc.names))
				stop("Wrong Gene Set Collection name(s) in 'gscs'! \n")
			if(!is.null(ntop))
				paraCheck(name="ntop",para=ntop)
			paraCheck(name="allSig",para=allSig)
			if((is.null(ntop) && !allSig)||(!is.null(ntop) && allSig))
				stop("Either specify 'ntop' or set 'allSig' to be TRUE!\n")
			paraCheck(name="resultName",para=resultName)
			filenames<-list()
			for(gsc in gscs) {
				all.gs.names <- rownames(object@result[[resultName]][[gsc]])
				##if(nrow(object@result[[resultName]][[gsc]])==0)
				##	stop("No gene sets in this gene set collection!\n")
				##find all gene sets that are going to be plot
				gs.names<-NULL
				if(allSig) {
					gs.names<-all.gs.names[object@result[[resultName]][[gsc]][,"Adjusted.Pvalue"]<object@para$pValueCutoff]
				##	if(is.null(gs.names))
				##		stop("No significant gene sets found in specified gene set collection!\n")
				} else {
					if(ntop>nrow(object@result[[resultName]][[gsc]]))
						stop("'ntop' is larger than the number of gene sets in specified gene set collection!\n")
					gs.names<-all.gs.names[1:ntop]
				} 		
				filenames[[gsc]]<-unlist(lapply(list(gs.names),gsub, pattern="/",replacement="_"))
				names(filenames[[gsc]])<-gs.names
			}
			return(filenames)
		}
)
##write observed hits in gene sets for HyperGeo
setMethod(
		"writeHits",
		"GSCA",
		function(object, gscs, species=NULL, ntop=NULL, allSig=FALSE, filepath=".") {
			##check arguments
			paraCheck(name="filepath",para=filepath)
			if(!is.null(species)) {
				paraCheck(name="species",para=species)
				anno.db.species<-paste("org",species,"eg","db",sep=".")
				if(! (paste("package:",anno.db.species, sep="") %in% search()))
					library(anno.db.species,character.only=TRUE)
				mapID <- as.list(get(paste("org", species, "egSYMBOL", sep=".")))
			} else 
				mapID<-NULL
			filenames<-getTopGeneSets(object, "HyperGeo.results", gscs, ntop, allSig)
			for(gsc in gscs) {
				##write for all gs.names	
				gs.names<-filenames[[gsc]]
				if(!is.null(gs.names)) {
					for(gs.name in gs.names){
						makeOverlapTable(geneSet=object@listOfGeneSetCollections[[gsc]][[gs.name]], hits=object@hits,
								mapID=mapID, filepath=filepath, filename=gs.name)
					}
				}
			}
		}		
)
##plot GSEA for GSCA
setMethod(
		"plotGSEA",
		"GSCA",
		function(object, gscs, ntop=NULL, allSig=FALSE, filepath=".", output="png", ...) {
			##check arguments
			paraCheck(name="filepath", para=filepath)
			paraCheck(name="output", para=output)
			paraCheck("allSig", allSig)
			if(!is.null(ntop))
				paraCheck("ntop", ntop)
			paraCheck("gscs.names", gscs)
			filenames<-getTopGeneSets(object, "GSEA.results", gscs, ntop, allSig)				
			for(gsc in gscs) {
				##plot for all gs.names	
				gs.names<-filenames[[gsc]]
				if(!is.null(gs.names)) {
					for(gs.name in gs.names){
						makeGSEAplots(geneList=object@geneList, geneSet=object@listOfGeneSetCollections[[gsc]][[gs.name]],
								exponent=object@para$exponent, filepath=filepath,
								filename=gsub("/", "_", gs.name), output=output, ...=...)
					}	
				}
			}
		}
)
setMethod(
		"viewGSEA",
		"GSCA",
		function(object, gscName, gsName) {
			##check argument
			paraCheck("gs.single", gsName)
			paraCheck("gsc.name", gscName)
			if(!("GSEA.results" %in% names(object@result)))
				stop("GSEA not performed!\n")
			gs.all<-lapply(object@result[["GSEA.results"]][1:length(object@listOfGeneSetCollections)], rownames)
			if(length(unlist(gs.all))==0)
				stop("No gene sets in GSEA results!\n")
			if(!(gsName %in% unlist(gs.all)))
				stop("'gs' is not a gene set that passes the 'minGeneSetSize'! \n")
			if(!(gscName %in% names(object@listOfGeneSetCollections)))
				stop("'gsc' is not a gene set collection in 'listOfGeneSetCollections'!\n")
			##gsc.name<-NULL
			##sapply(names(object@listOfGeneSetCollections), function(gsc) {
			##	gs.id<-which(gs.all[[gsc]]==gs)
			##	gsc.name<<-ifelse(length(gs.id)>0, gsc, gsc.name)
			##})
			test <- gseaScores(geneList = object@geneList, geneSet = object@listOfGeneSetCollections[[gscName]][[gsName]],
					exponent = object@para$exponent, mode = "graph")
			gseaPlots(runningScore = test[['runningScore']],
					enrichmentScore = test[['enrichmentScore']],
					positions = test[['positions']], geneList = object@geneList)
		}
)
##generate html report for GSCA
setMethod(
		"report",
		"GSCA",
		function(object, experimentName="Unknown", species=NULL, ntop=NULL, allSig=FALSE, keggGSCs=NULL, goGSCs=NULL, reportDir="HTSanalyzerReport") {
			##call writeReportHTSA
			writeReportHTSA(gsca=object, experimentName=experimentName, species=species, ntop=ntop, allSig=allSig, 
					keggGSCs=keggGSCs, goGSCs=goGSCs, reportDir=reportDir)
		}
)

















