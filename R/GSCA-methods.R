
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
##map gene set ids to terms, and insert terms to result data frames
##and merge hyperGeo results and GSEA results together to create a integrated data frame
##and search gene sets that are both significant in hyperGeo and GSEA analyses
setMethod(
		"appendGSTerms",
		"GSCA",
		function(object, keggGSCs=NULL, goGSCs=NULL) {
			
			if(length(object@result)==0)
				stop("No results generated!\n")
			gsc.names<-names(object@listOfGeneSetCollections)	
			if(!is.null(keggGSCs)) {
				paraCheck(name = "keggGSCs", para = keggGSCs)
				if(!all(keggGSCs %in% gsc.names))
					stop("Wrong gene set collection names specified in 'keggGSCs'!\n")
			}
			if(!is.null(goGSCs)) {
				paraCheck(name = "goGSCs", para = goGSCs)
				if(!all(goGSCs %in% gsc.names))
					stop("Wrong gene set collection names specified in 'goGSCs'!\n")
			}
			appendGOTerm<-function(df) {
				data.frame(Gene.Set.Term=Term(rownames(df)),df, stringsAsFactors=FALSE)
			}
			appendKEGGTerm<-function(df) {
				dfRownames<-rownames(df)
				gsKEGG<-sapply(dfRownames, function(c) substr(c,4,nchar(c)))
				newdf<-data.frame(Gene.Set.Term=unlist(mget(gsKEGG, env=KEGGPATHID2NAME)),df, stringsAsFactors=FALSE)
				rownames(newdf)<-rownames(df)
				return(newdf)
			}	

			##add gene set terms if possible
			for(rs in 1:length(object@result)) {
				if(names(object@result)[rs]%in%c("HyperGeo.results","GSEA.results","Sig.pvals.in.both","Sig.adj.pvals.in.both")) {
					sapply(names(object@result[[rs]]), function(gsc) {
						if(gsc %in% names(object@listOfGeneSetCollections)) {
							if(nrow(object@result[[rs]][[gsc]])>=1) {
								if(gsc %in%keggGSCs)
									object@result[[rs]][[gsc]]<<-appendKEGGTerm(object@result[[rs]][[gsc]])
								else if(gsc %in%goGSCs)
									object@result[[rs]][[gsc]]<<-appendGOTerm(object@result[[rs]][[gsc]])
								else 
									object@result[[rs]][[gsc]]<<-data.frame(Gene.Set.Term="--", object@result[[rs]][[gsc]], stringsAsFactors=FALSE)	
							} else {
								object@result[[rs]][[gsc]]<<-cbind(Gene.Set.Term=NULL,object@result[[rs]][[gsc]])	
							}
						}
					})
				}
			}
			object
		}
)

##show summary information on screen
setMethod(
		"show",
		"GSCA",
		function(object) {
			cat("A GSCA (Gene Set Collection Analysis) object:\n")
			summarize(object, what=c("GSC","GeneList","Hits","Para"))
		}
)
##print summary information on screen
setMethod(
		"summarize",
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
			paraCheck(name="resultName",para=resultName)
			if(!(resultName %in% names(object@result)))
				stop("Please input 'HyperGeo.results' or 'GSEA.results'!\n")
			gsc.names<-names(object@result[[resultName]])
			if(!all(gscs %in% gsc.names))
				stop("Wrong Gene Set Collection name(s) in 'gscs'! \n")
			if(!is.null(ntop))
				paraCheck(name="ntop",para=ntop)
			paraCheck(name="allSig",para=allSig)
			if((is.null(ntop) && !allSig)||(!is.null(ntop) && allSig))
				stop("Either specify 'ntop' or set 'allSig' to be TRUE!\n")	
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
##plot GSEA for GSCA
setMethod(
		"plotEnrichMap",
		"GSCA",
		function(object, gscs, ntop=NULL, allSig=TRUE, gsNameType="id", displayEdgeLabel=TRUE, 
				layout="layout.fruchterman.reingold", filepath=".", filename="test.png",output="png", ...) {
			##check arguments
			paraCheck(name="filepath", para=filepath)
			paraCheck(name="output", para=output)
			paraCheck(name="filename", para=filename)
			if(output == "pdf" ) 
				pdf(file.path(filepath, filename), ...=...)
			if(output == "png" ) 
				png(file.path(filepath, filename), ...=...)
			viewEnrichMap(object, gscs, ntop, allSig, gsNameType, displayEdgeLabel, layout)
			dev.off()
		}
)
setMethod(
		"viewEnrichMap",
		"GSCA",
		function(object, gscs, ntop=NULL, allSig=TRUE, gsNameType="id", displayEdgeLabel=TRUE, layout="layout.fruchterman.reingold") {
			##check arguments			
			if(missing(gscs))
				stop("Please specify the name(s) of Gene Set Collections in 'gscs'! \n")
			paraCheck(name="gscs.names",para=gscs)
			resultName<-"GSEA.results"
			if(!(resultName %in% names(object@result)))
				stop("No GSEA results found in object!\n")
			gsc.names<-names(object@result[[resultName]])
			if(!all(gscs %in% gsc.names))
				stop("Wrong Gene Set Collection name(s) in 'gscs'! \n")
			if(!is.null(ntop))
				paraCheck(name="ntop",para=ntop)
			paraCheck(name="allSig",para=allSig)
			if((is.null(ntop) && !allSig)||(!is.null(ntop) && allSig))
				stop("Either specify 'ntop' or set 'allSig' to be TRUE!\n")	
			paraCheck(name="gsNameType", gsNameType)
			paraCheck(name="displayEdgeLabel", displayEdgeLabel)
			paraCheck("layout", layout)
			##get top gene sets
			topGS<-getTopGeneSets(object, resultName, gscs, ntop, allSig)
			if(sum(unlist(lapply(topGS, length)))==0)
				stop("No significant gene sets found!\n")
			gsInUni<-list()
			uniIDs<-names(object@geneList)
			tempList<-list()
			junk<-sapply(1:length(topGS), function(i) {
				if(length(topGS[[i]])>0) {
					gsc.name<-names(topGS)[i]
					##compute overlapped genes between gene sets and universe
					gsInUni[[i]]<<-list()
					gsInUni[[i]]<<-sapply(topGS[[i]], function(j) intersect(object@listOfGeneSetCollections[[gsc.name]][[j]], uniIDs),simplify=FALSE)
					names(gsInUni)[i]<<-gsc.name
					##if(!is.null(keggGSCs) && gsc.name %in% keggGSCs) {
					##	thisGS<-sapply(topGS[[i]], function(c) substr(c,4,nchar(c)))
					##	gsTerms<-unlist(mget(thisGS, env=KEGGPATHID2NAME))
					##} else if(!is.null(goGSCs) && gsc.name %in% goGSCs) {
					##	gsTerms<-Term(topGS[[i]])
					##}
					tempList[[i]]<<-data.frame(gsID=topGS[[i]], gscID=gsc.name, object@result[[resultName]][[gsc.name]][topGS[[i]],,drop=FALSE])
					names(tempList)[i]<<-gsc.name
				}
			})		
			##collapse to a data frame
			tempdf<-do.call("rbind", lapply(tempList, data.frame, stringsAsFactors = FALSE))
			if(gsNameType=="term" && !("Gene.Set.Term" %in% colnames(tempdf)))
				stop("No gene set terms found in results!\n Please use the method 'appendGSTerms' or add a column named 'Gene.Set.Term' to the results!\n")
			##function to compute overlapped genes
			map.mat<-matrix(0,nrow(tempdf),nrow(tempdf))
			diag(map.mat)<-1
			map.diag<-sapply(1:nrow(tempdf), function(i) length(gsInUni[[as.character(tempdf[i,"gscID"])]][[as.character(tempdf[i,"gsID"])]]))
			if(nrow(tempdf)>=2) {
				sapply(1:(nrow(tempdf)-1), function(i) {	
							map.mat[i, (i+1):nrow(tempdf)]<<-sapply((i+1):nrow(tempdf), function(j) {
										length(intersect(gsInUni[[as.character(tempdf[i,"gscID"])]][[as.character(tempdf[i,"gsID"])]],
																gsInUni[[as.character(tempdf[j,"gscID"])]][[as.character(tempdf[j,"gsID"])]]))/
												length(union(gsInUni[[as.character(tempdf[i,"gscID"])]][[as.character(tempdf[i,"gsID"])]],
																gsInUni[[as.character(tempdf[j,"gscID"])]][[as.character(tempdf[j,"gsID"])]]))
									})	
							map.mat[(i+1):nrow(tempdf),i]<<-map.mat[i, (i+1):nrow(tempdf)]
				})		
				rownames(map.mat)<-rownames(tempdf)
				colnames(map.mat)<-rownames(tempdf)	
				##generate igraph from adjacency matrix
				###'Node name' controlled by the rownames of tempList
				g<-graph.adjacency(adjmatrix=map.mat, mode="undirected", weighted=TRUE, diag=TRUE)
				g<-simplify(g, remove.loops = TRUE)
			} else if(nrow(tempdf)==1) {
				diag(map.mat)<-0
				rownames(map.mat)<-rownames(tempdf)
				colnames(map.mat)<-rownames(tempdf)	
				##generate igraph from adjacency matrix
				###'Node name' controlled by the rownames of tempList
				g<-graph.adjacency(adjmatrix=map.mat, mode="undirected", weighted=NULL, diag=FALSE)
			}
			if(length(V(g))>=2) {
				###'Node size' controlled by the 'size of gene set'
				v.max.size<-18
				v.min.size<-4
				if(max(map.diag)!=min(map.diag))
					V(g)$size<-v.min.size+(v.max.size-v.min.size)*(map.diag-min(map.diag))/(max(map.diag)-min(map.diag))
				else
					V(g)$size<-6
			} else if(length(V(g))==1) {
				V(g)$size<-4
			}
			p.vec<-tempdf[,"Adjusted.Pvalue"]
			p.cutoff.vec<-c(0, 10^c(-3, -2.5), 0.01, 10^(-1.5), 0.05, 10^(-c(1.0, 0.5, 0)))
			
			posids<-which(tempdf[,"Observed.score"]>=0)
			negids<-which(tempdf[,"Observed.score"]<=0)
			
			redCols<-colorRampPalette(colors = c("red", "white"))
			redVec<-redCols(length(p.cutoff.vec))

			blueCols<-colorRampPalette(colors = c("blue", "white"))
			blueVec<-blueCols(length(p.cutoff.vec))
			V(g)$color<-""
			if(length(posids)>0)
				V(g)$color[posids]<-redVec[as.integer(cut(x=p.vec[posids],breaks=c(-1,p.cutoff.vec), labels=1:(length(p.cutoff.vec))))]
			if(length(negids)>0)
				V(g)$color[negids]<-blueVec[as.integer(cut(x=p.vec[negids],breaks=c(-1,p.cutoff.vec), labels=1:(length(p.cutoff.vec))))]
			##labels attributes
			graphLabelWrapper<-function(x, width=32) {paste(strwrap(x,width=width),collapse="\n")}
			if(gsNameType=="id") {
				V(g)$label<-as.character(tempdf[,"gsID"])
			} else if(gsNameType=="term") {
				templabels<-as.character(tempdf[,"Gene.Set.Term"])
				V(g)$label<-sapply(templabels, graphLabelWrapper)
			} 
				
			V(g)$label.dist<-0.4
			V(g)$label.cex<-0.75
			V(g)$label.font<-3
			V(g)$label.color<-"black"
			###'Node color' controlled by the 'adjusted pvalue'
			if(length(V(g))>=2)
				E(g)$color<-grey(0.7)
			if(displayEdgeLabel) {
				edgeWeights<-round(E(g)$weight*100)
				edgeWeights[edgeWeights==0]<-""
				E(g)$label<-edgeWeights
			}
			###'Edge thickness' controlled by the 'size of overlapped genes' between two gene sets
			edge.max.w<-14
			edge.min.w<-1
			if(length(V(g))>=2) {
				edgeWeightVec<-round(edge.min.w+(edge.max.w-edge.min.w)*(E(g)$weight))	
				E(g)$width<-edgeWeightVec
			} 
			##plot graph
			plot(g,layout=eval(parse(text=layout)))
			##title
			title(main=paste("Enrichment Map of gene set collection(s)", lapply(list(gscs), paste, collapse=",")[[1]], sep="--"))
			##p-value color legend
			colVec<-c(redVec[1:(length(redVec)-1)],rev(blueVec))
			points(
					x = rep(-1.2, length(colVec)), 
					y = seq(0.5, (0.5-(0.05*length(colVec))), 
							length.out = length(colVec)), 
					pch = 15, col = colVec
			)
			p.cutoff.labels<-rep("",length(colVec))
			p.cutoff.labels[c(1,4,6,9,11,14,17)]<-c(0,0.01,0.05,1,0.05,0.01,0)
			text(
					x = rep(-1.3, length(colVec)), 
					y = seq(0.5, (0.5-(0.05*length(colVec))), 
							length.out = length(colVec)),
					labels = p.cutoff.labels, 
					cex = 0.8,
					adj=1
			)
			text(
					x = -1.25,
					y = 0.7,
					labels = "Adjusted\np-values",
					cex=0.8,
					adj= 0.5,
					font=2
			)
			return(g)
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

















