#This function writes an html report following a complete analysis of a dataset 
#with the two wrapper functions geneSetAnalyze and networkAnalysis. 
#Most of the parameters of this function can simply be copied and paste 
#from the parameters of the two above functions.
writeReportHTSA <-
function(experimentName,enrichmentAnalysis,cutoffHits,hits,listOfGeneSetCollections,geneList,geneListName,p.adj.method="BH",nPermutations=1000,min.gene.set.size=15,exponent=1,nwAnalysisOutput,nwAnalysisGraphFile="EnrichedSubNw.png",controls="neg",alternative="two.sided",tests=c("T-test"),columns=c("t.test.pvalues.two.samples","t.test.pvalues.one.sample"),species=c("Drosophila_melanogaster"),fdr=0.001,genetic=FALSE,biogridObject=NA,nGseaPlots=10,whichSetIsKEGGIds,whichSetIsGOIds){
#Get the indices of the gene set collections for which there will be links to the relevant databases
	if(length(whichSetIsKEGGIds) == 1) {
		if(is.na(whichSetIsKEGGIds)){
		keggGS<-0
		}else{
			keggGS<-whichSetIsKEGGIds
			}
		}else{	
			keggGS<-whichSetIsKEGGIds
			}
	if(length(whichSetIsGOIds) == 1) {
		if(is.na(whichSetIsGOIds)){
		GOGS<-0
		}else{
			GOGS<-whichSetIsGOIds
			}
		}else{	
			GOGS<-whichSetIsGOIds
			}	
#Check that listOfGeneSetCollections is a named list
	if(is.list(listOfGeneSetCollections) == FALSE | is.null(names(listOfGeneSetCollections))) stop("The listOfGeneSetCollections should be a named list")
#Check that the enrichmentAnalysis object is a list resulting from the geneSetAnalyze function
	if(is.list(enrichmentAnalysis) == FALSE) stop("The enrichmentAnalysis object should be a list")
	if(length(intersect(names(enrichmentAnalysis),c("HyperGeo.results","GSEA.results","Sig.adj.pvals.in.both"))) != 3){
		stop("The enrichmentAnalysis object does not contain the right elements (see help(geneSetAnalyze))")
		}
#Check that nGseaPlots is a single number
	if(length(nGseaPlots) != 1 | (class(nGseaPlots) != "integer" && class(nGseaPlots) != "numeric")) stop("nGseaPlots should be a single numerical value")
#Check that nGseaPlots is not zero and is not bigger than the smaller gene set collection
	minGscSize<-min(sapply(enrichmentAnalysis$GSEA.results,FUN=dim)[1])
	if(nGseaPlots == 0 | nGseaPlots > minGscSize) stop(paste("The value of nGseaPlots should be between 1 and ",minGscSize))
#Produce the GSEA plots, separately for the GO and KEGG gene sets because the naming of the graphs is slightly different
#since the enrichmentAnalysis wrapper function adds the gene set name to the GO/KEGG IDs
	number.gsc<-length(listOfGeneSetCollections)
	for(s in 1:number.gsc){
#1.for the non-GO and non-KEGG gene set collections	
		if(length(intersect(s,keggGS)) == 0 && length(intersect(s,GOGS)) == 0){
			for(p in 1:nGseaPlots){
				t.name<-rownames(enrichmentAnalysis$GSEA.results[[names(listOfGeneSetCollections)[s]]])[p]
				while(length(grep(pattern="/",x=t.name)) != 0){
					t.name<-sub(pattern="/",replacement="_",x=t.name,perl=TRUE)
					}
				test<-gseaScores(geneList=geneList,geneSet=geneIds(listOfGeneSetCollections[[s]][[rownames(enrichmentAnalysis$GSEA.results[[names(listOfGeneSetCollections)[s]]])[p]]]),exponent=1,mode="graph")
				gseaPlots(Running.Score=test$Running.Score,Enrichment.Score=test$Enrichment.Score,Positions=test$Positions,geneList=geneList,output="png",name=t.name)
			}
		}	
#2.for the GO gene set collections		
		if(length(intersect(s,GOGS)) != 0){
			for(p in 1:nGseaPlots){
				GOgsname<-sub(pattern="(\\D*$)",replacement="",x=rownames(enrichmentAnalysis$GSEA.results[[names(listOfGeneSetCollections)[s]]])[p],perl=TRUE)
				test<-gseaScores(geneList=geneList,geneSet=geneIds(listOfGeneSetCollections[[s]][[GOgsname]]),exponent=1,mode="graph");
				gseaPlots(Running.Score=test$Running.Score,Enrichment.Score=test$Enrichment.Score,Positions=test$Positions,geneList=geneList,output="png",name=GOgsname);
				}
		}	
#3.for the KEGG gene set collections		
		if(length(intersect(s,keggGS)) != 0){
			for(p in 1:nGseaPlots){
				kegggsname<-sub(pattern="(\\d*$)",replacement="",x=rownames(enrichmentAnalysis$GSEA.results[[names(listOfGeneSetCollections)[s]]])[p],perl=TRUE)
				kegggsname<-sub(pattern="(\\D*$)",replacement="",x=kegggsname,perl=TRUE)
				kegggsname<-sub(pattern="( :: \\D*$)",replacement="",x=kegggsname,perl=TRUE)
				kegggsname<-sub(pattern="( :: \\w*$)",replacement="",x=kegggsname,perl=TRUE)
				test<-gseaScores(geneList=geneList,geneSet=geneIds(listOfGeneSetCollections[[s]][[kegggsname]]),exponent=1,mode="graph");
				gseaPlots(Running.Score=test$Running.Score,Enrichment.Score=test$Enrichment.Score,Positions=test$Positions,geneList=geneList,output="png",name=kegggsname);
			}
		}	
}
#Also, produce files containing the overlap between hits and gene sets, to be linked to the hypergeometric results pages
	if(species == "Drosophila_melanogaster"){mapID<-as.list(org.Dm.egSYMBOL)}
	if(species == "Homo_sapiens"){mapID<-as.list(org.Hs.egSYMBOL)}
	if(species == "Rattus_norvegicus"){mapID<-as.list(org.Rn.egSYMBOL)}
	if(species == "Mus_musculus"){mapID<-as.list(org.Mm.egSYMBOL)}
	if(species == "Caenorhabditis_elegans"){mapID<-as.list(org.Ce.egSYMBOL)}
	for(s in 1:number.gsc){
#1.for the non-GO and non-KEGG gene set collections	
		if(length(intersect(s,keggGS)) == 0 && length(intersect(s,GOGS)) == 0){
			for(p in 1:nGseaPlots){
				t.name<-rownames(enrichmentAnalysis$HyperGeo.results[[names(listOfGeneSetCollections)[s]]])[p]
				while(length(grep(pattern="/",x=t.name)) != 0){
					t.name<-sub(pattern="/",replacement="_",x=t.name,perl=TRUE)
					}
				overlap<-intersect(geneIds(listOfGeneSetCollections[[s]][[rownames(enrichmentAnalysis$HyperGeo.results[[names(listOfGeneSetCollections)[s]]])[p]]]),hits)
				overlapSymbols<-rep(0,length(overlap))
				if(is.list(mapID) == TRUE) {
					for(g in 1:length(overlapSymbols)){overlapSymbols[g]<-mapID[[overlap[g]]]}
					}
				filename<-paste(t.name,".txt",sep="")
				overlap<-cbind(overlap,overlapSymbols)
				colnames(overlap)<-c("EntrezGene","Symbols")
				write.table(overlap,file=filename,row.names=FALSE,quote=FALSE)
			}
		}	
#2.for the GO gene set collections		
		if(length(intersect(s,GOGS)) != 0){
			for(p in 1:nGseaPlots){
				GOgsname<-sub(pattern="(\\D*$)",replacement="",x=rownames(enrichmentAnalysis$HyperGeo.results[[names(listOfGeneSetCollections)[s]]])[p],perl=TRUE)
				overlap<-intersect(geneIds(listOfGeneSetCollections[[s]][[GOgsname]]),hits)
				overlapSymbols<-rep(0,length(overlap))
				if(is.list(mapID) == TRUE) {
					for(g in 1:length(overlapSymbols)){overlapSymbols[g]<-mapID[[overlap[g]]]}
					}
				filename<-paste(GOgsname,".txt",sep="")
				overlap<-cbind(overlap,overlapSymbols)
				colnames(overlap)<-c("EntrezGene","Symbols")
				write.table(overlap,file=filename,row.names=FALSE,quote=FALSE)
			}
		}	
#3.for the KEGG gene set collections		
		if(length(intersect(s,keggGS)) != 0){
			for(p in 1:nGseaPlots){
				kegggsname<-sub(pattern="(\\d*$)",replacement="",x=rownames(enrichmentAnalysis$HyperGeo.results[[names(listOfGeneSetCollections)[s]]])[p],perl=TRUE)
				kegggsname<-sub(pattern="(\\D*$)",replacement="",x=kegggsname,perl=TRUE)
				kegggsname<-sub(pattern="( :: \\D*$)",replacement="",x=kegggsname,perl=TRUE)
				kegggsname<-sub(pattern="( :: \\w*$)",replacement="",x=kegggsname,perl=TRUE)
				overlap<-intersect(geneIds(listOfGeneSetCollections[[s]][[kegggsname]]),hits)
				overlapSymbols<-rep(0,length(overlap))
				if(is.list(mapID) == TRUE) {
					for(g in 1:length(overlapSymbols)){overlapSymbols[g]<-mapID[[overlap[g]]]}
					}
				filename<-paste(kegggsname,".txt",sep="")
				overlap<-cbind(overlap,overlapSymbols)
				colnames(overlap)<-c("EntrezGene","Symbols")
				write.table(overlap,file=filename,row.names=FALSE,quote=FALSE)
			}
		}	
	}	
#The gseaScores function will check:
	#that the geneList has the right format
	#that the gene list has been named properly
	#that the exponent is a single integer
	#that the geneSet is a single vector	
#Create a report directory and copy css and logos in there
	dir.create("HTSanalyzerReport")
	cpfile<-dir(system.file("templates",package="HTSanalyzeR"),full=TRUE)
	reportdir<-file.path(getwd(),"HTSanalyzerReport")
	file.copy(from=cpfile,to=reportdir,overwrite=TRUE)
#Produce the index html page
	htmlfile=file.path(reportdir,"index.html")
	cat('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"http://www.w3.org/TR/html4/loose.dtd">', file = htmlfile)
	cat('\n <html> \n <link rel="stylesheet" href="./htsanalyzer.css" type="text/css"> ',append = TRUE, file = htmlfile)
	cat('\n <head> <title> HTSanalyzeR Experiment Report </title> </head>',append = TRUE, file = htmlfile)
	cat('\n <body> \n  <table class="border"> \n <tr class="border top"> \n <td class="border corner"> &nbsp&nbsp&nbsp&nbsp </td>',append = TRUE,file=htmlfile)
	cat(paste('\n <td class="border top"> <div class="header"> Report for Experiment <span class="header"> ',experimentName,'</span> </div>',sep=""),append = TRUE,file=htmlfile)
	cat(paste('<div class="timestamp">generated  ',date(),'(<small>version 1.1.7</small>) </div>',sep=""),append = TRUE,file=htmlfile)
	cat('\n <div class="HTSheader"> HTSanalyzeR </div> </td> </tr> <tr class="border middle"> <td class="border left"></td>',append = TRUE,file=htmlfile)
	cat('<td class="main"> <table> <tr> <div class="HTSlogos"> <img src="./Rlogo.png" width="50" height="40"/>',append = TRUE,file=htmlfile)
	cat('&nbsp <img src="./blue_cruklogo.gif" width="120" height="50"/>',append = TRUE,file=htmlfile)
	cat('&nbsp <img src="./goatcomputer.png" width="85" height="60"/> </div></tr><tr>',append = TRUE,file=htmlfile)
#Produce the tabs	
	cat('\n <table class="noframe"> <tr>  ',append = TRUE, file = htmlfile)
	l.HyperGeo.results<-length(enrichmentAnalysis$HyperGeo.results)
	l.GSEA.results<-length(enrichmentAnalysis$GSEA.results)
	l.Sig.adj.pvals.in.both<-length(enrichmentAnalysis$Sig.adj.pvals.in.both)
	cat('<td class="tabs"> <h3>  </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.HyperGeo.results){
		cat(paste('\n <td class="tabs"> <h3>',names(enrichmentAnalysis$HyperGeo.results)[i],'</h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)	
	cat('<td class="tabs"> <h3> Hypergeometric tests </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.HyperGeo.results){
		cat(paste('\n <td class="tabs"> <h3><a href="./hyperg',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
		cat(paste(names(enrichmentAnalysis$HyperGeo.results)[i],' Hyperg. Tests">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)
	cat('<td class="tabs"> <h3> GSEA </h3> </td>', append = TRUE, file = htmlfile)	
	for(i in 1:l.GSEA.results){
		cat(paste('\n <td class="tabs"> <h3><a href="./gsea',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
		cat(paste(names(enrichmentAnalysis$GSEA.results)[i],' GSEA">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)	
	cat('<td class="tabs"> <h3> Enrichment Summary </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.Sig.adj.pvals.in.both){
			cat(paste('\n <td class="tabs"> <h3><a href="./enrichment',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
			cat(paste(names(enrichmentAnalysis$Sig.adj.pvals.in.both)[i],' Enrichment.summary">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
			}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)
	cat(paste('\n  <td class="tabs"> <h3><a href="./network.html" title="network">Network Analysis</a></h3> </td> </tr> </table>',sep=""), append = TRUE, file = htmlfile)
#check that the arguments experimentName, geneListName have length one
	if(length(experimentName) != 1  | length(geneListName) != 1 ) warning("The experimentName and geneListName should be of length one")
	cat("\n <hr/> \n <br> The enrichment analysis was performed using the data: ", append = TRUE, file = htmlfile)
	cat(geneListName, append = TRUE, file = htmlfile)
	cat(paste(" ( ",length(geneList), " genes)",sep=""), append = TRUE, file = htmlfile)
	cat("\n <br> This analysis was performed using the gene set collection(s): ", append = TRUE, file = htmlfile)
	cat("\n \t <UL>", append = TRUE, file = htmlfile)
	for(i in 1:number.gsc){
		cat(paste(" \n \t \t <LI>",names(listOfGeneSetCollections)[i],sep=""), append = TRUE, file = htmlfile)
		cat(paste(" ( ",length(listOfGeneSetCollections[[i]]), " gene sets, of which ",dim(enrichmentAnalysis$GSEA.results[[i]])[1] ," were above the minimum size )",sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n \t </UL>", append = TRUE, file = htmlfile)
	cat("\n <br> The following methods were used: ", append = TRUE, file = htmlfile)
	cat("\n \t <UL> \n \t \t <LI>", append = TRUE, file = htmlfile)
	cat("Hypergeometric test", append = TRUE, file = htmlfile)
	cat("\n \t \t <UL> \n \t \t \t <LI>", append = TRUE, file = htmlfile)
#check that the arguments cutoffHits, p.adj.method, min.gene.set.size, nPermutations have length one	
	if(length(cutoffHits) != 1  | length(p.adj.method) != 1 | length(min.gene.set.size) != 1 | length(nPermutations) != 1) warning("Thearguments cutoffHits, p.adj.method, min.gene.set.size, nPermutations should be of length one")
	cat(paste("Cutoff for hits: ",cutoffHits,sep=""), append = TRUE, file = htmlfile)
	cat(paste("\n \t \t \t <LI> MHT correction method: ",p.adj.method), append = TRUE, file = htmlfile)
	cat(paste("\n \t \t \t <LI> Minimum gene set size: ",min.gene.set.size), append = TRUE, file = htmlfile)
	cat("\n \t \t </UL>", append = TRUE, file = htmlfile)
	cat("\n \t \t <LI>", append = TRUE, file = htmlfile)
	cat("Gene Set Enrichment Analysis", append = TRUE, file = htmlfile)
	cat("\n \t \t <UL> \n \t \t \t <LI>", append = TRUE, file = htmlfile)
	cat(paste("Minimum gene set size: ",min.gene.set.size,sep=""), append = TRUE, file = htmlfile)
	cat(paste("\n \t \t \t <LI> MHT correction method: ",p.adj.method), append = TRUE, file = htmlfile)
	cat(paste("\n \t \t \t <LI> Number of permutations: ",nPermutations), append = TRUE, file = htmlfile)
	cat(paste("\n \t \t \t <LI> Exponent: ",exponent), append = TRUE, file = htmlfile)
	cat("\n \t \t </UL>", append = TRUE, file = htmlfile)
	cat("\n \t \t <LI>", append = TRUE, file = htmlfile)
	cat("Network Analysis", append = TRUE, file = htmlfile)
	cat("\n \t \t <UL> \n \t \t \t <LI>", append = TRUE, file = htmlfile)
	cat(paste("Tests performed: ",tests,"( ",alternative," )",sep=""), append = TRUE, file = htmlfile)
	cat(paste("\n \t \t \t <LI> P-values used: ",columns), append = TRUE, file = htmlfile)
	if(is.na(biogridObject)) {
		cat(paste("\n \t \t \t <LI> Interaction dataset: The Biogrid organism: ",species), append = TRUE, file = htmlfile)
		}else{
		cat(paste("\n \t \t \t <LI> Interaction dataset: Biogrid object: ",deparse(substitute(biogridObject))), append = TRUE, file = htmlfile)
		}
	if(genetic == FALSE) {
		cat(" (excluding genetic interactions)", append = TRUE, file = htmlfile)
		}else{
		cat(" (including genetic interactions)", append = TRUE, file = htmlfile)
		}	
	cat(paste("\n \t \t \t <LI> FDR for score calculation: ",fdr), append = TRUE, file = htmlfile)
	cat("\n \t \t </UL>", append = TRUE, file = htmlfile)
	cat("\n \t </UL>", append = TRUE, file = htmlfile)
	cat('\n </tr> \n </table> \n </td> \n </tr> \n </table> \n </body> \n </html>',append = TRUE, file = htmlfile)
#Produce the HyperG results pages: one page is produced for each element of the enrichmentAnalysis$HyperGeo.results list (i.e. for each gene set collection)
	for(n in 1:l.HyperGeo.results){
		htmlfile =  file.path(reportdir,paste("hyperg",n,".html",sep=""))
#First, take care of the non-GO and non-KEGG gene set collections		
		if(length(intersect(n,keggGS)) == 0 && length(intersect(n,GOGS)) == 0){
			cat('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"http://www.w3.org/TR/html4/loose.dtd">', file = htmlfile)
			cat('\n <html> \n <link rel="stylesheet" href="./htsanalyzer.css" type="text/css"> ',append = TRUE, file = htmlfile)
			cat('\n <head> <title> HTSanalyzeR Experiment Report </title> </head>',append = TRUE, file = htmlfile)
			cat('\n <body> \n  <table class="border"> \n <tr class="border top"> \n <td class="border corner"> &nbsp&nbsp&nbsp&nbsp </td>',append = TRUE,file=htmlfile)
			cat(paste('\n <td class="border top"> <div class="header"> Report for Experiment <span class="header"> ',experimentName,'</span> </div>',sep=""),append = TRUE,file=htmlfile)
			cat(paste('<div class="timestamp">generated  ',date(),'(<small>version 1.1.7</small>) </div>',sep=""),append = TRUE,file=htmlfile)
			cat('\n <div class="HTSheader"> HTSanalyzeR </div> </td> </tr> <tr class="border middle"> <td class="border left"></td>',append = TRUE,file=htmlfile)
			cat('<td class="main"> <table> <tr> <div class="HTSlogos"> <img src="./Rlogo.png" width="50" height="40"/>',append = TRUE,file=htmlfile)
			cat('&nbsp <img src="./blue_cruklogo.gif" width="120" height="50"/>',append = TRUE,file=htmlfile)
			cat('&nbsp <img src="./goatcomputer.png" width="85" height="60"/> </div></tr><tr>',append = TRUE,file=htmlfile)
#Produce the tabs
	cat(paste('\n <table class="noframe"> <tr> <td class="tabs"> <h3><a href="./index.html" title="index">Index</a></h3> </td> <tr>',sep=""), append = TRUE, file = htmlfile)
	cat("\n <tr>", append = TRUE, file = htmlfile)
	l.HyperGeo.results<-length(enrichmentAnalysis$HyperGeo.results)
	l.GSEA.results<-length(enrichmentAnalysis$GSEA.results)
	l.Sig.adj.pvals.in.both<-length(enrichmentAnalysis$Sig.adj.pvals.in.both)
	cat('<td class="tabs"> <h3>  </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.HyperGeo.results){
		cat(paste('\n <td class="tabs"> <h3>',names(enrichmentAnalysis$HyperGeo.results)[i],'</h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)	
	cat('<td class="tabs"> <h3> Hypergeometric tests </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.HyperGeo.results){
		cat(paste('\n <td class="tabs"> <h3><a href="./hyperg',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
		cat(paste(names(enrichmentAnalysis$HyperGeo.results)[i],' Hyperg. Tests">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)
	cat('<td class="tabs"> <h3> GSEA </h3> </td>', append = TRUE, file = htmlfile)	
	for(i in 1:l.GSEA.results){
		cat(paste('\n <td class="tabs"> <h3><a href="./gsea',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
		cat(paste(names(enrichmentAnalysis$GSEA.results)[i],' GSEA">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)	
	cat('<td class="tabs"> <h3> Enrichment Summary </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.Sig.adj.pvals.in.both){
			cat(paste('\n <td class="tabs"> <h3><a href="./enrichment',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
			cat(paste(names(enrichmentAnalysis$Sig.adj.pvals.in.both)[i],' Enrichment.summary">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
			}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)
	cat(paste('\n  <td class="tabs"> <h3><a href="./network.html" title="network">Network Analysis</a></h3> </td> </tr> </table>',sep=""), append = TRUE, file = htmlfile)
#Produce the main data frame: headers and first row	
		t.name<-rownames(enrichmentAnalysis$HyperGeo.results[[n]])[1]
		while(length(grep(pattern="/",x=t.name)) != 0){
					t.name<-sub(pattern="/",replacement="_",x=t.name,perl=TRUE)
					}	
		cat(paste('\n <hr/> \n <br>', names(enrichmentAnalysis$HyperGeo.results)[n],' Hyperg. Tests <br>','\n <table class="results"> <tr class="head"> <th>Gene Set name </th> \n <th>',
		colnames(enrichmentAnalysis$HyperGeo.results[[n]])[1],"</th> \n <th>",
		colnames(enrichmentAnalysis$HyperGeo.results[[n]])[2],"</th> \n <th>",
		colnames(enrichmentAnalysis$HyperGeo.results[[n]])[3],"</th> \n <th>",
		colnames(enrichmentAnalysis$HyperGeo.results[[n]])[4],"</th> \n <th>",
		colnames(enrichmentAnalysis$HyperGeo.results[[n]])[5],"</th> \n <th>",
		colnames(enrichmentAnalysis$HyperGeo.results[[n]])[6],"</th> \n <th>",
		colnames(enrichmentAnalysis$HyperGeo.results[[n]])[8],"</th> \n </tr>",
		'<tr class="odd"> <td>',rownames(enrichmentAnalysis$HyperGeo.results[[n]])[1],"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][1,1], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][1,2], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][1,3], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][1,4], digits=4),"</td> \n <td>", sep=""), append = TRUE, file = htmlfile)
		cat(paste('<a href="../',t.name,'.txt" target="_blank" title="Observed.hits">',signif(enrichmentAnalysis$HyperGeo.results[[n]][1,5], digits=4),'</a> </td> \n <td>', sep=""), append = TRUE, file = htmlfile)
		cat(paste(signif(enrichmentAnalysis$HyperGeo.results[[n]][1,6], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][1,8], digits=4),"</td> \n </tr>", sep=""), append = TRUE, file = htmlfile)
#Produce the main dataframe: rows with plots
	for(i in 2:nGseaPlots){
			t.name<-rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i]
			while(length(grep(pattern="/",x=t.name)) != 0){
					t.name<-sub(pattern="/",replacement="_",x=t.name,perl=TRUE)
					}
	if((i%%2) == 0){
		cat(paste('<tr class="even"> <td>',rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,1], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,2], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,3], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,4], digits=4),"</td> \n <td>", sep=""), append = TRUE, file = htmlfile)
		cat(paste('<a href="../',t.name,'.txt" target="_blank" title="Observed.hits">',signif(enrichmentAnalysis$HyperGeo.results[[n]][i,5], digits=4),'</a> </td> \n <td>', sep=""), append = TRUE, file = htmlfile)
		cat(paste(signif(enrichmentAnalysis$HyperGeo.results[[n]][i,6], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,8], digits=4),"</td> \n </tr>", sep=""), append = TRUE, file = htmlfile)
		}else{
		cat(paste('<tr class="odd"> <td>',rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,1], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,2], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,3], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,4], digits=4),"</td> \n <td>", sep=""), append = TRUE, file = htmlfile)
		cat(paste('<a href="../',t.name,'.txt" target="_blank" title="Observed.hits">',signif(enrichmentAnalysis$HyperGeo.results[[n]][i,5], digits=4),'</a> </td> \n <td>', sep=""), append = TRUE, file = htmlfile)
		cat(paste(signif(enrichmentAnalysis$HyperGeo.results[[n]][i,6], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,8], digits=4),"</td> \n </tr>", sep=""), append = TRUE, file = htmlfile)
		}	
	}
#Produce the main data frame: rows without plots			
	for(i in (nGseaPlots+1):dim(enrichmentAnalysis$HyperGeo.results[[n]])[1]){
	if((i%%2) == 0){
		cat(paste('<tr class="even"> <td>',rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,1], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,2], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,3], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,4], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,5], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,6], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,8], digits=4),"</td> </tr>", sep=""),append = TRUE, file = htmlfile)
		}else{
		cat(paste('<tr class="odd"> <td>',rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,1], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,2], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,3], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,4], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,5], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,6], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,8], digits=4),"</td> </tr>", sep=""),append = TRUE, file = htmlfile)
		}	
	}
	cat("</table>",append = TRUE, file = htmlfile)	
	cat('\n </tr> \n </table> \n </td> \n </tr> \n </table> \n </body> \n </html>',append = TRUE, file = htmlfile)
	}
#Second, take care of the KEGG gene set collections	
	if(length(intersect(n,keggGS)) != 0){
			cat('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"http://www.w3.org/TR/html4/loose.dtd">', file = htmlfile)
			cat('\n <html> \n <link rel="stylesheet" href="./htsanalyzer.css" type="text/css"> ',append = TRUE, file = htmlfile)
			cat('\n <head> <title> HTSanalyzeR Experiment Report </title> </head>',append = TRUE, file = htmlfile)
			cat('\n <body> \n  <table class="border"> \n <tr class="border top"> \n <td class="border corner"> &nbsp&nbsp&nbsp&nbsp </td>',append = TRUE,file=htmlfile)
			cat(paste('\n <td class="border top"> <div class="header"> Report for Experiment <span class="header"> ',experimentName,'</span> </div>',sep=""),append = TRUE,file=htmlfile)
			cat(paste('<div class="timestamp">generated  ',date(),'(<small>version 1.1.7</small>) </div>',sep=""),append = TRUE,file=htmlfile)
			cat('\n <div class="HTSheader"> HTSanalyzeR </div> </td> </tr> <tr class="border middle"> <td class="border left"></td>',append = TRUE,file=htmlfile)
			cat('<td class="main"> <table> <tr> <div class="HTSlogos"> <img src="./Rlogo.png" width="50" height="40"/>',append = TRUE,file=htmlfile)
			cat('&nbsp <img src="./blue_cruklogo.gif" width="120" height="50"/>',append = TRUE,file=htmlfile)
			cat('&nbsp <img src="./goatcomputer.png" width="85" height="60"/> </div></tr><tr>',append = TRUE,file=htmlfile)
#Produce the tabs
	cat(paste('\n <table class="noframe"> <tr> <td class="tabs"> <h3><a href="./index.html" title="index">Index</a></h3> </td> <tr>',sep=""), append = TRUE, file = htmlfile)
	cat("\n <tr>", append = TRUE, file = htmlfile)
	l.HyperGeo.results<-length(enrichmentAnalysis$HyperGeo.results)
	l.GSEA.results<-length(enrichmentAnalysis$GSEA.results)
	l.Sig.adj.pvals.in.both<-length(enrichmentAnalysis$Sig.adj.pvals.in.both)
	cat('<td class="tabs"> <h3>  </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.HyperGeo.results){
		cat(paste('\n <td class="tabs"> <h3>',names(enrichmentAnalysis$HyperGeo.results)[i],'</h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)	
	cat('<td class="tabs"> <h3> Hypergeometric tests </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.HyperGeo.results){
		cat(paste('\n <td class="tabs"> <h3><a href="./hyperg',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
		cat(paste(names(enrichmentAnalysis$HyperGeo.results)[i],' Hyperg. Tests">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)
	cat('<td class="tabs"> <h3> GSEA </h3> </td>', append = TRUE, file = htmlfile)	
	for(i in 1:l.GSEA.results){
		cat(paste('\n <td class="tabs"> <h3><a href="./gsea',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
		cat(paste(names(enrichmentAnalysis$GSEA.results)[i],' GSEA">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)	
	cat('<td class="tabs"> <h3> Enrichment Summary </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.Sig.adj.pvals.in.both){
			cat(paste('\n <td class="tabs"> <h3><a href="./enrichment',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
			cat(paste(names(enrichmentAnalysis$Sig.adj.pvals.in.both)[i],' Enrichment.summary">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
			}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)
	cat(paste('\n  <td class="tabs"> <h3><a href="./network.html" title="network">Network Analysis</a></h3> </td> </tr> </table>',sep=""), append = TRUE, file = htmlfile)
#Produce the main data frame
		kegggsname<-sub(pattern="(\\d*$)",replacement="",x=rownames(enrichmentAnalysis$HyperGeo.results[[n]])[1],perl=TRUE)
		kegggsname<-sub(pattern="(\\D*$)",replacement="",x=kegggsname,perl=TRUE)
		kegggsname<-sub(pattern="( :: \\D*$)",replacement="",x=kegggsname,perl=TRUE)
 		kegggsname<-sub(pattern="( :: \\w*$)",replacement="",x=kegggsname,perl=TRUE)
		cat(paste('\n <hr/> \n <br>', names(enrichmentAnalysis$HyperGeo.results)[n],' Hyperg. Tests <br>','\n <table class="results"> <tr class="head"> <th>Gene Set name </th> \n <th>',
		colnames(enrichmentAnalysis$HyperGeo.results[[n]])[1],"</th> \n <th>",
		colnames(enrichmentAnalysis$HyperGeo.results[[n]])[2],"</th> \n <th>",
		colnames(enrichmentAnalysis$HyperGeo.results[[n]])[3],"</th> \n <th>",
		colnames(enrichmentAnalysis$HyperGeo.results[[n]])[4],"</th> \n <th>",
		colnames(enrichmentAnalysis$HyperGeo.results[[n]])[5],"</th> \n <th>",
		colnames(enrichmentAnalysis$HyperGeo.results[[n]])[6],"</th> \n <th>",
		colnames(enrichmentAnalysis$HyperGeo.results[[n]])[8],"</th> \n </tr>",
		'<tr class="odd"> <td> <a href="http://www.genome.jp/dbget-bin/www_bget?pathway:', kegggsname,'" target="_blank" title="',rownames(enrichmentAnalysis$HyperGeo.results[[n]])[1],'">',rownames(enrichmentAnalysis$HyperGeo.results[[n]])[1],'</a></td> \n <td>',
		signif(enrichmentAnalysis$HyperGeo.results[[n]][1,1], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][1,2], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][1,3], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][1,4], digits=4),"</td> \n <td>", sep=""), append = TRUE, file = htmlfile)
		cat(paste('<a href="../',kegggsname,'.txt" target="_blank" title="Observed.hits">',signif(enrichmentAnalysis$HyperGeo.results[[n]][1,5], digits=4),'</a> </td> \n <td>', sep=""), append = TRUE, file = htmlfile)
		cat(paste(signif(enrichmentAnalysis$HyperGeo.results[[n]][1,6], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][1,8], digits=4),"</td> \n </tr>", sep=""), append = TRUE, file = htmlfile)
#Rows with a plot		
	for(i in 2:nGseaPlots){
	kegggsname<-sub(pattern="(\\d*$)",replacement="",x=rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],perl=TRUE)
	kegggsname<-sub(pattern="(\\D*$)",replacement="",x=kegggsname,perl=TRUE)
	kegggsname<-sub(pattern="( :: \\D*$)",replacement="",x=kegggsname,perl=TRUE)
	kegggsname<-sub(pattern="( :: \\w*$)",replacement="",x=kegggsname,perl=TRUE)
	if((i%%2) == 0){
		cat(paste('<tr class="even"> <td> <a href="http://www.genome.jp/dbget-bin/www_bget?pathway:',kegggsname,'" target="_blank" title="',rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],'">',rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],'</a> </td> \n <td>',
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,1], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,2], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,3], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,4], digits=4),"</td> \n <td>", sep=""), append = TRUE, file = htmlfile)
		cat(paste('<a href="../',kegggsname,'.txt" target="_blank" title="Observed.hits">',signif(enrichmentAnalysis$HyperGeo.results[[n]][i,5], digits=4),'</a> </td> \n <td>', sep=""), append = TRUE, file = htmlfile)
		cat(paste(signif(enrichmentAnalysis$HyperGeo.results[[n]][i,6], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,8], digits=4),"</td> \n </tr>", sep=""), append = TRUE, file = htmlfile)
		}else{
		cat(paste('<tr class="odd"> <td><a href="http://www.genome.jp/dbget-bin/www_bget?pathway:',kegggsname,'" target="_blank" title="',rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],'">',rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],'</a> </td> \n <td>',
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,1], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,2], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,3], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,4], digits=4),"</td> \n <td>", sep=""), append = TRUE, file = htmlfile)
		cat(paste('<a href="../',kegggsname,'.txt" target="_blank" title="Observed.hits">',signif(enrichmentAnalysis$HyperGeo.results[[n]][i,5], digits=4),'</a> </td> \n <td>', sep=""), append = TRUE, file = htmlfile)
		cat(paste(signif(enrichmentAnalysis$HyperGeo.results[[n]][i,6], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,8], digits=4),"</td> \n </tr>", sep=""), append = TRUE, file = htmlfile)
		}	
	}
#Rows without plot		
	for(i in (nGseaPlots+1):dim(enrichmentAnalysis$HyperGeo.results[[n]])[1]){
	kegggsname<-sub(pattern="(\\d*$)",replacement="",x=rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],perl=TRUE)
	kegggsname<-sub(pattern="(\\D*$)",replacement="",x=kegggsname,perl=TRUE)
	kegggsname<-sub(pattern="( :: \\D*$)",replacement="",x=kegggsname,perl=TRUE)
	kegggsname<-sub(pattern="( :: \\w*$)",replacement="",x=kegggsname,perl=TRUE)
	if((i%%2) == 0){
		cat(paste('<tr class="even"> <td> <a href="http://www.genome.jp/dbget-bin/www_bget?pathway:',kegggsname,'" target="_blank" title="',rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],'">',rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],'</a> </td> \n <td>',
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,1], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,2], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,3], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,4], digits=4),"</td> \n <td>", 
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,5], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,6], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,8], digits=4),"</td> \n </tr>", sep=""), append = TRUE, file = htmlfile)
		}else{
		cat(paste('<tr class="odd"> <td><a href="http://www.genome.jp/dbget-bin/www_bget?pathway:',kegggsname,'" target="_blank" title="',rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],'">',rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],'</a> </td> \n <td>',
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,1], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,2], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,3], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,4], digits=4),"</td> \n <td>", 
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,5], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,6], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,8], digits=4),"</td> \n </tr>", sep=""), append = TRUE, file = htmlfile)
		}	
	}	
	cat("</table>",append = TRUE, file = htmlfile)	
	cat('\n </tr> \n </table> \n </td> \n </tr> \n </table> \n </body> \n </html>',append = TRUE, file = htmlfile)
	}
#Third, take care of the GO gene set collections	
	if(length(intersect(n,GOGS)) != 0){
			cat('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"http://www.w3.org/TR/html4/loose.dtd">', file = htmlfile)
			cat('\n <html> \n <link rel="stylesheet" href="./htsanalyzer.css" type="text/css"> ',append = TRUE, file = htmlfile)
			cat('\n <head> <title> HTSanalyzeR Experiment Report </title> </head>',append = TRUE, file = htmlfile)
			cat('\n <body> \n  <table class="border"> \n <tr class="border top"> \n <td class="border corner"> &nbsp&nbsp&nbsp&nbsp </td>',append = TRUE,file=htmlfile)
			cat(paste('\n <td class="border top"> <div class="header"> Report for Experiment <span class="header"> ',experimentName,'</span> </div>',sep=""),append = TRUE,file=htmlfile)
			cat(paste('<div class="timestamp">generated  ',date(),'(<small>version 1.1.7</small>) </div>',sep=""),append = TRUE,file=htmlfile)
			cat('\n <div class="HTSheader"> HTSanalyzeR </div> </td> </tr> <tr class="border middle"> <td class="border left"></td>',append = TRUE,file=htmlfile)
			cat('<td class="main"> <table> <tr> <div class="HTSlogos"> <img src="./Rlogo.png" width="50" height="40"/>',append = TRUE,file=htmlfile)
			cat('&nbsp <img src="./blue_cruklogo.gif" width="120" height="50"/>',append = TRUE,file=htmlfile)
			cat('&nbsp <img src="./goatcomputer.png" width="85" height="60"/> </div></tr><tr>',append = TRUE,file=htmlfile)
#Produce the tabs
	cat(paste('\n <table class="noframe"> <tr> <td class="tabs"> <h3><a href="./index.html" title="index">Index</a></h3> </td> <tr>',sep=""), append = TRUE, file = htmlfile)
	cat("\n <tr>", append = TRUE, file = htmlfile)
	l.HyperGeo.results<-length(enrichmentAnalysis$HyperGeo.results)
	l.GSEA.results<-length(enrichmentAnalysis$GSEA.results)
	l.Sig.adj.pvals.in.both<-length(enrichmentAnalysis$Sig.adj.pvals.in.both)
	cat('<td class="tabs"> <h3>  </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.HyperGeo.results){
		cat(paste('\n <td class="tabs"> <h3>',names(enrichmentAnalysis$HyperGeo.results)[i],'</h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)	
	cat('<td class="tabs"> <h3> Hypergeometric tests </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.HyperGeo.results){
		cat(paste('\n <td class="tabs"> <h3><a href="./hyperg',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
		cat(paste(names(enrichmentAnalysis$HyperGeo.results)[i],' Hyperg. Tests">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)
	cat('<td class="tabs"> <h3> GSEA </h3> </td>', append = TRUE, file = htmlfile)	
	for(i in 1:l.GSEA.results){
		cat(paste('\n <td class="tabs"> <h3><a href="./gsea',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
		cat(paste(names(enrichmentAnalysis$GSEA.results)[i],' GSEA">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)	
	cat('<td class="tabs"> <h3> Enrichment Summary </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.Sig.adj.pvals.in.both){
			cat(paste('\n <td class="tabs"> <h3><a href="./enrichment',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
			cat(paste(names(enrichmentAnalysis$Sig.adj.pvals.in.both)[i],' Enrichment.summary">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
			}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)
	cat(paste('\n  <td class="tabs"> <h3><a href="./network.html" title="network">Network Analysis</a></h3> </td> </tr> </table>',sep=""), append = TRUE, file = htmlfile)
#Produce the main dataframe: headers and first row		
		cat(paste('\n <hr/> \n <br>', names(enrichmentAnalysis$HyperGeo.results)[n],' Hyperg. Tests <br>','\n <table class="results"> <tr class="head"> <th>Gene Set name </th> \n <th>',
		colnames(enrichmentAnalysis$HyperGeo.results[[n]])[1],"</th> \n <th>",
		colnames(enrichmentAnalysis$HyperGeo.results[[n]])[2],"</th> \n <th>",
		colnames(enrichmentAnalysis$HyperGeo.results[[n]])[3],"</th> \n <th>",
		colnames(enrichmentAnalysis$HyperGeo.results[[n]])[4],"</th> \n <th>",
		colnames(enrichmentAnalysis$HyperGeo.results[[n]])[5],"</th> \n <th>",
		colnames(enrichmentAnalysis$HyperGeo.results[[n]])[6],"</th> \n <th>",
		colnames(enrichmentAnalysis$HyperGeo.results[[n]])[8],"</th> \n </tr>",
		'<tr class="odd"> <td> <a href="http://www.ebi.ac.uk/QuickGO/GTerm?id=GO:',sub(pattern="(\\D*$)",replacement="",x=sub(pattern="(\\D*)",replacement="",x=rownames(enrichmentAnalysis$HyperGeo.results[[n]])[1],perl=TRUE),perl=TRUE),'" target="_blank" title="',rownames(enrichmentAnalysis$HyperGeo.results[[n]])[1],'">',rownames(enrichmentAnalysis$HyperGeo.results[[n]])[1],'</a></td> \n <td>',
		signif(enrichmentAnalysis$HyperGeo.results[[n]][1,1], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][1,2], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][1,3], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][1,4], digits=4),"</td> \n <td>", sep=""), append = TRUE, file = htmlfile)
		cat(paste('<a href="../',
		sub(pattern="(\\D*$)",replacement="",x=rownames(enrichmentAnalysis$HyperGeo.results[[n]])[1],perl=TRUE),
		'.txt" target="_blank" title="Observed.hits">',signif(enrichmentAnalysis$HyperGeo.results[[n]][1,5], digits=4),'</a> </td> \n <td>', sep=""), append = TRUE, file = htmlfile)
		cat(paste(signif(enrichmentAnalysis$HyperGeo.results[[n]][1,6], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][1,8], digits=4),"</td> \n </tr>", sep=""), append = TRUE, file = htmlfile)
#Produce the main dataframe: rows with a plot
	for(i in 2:nGseaPlots){
	if((i%%2) == 0){
		cat(paste('<tr class="even"> <td> <a href="http://www.ebi.ac.uk/QuickGO/GTerm?id=GO:',sub(pattern="(\\D*$)",replacement="",x=sub(pattern="(\\D*)",replacement="",x=rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],perl=TRUE),perl=TRUE),'" target="_blank" title="',rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],'">',rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],'</a> </td> \n <td>',
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,1], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,2], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,3], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,4], digits=4),"</td> \n <td>", sep=""), append = TRUE, file = htmlfile)
		cat(paste('<a href="../',
		sub(pattern="(\\D*$)",replacement="",x=rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],perl=TRUE),
		'.txt" target="_blank" title="Observed.hits">',signif(enrichmentAnalysis$HyperGeo.results[[n]][i,5], digits=4),'</a> </td> \n <td>', sep=""), append = TRUE, file = htmlfile)
		cat(paste(signif(enrichmentAnalysis$HyperGeo.results[[n]][i,6], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,8], digits=4),"</td> \n </tr>", sep=""), append = TRUE, file = htmlfile)
		}else{
		cat(paste('<tr class="odd"> <td><a href="http://www.ebi.ac.uk/QuickGO/GTerm?id=GO:',sub(pattern="(\\D*$)",replacement="",x=sub(pattern="(\\D*)",replacement="",x=rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],perl=TRUE),perl=TRUE),'" target="_blank" title="',rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],'">',rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],'</a> </td> \n <td>',
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,1], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,2], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,3], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,4], digits=4),"</td> \n <td>", sep=""), append = TRUE, file = htmlfile)
		cat(paste('<a href="../',
		sub(pattern="(\\D*$)",replacement="",x=rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],perl=TRUE),
		'.txt" target="_blank" title="Observed.hits">',signif(enrichmentAnalysis$HyperGeo.results[[n]][i,5], digits=4),'</a> </td> \n <td>', sep=""), append = TRUE, file = htmlfile)
		cat(paste(signif(enrichmentAnalysis$HyperGeo.results[[n]][i,6], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,8], digits=4),"</td> \n </tr>", sep=""), append = TRUE, file = htmlfile)
		}	
	}
#Produce the main dataframe: rows without plot
	for(i in (nGseaPlots+1): dim(enrichmentAnalysis$HyperGeo.results[[n]])[1]){
	if((i%%2) == 0){
		cat(paste('<tr class="even"> <td> <a href="http://www.ebi.ac.uk/QuickGO/GTerm?id=GO:',sub(pattern="(\\D*$)",replacement="",x=sub(pattern="(\\D*)",replacement="",x=rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],perl=TRUE),perl=TRUE),'" target="_blank" title="',rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],'">',rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],'</a> </td> \n <td>',
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,1], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,2], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,3], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,4], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,5], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,6], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,8], digits=4),"</td> \n </tr>", sep=""), append = TRUE, file = htmlfile)
		}else{
		cat(paste('<tr class="odd"> <td><a href="http://www.ebi.ac.uk/QuickGO/GTerm?id=GO:',sub(pattern="(\\D*$)",replacement="",x=sub(pattern="(\\D*)",replacement="",x=rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],perl=TRUE),perl=TRUE),'" target="_blank" title="',rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],'">',rownames(enrichmentAnalysis$HyperGeo.results[[n]])[i],'</a> </td> \n <td>',
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,1], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,2], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,3], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,4], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,5], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,6], digits=4),"</td> \n <td>",
		signif(enrichmentAnalysis$HyperGeo.results[[n]][i,8], digits=4),"</td> \n </tr>", sep=""), append = TRUE, file = htmlfile)
		}	
	}	
	cat("</table>",append = TRUE, file = htmlfile)	
	cat('\n </tr> \n </table> \n </td> \n </tr> \n </table> \n </body> \n </html>',append = TRUE, file = htmlfile)
	}
	}
#Produce the gsea results pages:  one page is produced for each element of the enrichmentAnalysis$GSEA.results list (i.e. for each gene set collection)
	for(n in 1:l.GSEA.results){
		htmlfile =  file.path(reportdir,paste("gsea",n,".html",sep=""))
#First, take care of the non-GO and non-KEGG gene set collections		
		if(length(intersect(n,keggGS)) == 0 && length(intersect(n,GOGS)) == 0){
			cat('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"http://www.w3.org/TR/html4/loose.dtd">', file = htmlfile)
			cat('\n <html> \n <link rel="stylesheet" href="./htsanalyzer.css" type="text/css"> ',append = TRUE, file = htmlfile)
			cat('\n <head> <title> HTSanalyzeR Experiment Report </title> </head>',append = TRUE, file = htmlfile)
			cat('\n <body> \n  <table class="border"> \n <tr class="border top"> \n <td class="border corner"> &nbsp&nbsp&nbsp&nbsp </td>',append = TRUE,file=htmlfile)
			cat(paste('\n <td class="border top"> <div class="header"> Report for Experiment <span class="header"> ',experimentName,'</span> </div>',sep=""),append = TRUE,file=htmlfile)
			cat(paste('<div class="timestamp">generated  ',date(),'(<small>version 1.1.7</small>) </div>',sep=""),append = TRUE,file=htmlfile)
			cat('\n <div class="HTSheader"> HTSanalyzeR </div> </td> </tr> <tr class="border middle"> <td class="border left"></td>',append = TRUE,file=htmlfile)
			cat('<td class="main"> <table> <tr> <div class="HTSlogos"> <img src="./Rlogo.png" width="50" height="40"/>',append = TRUE,file=htmlfile)
			cat('&nbsp <img src="./blue_cruklogo.gif" width="120" height="50"/>',append = TRUE,file=htmlfile)
			cat('&nbsp <img src="./goatcomputer.png" width="85" height="60"/> </div></tr><tr>',append = TRUE,file=htmlfile)
#Produce the tabs
	cat(paste('\n <table class="noframe"> <tr> <td class="tabs"> <h3><a href="./index.html" title="index">Index</a></h3> </td> <tr>',sep=""), append = TRUE, file = htmlfile)
	cat("\n <tr>", append = TRUE, file = htmlfile)
	l.HyperGeo.results<-length(enrichmentAnalysis$HyperGeo.results)
	l.GSEA.results<-length(enrichmentAnalysis$GSEA.results)
	l.Sig.adj.pvals.in.both<-length(enrichmentAnalysis$Sig.adj.pvals.in.both)
	cat('<td class="tabs"> <h3>  </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.HyperGeo.results){
		cat(paste('\n <td class="tabs"> <h3>',names(enrichmentAnalysis$HyperGeo.results)[i],'</h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)	
	cat('<td class="tabs"> <h3> Hypergeometric tests </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.HyperGeo.results){
		cat(paste('\n <td class="tabs"> <h3><a href="./hyperg',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
		cat(paste(names(enrichmentAnalysis$HyperGeo.results)[i],' Hyperg. Tests">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)
	cat('<td class="tabs"> <h3> GSEA </h3> </td>', append = TRUE, file = htmlfile)	
	for(i in 1:l.GSEA.results){
		cat(paste('\n <td class="tabs"> <h3><a href="./gsea',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
		cat(paste(names(enrichmentAnalysis$GSEA.results)[i],' GSEA">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)	
	cat('<td class="tabs"> <h3> Enrichment Summary </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.Sig.adj.pvals.in.both){
			cat(paste('\n <td class="tabs"> <h3><a href="./enrichment',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
			cat(paste(names(enrichmentAnalysis$Sig.adj.pvals.in.both)[i],' Enrichment.summary">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
			}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)
	cat(paste('\n  <td class="tabs"> <h3><a href="./network.html" title="network">Network Analysis</a></h3> </td> </tr> </table>',sep=""), append = TRUE, file = htmlfile)
#Produce the main dataframe: headers and first row	
		t.name<-rownames(enrichmentAnalysis$GSEA.results[[n]])[1]
		while(length(grep(pattern="/",x=t.name)) != 0){
					t.name<-sub(pattern="/",replacement="_",x=t.name,perl=TRUE)
					}
		cat(paste('\n <hr/> \n <br>',names(enrichmentAnalysis$GSEA.results)[n],' GSEA ' ,'\n <br> \n <table class="results"> <tr class="head"> <th>Gene Set name</th> <th>',
			colnames(enrichmentAnalysis$GSEA.results[[n]])[1],"</th> \n <th>",
			colnames(enrichmentAnalysis$GSEA.results[[n]])[2],"</th> \n <th>",
			colnames(enrichmentAnalysis$GSEA.results[[n]])[3],"</th> \n <th>",
			colnames(enrichmentAnalysis$GSEA.results[[n]])[4],"</th> \n <th>",
				"Plots","</th> \n </tr>",
				'<tr class="odd"> <td>',rownames(enrichmentAnalysis$GSEA.results[[n]])[1],"</td> \n <td>",
			signif(enrichmentAnalysis$GSEA.results[[n]][1,1], digits=4),"</td> \n <td>",
			signif(enrichmentAnalysis$GSEA.results[[n]][1,2], digits=4),"</td> \n <td>",
			signif(enrichmentAnalysis$GSEA.results[[n]][1,3], digits=4),"</td> \n <td>",
			signif(enrichmentAnalysis$GSEA.results[[n]][1,4], digits=4),"</td> \n <td>",
				'<a href="../gsea_plots',t.name,'.png" target="_blank" title="gseaplots">plot</a> ',"</td> </tr>", sep=""), append = TRUE, file = htmlfile)
#Produce the main dataframe: rows with plots	
		for(i in 2:nGseaPlots){
			t.name<-rownames(enrichmentAnalysis$GSEA.results[[n]])[i]
			while(length(grep(pattern="/",x=t.name)) != 0){
					t.name<-sub(pattern="/",replacement="_",x=t.name,perl=TRUE)
					}
			if((i%%2) == 0){
				cat(paste('<tr class="even"> <td>',rownames(enrichmentAnalysis$GSEA.results[[n]])[i],"</td> \n <td>",
				signif(enrichmentAnalysis$GSEA.results[[n]][i,1], digits=4),"</td> \n <td>",
				signif(enrichmentAnalysis$GSEA.results[[n]][i,2], digits=4),"</td> \n <td>",
				signif(enrichmentAnalysis$GSEA.results[[n]][i,3], digits=4),"</td> \n <td>",
				signif(enrichmentAnalysis$GSEA.results[[n]][i,4], digits=4),"</td> \n <td>",
					'<a href="../gsea_plots',t.name,'.png" target="_blank" title="gseaplots">plot</a> ',"</td> </tr>", sep=""),append = TRUE, file = htmlfile)
				}else{
					cat(paste('<tr class="odd"> <td>',rownames(enrichmentAnalysis$GSEA.results[[n]])[i],"</td> \n <td>",
					signif(enrichmentAnalysis$GSEA.results[[n]][i,1], digits=4),"</td> \n <td>",
					signif(enrichmentAnalysis$GSEA.results[[n]][i,2], digits=4),"</td> \n <td>",
					signif(enrichmentAnalysis$GSEA.results[[n]][i,3], digits=4),"</td> \n <td>",
					signif(enrichmentAnalysis$GSEA.results[[n]][i,4], digits=4),"</td> \n <td>",
						'<a href="../gsea_plots',t.name,'.png" target="_blank" title="gseaplots">plot</a> ',"</td> </tr>", sep=""),append = TRUE, file = htmlfile)
					}
			}
#Produce the main dataframe: rows without plots			
		for(i in (nGseaPlots+1):dim(enrichmentAnalysis$GSEA.results[[n]])[1]){
			if((i%%2) == 0){
				cat(paste('<tr class="even"> <td>',rownames(enrichmentAnalysis$GSEA.results[[n]])[i],"</td> \n <td>",
				signif(enrichmentAnalysis$GSEA.results[[n]][i,1], digits=4),"</td> \n <td>",
				signif(enrichmentAnalysis$GSEA.results[[n]][i,2], digits=4),"</td> \n <td>",
				signif(enrichmentAnalysis$GSEA.results[[n]][i,3], digits=4),"</td> \n <td>",
				signif(enrichmentAnalysis$GSEA.results[[n]][i,4], digits=4),"</td> \n <td>","</td> \n </tr>", sep=""),append = TRUE, file = htmlfile)
					}else{
					cat(paste('<tr class="odd"> <td>',rownames(enrichmentAnalysis$GSEA.results[[n]])[i],"</td> \n <td>",
					signif(enrichmentAnalysis$GSEA.results[[n]][i,1], digits=4),"</td> \n <td>",
					signif(enrichmentAnalysis$GSEA.results[[n]][i,2], digits=4),"</td> \n <td>",
					signif(enrichmentAnalysis$GSEA.results[[n]][i,3], digits=4),"</td> \n <td>",
					signif(enrichmentAnalysis$GSEA.results[[n]][i,4], digits=4),"</td> \n <td>","</td> \n </tr>", sep=""),append = TRUE, file = htmlfile)
					}
			}	
		cat("</table>",append = TRUE, file = htmlfile)	
		cat('\n </tr> \n </table> \n </td> \n </tr> \n </table> \n </body> \n </html>',append = TRUE, file = htmlfile)
		}
#Second:take care of the Kegg gene sets		
	if(length(intersect(n,keggGS)) != 0){
			cat('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"http://www.w3.org/TR/html4/loose.dtd">', file = htmlfile)
			cat('\n <html> \n <link rel="stylesheet" href="./htsanalyzer.css" type="text/css"> ',append = TRUE, file = htmlfile)
			cat('\n <head> <title> HTSanalyzeR Experiment Report </title> </head>',append = TRUE, file = htmlfile)
			cat('\n <body> \n  <table class="border"> \n <tr class="border top"> \n <td class="border corner"> &nbsp&nbsp&nbsp&nbsp </td>',append = TRUE,file=htmlfile)
			cat(paste('\n <td class="border top"> <div class="header"> Report for Experiment <span class="header"> ',experimentName,'</span> </div>',sep=""),append = TRUE,file=htmlfile)
			cat(paste('<div class="timestamp">generated  ',date(),'(<small>version 1.1.7</small>) </div>',sep=""),append = TRUE,file=htmlfile)
			cat('\n <div class="HTSheader"> HTSanalyzeR </div> </td> </tr> <tr class="border middle"> <td class="border left"></td>',append = TRUE,file=htmlfile)
			cat('<td class="main"> <table> <tr> <div class="HTSlogos"> <img src="./Rlogo.png" width="50" height="40"/>',append = TRUE,file=htmlfile)
			cat('&nbsp <img src="./blue_cruklogo.gif" width="120" height="50"/>',append = TRUE,file=htmlfile)
			cat('&nbsp <img src="./goatcomputer.png" width="85" height="60"/> </div></tr><tr>',append = TRUE,file=htmlfile)
#Produce the tabs
	cat(paste('\n <table class="noframe"> <tr> <td class="tabs"> <h3><a href="./index.html" title="index">Index</a></h3> </td> <tr>',sep=""), append = TRUE, file = htmlfile)
	cat("\n <tr>", append = TRUE, file = htmlfile)
	l.HyperGeo.results<-length(enrichmentAnalysis$HyperGeo.results)
	l.GSEA.results<-length(enrichmentAnalysis$GSEA.results)
	l.Sig.adj.pvals.in.both<-length(enrichmentAnalysis$Sig.adj.pvals.in.both)
	cat('<td class="tabs"> <h3>  </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.HyperGeo.results){
		cat(paste('\n <td class="tabs"> <h3>',names(enrichmentAnalysis$HyperGeo.results)[i],'</h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)	
	cat('<td class="tabs"> <h3> Hypergeometric tests </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.HyperGeo.results){
		cat(paste('\n <td class="tabs"> <h3><a href="./hyperg',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
		cat(paste(names(enrichmentAnalysis$HyperGeo.results)[i],' Hyperg. Tests">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)
	cat('<td class="tabs"> <h3> GSEA </h3> </td>', append = TRUE, file = htmlfile)	
	for(i in 1:l.GSEA.results){
		cat(paste('\n <td class="tabs"> <h3><a href="./gsea',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
		cat(paste(names(enrichmentAnalysis$GSEA.results)[i],' GSEA">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)	
	cat('<td class="tabs"> <h3> Enrichment Summary </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.Sig.adj.pvals.in.both){
			cat(paste('\n <td class="tabs"> <h3><a href="./enrichment',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
			cat(paste(names(enrichmentAnalysis$Sig.adj.pvals.in.both)[i],' Enrichment.summary">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
			}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)
	cat(paste('\n  <td class="tabs"> <h3><a href="./network.html" title="network">Network Analysis</a></h3> </td> </tr> </table>',sep=""), append = TRUE, file = htmlfile)
#Produce the main dataframe: headers and first row	
		kegggsname<-sub(pattern="(\\d*$)",replacement="",x=rownames(enrichmentAnalysis$GSEA.results[[n]])[1],perl=TRUE)
		kegggsname<-sub(pattern="(\\D*$)",replacement="",x=kegggsname,perl=TRUE)
		kegggsname<-sub(pattern="( :: \\D*$)",replacement="",x=kegggsname,perl=TRUE)
 		kegggsname<-sub(pattern="( :: \\w*$)",replacement="",x=kegggsname,perl=TRUE)
		cat(paste('\n <hr/> \n <br>',names(enrichmentAnalysis$GSEA.results)[n],' GSEA ' ,'\n <br> \n <table class="results"> <tr class="head"> <th>Gene Set name</th> <th>',
			colnames(enrichmentAnalysis$GSEA.results[[n]])[1],"</th> \n <th>",
			colnames(enrichmentAnalysis$GSEA.results[[n]])[2],"</th> \n <th>",
			colnames(enrichmentAnalysis$GSEA.results[[n]])[3],"</th> \n <th>",
			colnames(enrichmentAnalysis$GSEA.results[[n]])[4],"</th> \n <th>",
				"Plots","</th> \n </tr>",
				'<tr class="odd"> <td><a href="http://www.genome.jp/dbget-bin/www_bget?pathway:',kegggsname,'" target="_blank" title="',rownames(enrichmentAnalysis$GSEA.results[[n]])[1],'">',rownames(enrichmentAnalysis$GSEA.results[[n]])[1],'</a> </td> \n <td>',
			signif(enrichmentAnalysis$GSEA.results[[n]][1,1], digits=4),"</td> \n <td>",
			signif(enrichmentAnalysis$GSEA.results[[n]][1,2], digits=4),"</td> \n <td>",
			signif(enrichmentAnalysis$GSEA.results[[n]][1,3], digits=4),"</td> \n <td>",
			signif(enrichmentAnalysis$GSEA.results[[n]][1,4], digits=4),"</td> \n <td>",
				'<a href="../gsea_plots',kegggsname,'.png" target="_blank" title="gseaplots">plot</a> ',"</td> </tr>", sep=""), append = TRUE, file = htmlfile)
#Produce the main dataframe: rows with plots	
		for(i in 2:nGseaPlots){
			kegggsname<-sub(pattern="(\\D*$)",replacement="",x=rownames(enrichmentAnalysis$GSEA.results[[n]])[i],perl=TRUE)
			kegggsname<-sub(pattern="(\\D*$)",replacement="",x=kegggsname,perl=TRUE)
			kegggsname<-sub(pattern="( :: \\D*$)",replacement="",x=kegggsname,perl=TRUE)
			kegggsname<-sub(pattern="( :: \\w*$)",replacement="",x=kegggsname,perl=TRUE)
			if((i%%2) == 0){
				cat(paste('<tr class="even"> <td><a href="http://www.genome.jp/dbget-bin/www_bget?pathway:',kegggsname,'" target="_blank" title="',rownames(enrichmentAnalysis$GSEA.results[[n]])[i],'">',rownames(enrichmentAnalysis$GSEA.results[[n]])[i],'</a> </td> \n <td>',
				signif(enrichmentAnalysis$GSEA.results[[n]][i,1], digits=4),"</td> \n <td>",
				signif(enrichmentAnalysis$GSEA.results[[n]][i,2], digits=4),"</td> \n <td>",
				signif(enrichmentAnalysis$GSEA.results[[n]][i,3], digits=4),"</td> \n <td>",
				signif(enrichmentAnalysis$GSEA.results[[n]][i,4], digits=4),"</td> \n <td>",
					'<a href="../gsea_plots',kegggsname,'.png" target="_blank" title="gseaplots">plot</a> ',"</td> </tr>", sep=""),append = TRUE, file = htmlfile)
				}else{
					cat(paste('<tr class="odd"> <td><a href="http://www.genome.jp/dbget-bin/www_bget?pathway:',kegggsname,'" target="_blank" title="',rownames(enrichmentAnalysis$GSEA.results[[n]])[i],'">',rownames(enrichmentAnalysis$GSEA.results[[n]])[i],'</a> </td> \n <td>',
					signif(enrichmentAnalysis$GSEA.results[[n]][i,1], digits=4),"</td> \n <td>",
					signif(enrichmentAnalysis$GSEA.results[[n]][i,2], digits=4),"</td> \n <td>",
					signif(enrichmentAnalysis$GSEA.results[[n]][i,3], digits=4),"</td> \n <td>",
					signif(enrichmentAnalysis$GSEA.results[[n]][i,4], digits=4),"</td> \n <td>",
						'<a href="../gsea_plots',kegggsname,'.png" target="_blank" title="gseaplots">plot</a> ',"</td> </tr>", sep=""),append = TRUE, file = htmlfile)
					}
			}
#Produce the main dataframe: rows without plot				
		for(i in (nGseaPlots+1):dim(enrichmentAnalysis$GSEA.results[[n]])[1]){
			kegggsname<-sub(pattern="(\\D*$)",replacement="",x=rownames(enrichmentAnalysis$GSEA.results[[n]])[i],perl=TRUE)
			kegggsname<-sub(pattern="(\\D*$)",replacement="",x=kegggsname,perl=TRUE)
			kegggsname<-sub(pattern="( :: \\D*$)",replacement="",x=kegggsname,perl=TRUE)
			kegggsname<-sub(pattern="( :: \\w*$)",replacement="",x=kegggsname,perl=TRUE)
			if((i%%2) == 0){
				cat(paste('<tr class="even"> <td><a href="http://www.genome.jp/dbget-bin/www_bget?pathway:',kegggsname,'" target="_blank" title="',rownames(enrichmentAnalysis$GSEA.results[[n]])[i],'">',rownames(enrichmentAnalysis$GSEA.results[[n]])[i],'</a> </td> \n <td>',
				signif(enrichmentAnalysis$GSEA.results[[n]][i,1], digits=4),"</td> \n <td>",
				signif(enrichmentAnalysis$GSEA.results[[n]][i,2], digits=4),"</td> \n <td>",
				signif(enrichmentAnalysis$GSEA.results[[n]][i,3], digits=4),"</td> \n <td>",
				signif(enrichmentAnalysis$GSEA.results[[n]][i,4], digits=4),"</td> \n <td>","</td> \n </tr>", sep=""),append = TRUE, file = htmlfile)
					}else{
					cat(paste('<tr class="odd"> <td><a href="http://www.genome.jp/dbget-bin/www_bget?pathway:',kegggsname,'" target="_blank" title="',rownames(enrichmentAnalysis$GSEA.results[[n]])[i],'">',rownames(enrichmentAnalysis$GSEA.results[[n]])[i],'</a> </td> \n <td>',
					signif(enrichmentAnalysis$GSEA.results[[n]][i,1], digits=4),"</td> \n <td>",
					signif(enrichmentAnalysis$GSEA.results[[n]][i,2], digits=4),"</td> \n <td>",
					signif(enrichmentAnalysis$GSEA.results[[n]][i,3], digits=4),"</td> \n <td>",
					signif(enrichmentAnalysis$GSEA.results[[n]][i,4], digits=4),"</td> \n <td>","</td> \n </tr>", sep=""),append = TRUE, file = htmlfile)
					}
			}	
		cat("</table>",append = TRUE, file = htmlfile)	
		cat('\n </tr> \n </table> \n </td> \n </tr> \n </table> \n </body> \n </html>',append = TRUE, file = htmlfile)
		}
#Third, take care of the GO gene sets		
	if(length(intersect(n,GOGS)) != 0){
			cat('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"http://www.w3.org/TR/html4/loose.dtd">', file = htmlfile)
			cat('\n <html> \n <link rel="stylesheet" href="./htsanalyzer.css" type="text/css"> ',append = TRUE, file = htmlfile)
			cat('\n <head> <title> HTSanalyzeR Experiment Report </title> </head>',append = TRUE, file = htmlfile)
			cat('\n <body> \n  <table class="border"> \n <tr class="border top"> \n <td class="border corner"> &nbsp&nbsp&nbsp&nbsp </td>',append = TRUE,file=htmlfile)
			cat(paste('\n <td class="border top"> <div class="header"> Report for Experiment <span class="header"> ',experimentName,'</span> </div>',sep=""),append = TRUE,file=htmlfile)
			cat(paste('<div class="timestamp">generated  ',date(),'(<small>version 1.1.7</small>) </div>',sep=""),append = TRUE,file=htmlfile)
			cat('\n <div class="HTSheader"> HTSanalyzeR </div> </td> </tr> <tr class="border middle"> <td class="border left"></td>',append = TRUE,file=htmlfile)
			cat('<td class="main"> <table> <tr> <div class="HTSlogos"> <img src="./Rlogo.png" width="50" height="40"/>',append = TRUE,file=htmlfile)
			cat('&nbsp <img src="./blue_cruklogo.gif" width="120" height="50"/>',append = TRUE,file=htmlfile)
			cat('&nbsp <img src="./goatcomputer.png" width="85" height="60"/> </div></tr><tr>',append = TRUE,file=htmlfile)		
#Produce the tabs
	cat(paste('\n <table class="noframe"> <tr> <td class="tabs"> <h3><a href="./index.html" title="index">Index</a></h3> </td> <tr>',sep=""), append = TRUE, file = htmlfile)
	cat("\n <tr>", append = TRUE, file = htmlfile)
	l.HyperGeo.results<-length(enrichmentAnalysis$HyperGeo.results)
	l.GSEA.results<-length(enrichmentAnalysis$GSEA.results)
	l.Sig.adj.pvals.in.both<-length(enrichmentAnalysis$Sig.adj.pvals.in.both)
	cat('<td class="tabs"> <h3>  </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.HyperGeo.results){
		cat(paste('\n <td class="tabs"> <h3>',names(enrichmentAnalysis$HyperGeo.results)[i],'</h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)	
	cat('<td class="tabs"> <h3> Hypergeometric tests </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.HyperGeo.results){
		cat(paste('\n <td class="tabs"> <h3><a href="./hyperg',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
		cat(paste(names(enrichmentAnalysis$HyperGeo.results)[i],' Hyperg. Tests">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)
	cat('<td class="tabs"> <h3> GSEA </h3> </td>', append = TRUE, file = htmlfile)	
	for(i in 1:l.GSEA.results){
		cat(paste('\n <td class="tabs"> <h3><a href="./gsea',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
		cat(paste(names(enrichmentAnalysis$GSEA.results)[i],' GSEA">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)	
	cat('<td class="tabs"> <h3> Enrichment Summary </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.Sig.adj.pvals.in.both){
			cat(paste('\n <td class="tabs"> <h3><a href="./enrichment',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
			cat(paste(names(enrichmentAnalysis$Sig.adj.pvals.in.both)[i],' Enrichment.summary">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
			}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)
	cat(paste('\n  <td class="tabs"> <h3><a href="./network.html" title="network">Network Analysis</a></h3> </td> </tr> </table>',sep=""), append = TRUE, file = htmlfile)
##Produce the main dataframe: headers and first row	
		cat(paste('\n <hr/> \n <br>',names(enrichmentAnalysis$GSEA.results)[n],' GSEA ' ,'\n <br> \n <table class="results"> <tr class="head"> <th>Gene Set name</th> <th>',
			colnames(enrichmentAnalysis$GSEA.results[[n]])[1],"</th> \n <th>",
			colnames(enrichmentAnalysis$GSEA.results[[n]])[2],"</th> \n <th>",
			colnames(enrichmentAnalysis$GSEA.results[[n]])[3],"</th> \n <th>",
			colnames(enrichmentAnalysis$GSEA.results[[n]])[4],"</th> \n <th>",
				"Plots","</th> \n </tr>",
				'<tr class="odd"> <td><a href="http://www.ebi.ac.uk/QuickGO/GTerm?id=GO:',sub(pattern="(\\D*$)",replacement="",x=sub(pattern="(\\D*)",replacement="",x=rownames(enrichmentAnalysis$GSEA.results[[n]])[1],perl=TRUE),perl=TRUE),'" target="_blank" title="',rownames(enrichmentAnalysis$GSEA.results[[n]])[1],'">',rownames(enrichmentAnalysis$GSEA.results[[n]])[1],'</a> </td> \n <td>',
			signif(enrichmentAnalysis$GSEA.results[[n]][1,1], digits=4),"</td> \n <td>",
			signif(enrichmentAnalysis$GSEA.results[[n]][1,2], digits=4),"</td> \n <td>",
			signif(enrichmentAnalysis$GSEA.results[[n]][1,3], digits=4),"</td> \n <td>",
			signif(enrichmentAnalysis$GSEA.results[[n]][1,4], digits=4),"</td> \n <td>",
				'<a href="../gsea_plots',sub(pattern="(\\D*$)",replacement="",x=rownames(enrichmentAnalysis$GSEA.results[[n]])[1],perl=TRUE),'.png" target="_blank" title="gseaplots">plot</a> ',"</td> </tr>", sep=""), append = TRUE, file = htmlfile)
#Produce the main dataframe: rows with a plot	
		for(i in 2:nGseaPlots){
			if((i%%2) == 0){
				cat(paste('<tr class="even"> <td><a href="http://www.ebi.ac.uk/QuickGO/GTerm?id=GO:',sub(pattern="(\\D*$)",replacement="",x=sub(pattern="(\\D*)",replacement="",x=rownames(enrichmentAnalysis$GSEA.results[[n]])[i],perl=TRUE),perl=TRUE),'" target="_blank" title="',rownames(enrichmentAnalysis$GSEA.results[[n]])[i],'">',rownames(enrichmentAnalysis$GSEA.results[[n]])[i],'</a> </td> \n <td>',
				signif(enrichmentAnalysis$GSEA.results[[n]][i,1], digits=4),"</td> \n <td>",
				signif(enrichmentAnalysis$GSEA.results[[n]][i,2], digits=4),"</td> \n <td>",
				signif(enrichmentAnalysis$GSEA.results[[n]][i,3], digits=4),"</td> \n <td>",
				signif(enrichmentAnalysis$GSEA.results[[n]][i,4], digits=4),"</td> \n <td>",
					'<a href="../gsea_plots',sub(pattern="(\\D*$)",replacement="",x=rownames(enrichmentAnalysis$GSEA.results[[n]])[i],perl=TRUE),'.png" target="_blank" title="gseaplots">plot</a> ',"</td> </tr>", sep=""),append = TRUE, file = htmlfile)
				}else{
					cat(paste('<tr class="odd"> <td><a href="http://www.ebi.ac.uk/QuickGO/GTerm?id=GO:',sub(pattern="(\\D*$)",replacement="",x=sub(pattern="(\\D*)",replacement="",x=rownames(enrichmentAnalysis$GSEA.results[[n]])[i],perl=TRUE),perl=TRUE),'" target="_blank" title="',rownames(enrichmentAnalysis$GSEA.results[[n]])[i],'">',rownames(enrichmentAnalysis$GSEA.results[[n]])[i],'</a> </td> \n <td>',
					signif(enrichmentAnalysis$GSEA.results[[n]][i,1], digits=4),"</td> \n <td>",
					signif(enrichmentAnalysis$GSEA.results[[n]][i,2], digits=4),"</td> \n <td>",
					signif(enrichmentAnalysis$GSEA.results[[n]][i,3], digits=4),"</td> \n <td>",
					signif(enrichmentAnalysis$GSEA.results[[n]][i,4], digits=4),"</td> \n <td>",
						'<a href="../gsea_plots',sub(pattern="(\\D*$)",replacement="",x=rownames(enrichmentAnalysis$GSEA.results[[n]])[i],perl=TRUE),'.png" target="_blank" title="gseaplots">plot</a> ',"</td> </tr>", sep=""),append = TRUE, file = htmlfile)
					}
			}
#Produce the main dataframe: rows without a plot
		for(i in (nGseaPlots+1):dim(enrichmentAnalysis$GSEA.results[[n]])[1]){
			if((i%%2) == 0){
				cat(paste('<tr class="even"> <td><a href="http://www.ebi.ac.uk/QuickGO/GTerm?id=GO:',sub(pattern="(\\D*$)",replacement="",x=sub(pattern="(\\D*)",replacement="",x=rownames(enrichmentAnalysis$GSEA.results[[n]])[i],perl=TRUE),perl=TRUE),'" target="_blank" title="',rownames(enrichmentAnalysis$GSEA.results[[n]])[i],'">',rownames(enrichmentAnalysis$GSEA.results[[n]])[i],'</a> </td> \n <td>',
				signif(enrichmentAnalysis$GSEA.results[[n]][i,1], digits=4),"</td> \n <td>",
				signif(enrichmentAnalysis$GSEA.results[[n]][i,2], digits=4),"</td> \n <td>",
				signif(enrichmentAnalysis$GSEA.results[[n]][i,3], digits=4),"</td> \n <td>",
				signif(enrichmentAnalysis$GSEA.results[[n]][i,4], digits=4),"</td> \n <td>","</td> \n </tr>", sep=""),append = TRUE, file = htmlfile)
					}else{
					cat(paste('<tr class="odd"> <td><a href="http://www.ebi.ac.uk/QuickGO/GTerm?id=GO:',sub(pattern="(\\D*$)",replacement="",x=sub(pattern="(\\D*)",replacement="",x=rownames(enrichmentAnalysis$GSEA.results[[n]])[i],perl=TRUE),perl=TRUE),'" target="_blank" title="',rownames(enrichmentAnalysis$GSEA.results[[n]])[i],'">',rownames(enrichmentAnalysis$GSEA.results[[n]])[i],'</a> </td> \n <td>',
					signif(enrichmentAnalysis$GSEA.results[[n]][i,1], digits=4),"</td> \n <td>",
					signif(enrichmentAnalysis$GSEA.results[[n]][i,2], digits=4),"</td> \n <td>",
					signif(enrichmentAnalysis$GSEA.results[[n]][i,3], digits=4),"</td> \n <td>",
					signif(enrichmentAnalysis$GSEA.results[[n]][i,4], digits=4),"</td> \n <td>","</td> \n </tr>", sep=""),append = TRUE, file = htmlfile)
					}
			}	
		cat("</table>",append = TRUE, file = htmlfile)	
		cat('\n </tr> \n </table> \n </td> \n </tr> \n </table> \n </body> \n </html>',append = TRUE, file = htmlfile)
		}	
	}	
#Produce the enrichment summary page
for(n in 1:l.Sig.adj.pvals.in.both){
	htmlfile =  file.path(reportdir,paste("enrichment",n,".html",sep=""))
	cat('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"http://www.w3.org/TR/html4/loose.dtd">', file = htmlfile)
	cat('\n <html> \n <link rel="stylesheet" href="./htsanalyzer.css" type="text/css"> ',append = TRUE, file = htmlfile)
	cat('\n <head> <title> HTSanalyzeR Experiment Report </title> </head>',append = TRUE, file = htmlfile)
	cat('\n <body> \n  <table class="border"> \n <tr class="border top"> \n <td class="border corner"> &nbsp&nbsp&nbsp&nbsp </td>',append = TRUE,file=htmlfile)
	cat(paste('\n <td class="border top"> <div class="header"> Report for Experiment <span class="header"> ',experimentName,'</span> </div>',sep=""),append = TRUE,file=htmlfile)
	cat(paste('<div class="timestamp">generated  ',date(),'(<small>version 1.1.7</small>) </div>',sep=""),append = TRUE,file=htmlfile)
	cat('\n <div class="HTSheader"> HTSanalyzeR </div> </td> </tr> <tr class="border middle"> <td class="border left"></td>',append = TRUE,file=htmlfile)
	cat('<td class="main"> <table> <tr> <div class="HTSlogos"> <img src="./Rlogo.png" width="50" height="40"/>',append = TRUE,file=htmlfile)
	cat('&nbsp <img src="./blue_cruklogo.gif" width="120" height="50"/>',append = TRUE,file=htmlfile)
	cat('&nbsp <img src="./goatcomputer.png" width="85" height="60"/> </div></tr><tr>',append = TRUE,file=htmlfile)
#Produce the tabs
	cat(paste('\n <table class="noframe"> <tr> <td class="tabs"> <h3><a href="./index.html" title="index">Index</a></h3> </td> <tr>',sep=""), append = TRUE, file = htmlfile)
	cat("\n <tr>", append = TRUE, file = htmlfile)
	cat('<td class="tabs"> <h3>  </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.HyperGeo.results){
		cat(paste('\n <td class="tabs"> <h3>',names(enrichmentAnalysis$HyperGeo.results)[i],'</h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)	
	cat('<td class="tabs"> <h3> Hypergeometric tests </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.HyperGeo.results){
		cat(paste('\n <td class="tabs"> <h3><a href="./hyperg',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
		cat(paste(names(enrichmentAnalysis$HyperGeo.results)[i],' Hyperg. Tests">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)
	cat('<td class="tabs"> <h3> GSEA </h3> </td>', append = TRUE, file = htmlfile)	
	for(i in 1:l.GSEA.results){
		cat(paste('\n <td class="tabs"> <h3><a href="./gsea',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
		cat(paste(names(enrichmentAnalysis$GSEA.results)[i],' GSEA">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)	
	cat('<td class="tabs"> <h3> Enrichment Summary </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.Sig.adj.pvals.in.both){
			cat(paste('\n <td class="tabs"> <h3><a href="./enrichment',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
			cat(paste(names(enrichmentAnalysis$Sig.adj.pvals.in.both)[i],' Enrichment.summary">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
			}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)
	cat(paste('\n  <td class="tabs"> <h3><a href="./network.html" title="network">Network Analysis</a></h3> </td> </tr> </table>',sep=""), append = TRUE, file = htmlfile)
#Produce the main dataframe	
	if(dim(enrichmentAnalysis$Sig.adj.pvals.in.both[[n]])[1] !=0){
		coll<-names(enrichmentAnalysis$Sig.adj.pvals.in.both)[n]
		cat(paste('\n <hr/> \n <br>', 'Gene sets with significant adjusted p value with both hypergeometric test and GSEA: ', coll ,'<br>','\n <table class="results"> <tr class="head"> <th>Gene Set name </th> \n <th>',
			colnames(enrichmentAnalysis$Sig.adj.pvals.in.both[[n]])[1],"</th> \n <th>",
			colnames(enrichmentAnalysis$Sig.adj.pvals.in.both[[n]])[2],"</th> \n </tr>",
			'<tr class="odd"> <td>',rownames(enrichmentAnalysis$Sig.adj.pvals.in.both[[n]])[1],"</td> \n <td>",
			signif(enrichmentAnalysis$Sig.adj.pvals.in.both[[n]][1,1], digits=4),"</td> \n <td>",
			signif(enrichmentAnalysis$Sig.adj.pvals.in.both[[n]][1,2], digits=4),"</td> \n </tr>", sep=""), append = TRUE, file = htmlfile)
		if(dim(enrichmentAnalysis$Sig.adj.pvals.in.both[[n]])[1] > 1){
			for(i in 2:dim(enrichmentAnalysis$Sig.adj.pvals.in.both[[n]])[1]){
				if((i%%2) == 0){
					cat(paste('<tr class="even"> <td>',rownames(enrichmentAnalysis$Sig.adj.pvals.in.both[[n]])[i],"</td> \n <td>",
					signif(enrichmentAnalysis$Sig.adj.pvals.in.both[[n]][i,1], digits=4),"</td> \n <td>",
					signif(enrichmentAnalysis$Sig.adj.pvals.in.both[[n]][i,2], digits=4),"</td> </tr>", sep=""),append = TRUE, file = htmlfile)
					}else{
						cat(paste('<tr class="odd"> <td>',rownames(enrichmentAnalysis$Sig.adj.pvals.in.both[[n]])[i],"</td> \n <td>",
						signif(enrichmentAnalysis$Sig.adj.pvals.in.both[[n]][i,1], digits=4),"</td> \n <td>",
						signif(enrichmentAnalysis$Sig.adj.pvals.in.both[[n]][i,2], digits=4),"</td> </tr>", sep=""),append = TRUE, file = htmlfile)
						}	
				}
			}	
		cat("</table>",append = TRUE, file = htmlfile)		
		cat('\n </tr> \n </table> \n </td> \n </tr> \n </table> \n </body> \n </html>',append = TRUE, file = htmlfile)
		}else{
			coll<-names(enrichmentAnalysis$Sig.adj.pvals.in.both)[n]
			cat(paste('\n <hr/> \n <br>', 'Gene sets with significant adjusted p value with both hypergeometric test and GSEA: ', coll ,'<br>',sep=""), append = TRUE, file = htmlfile)
			cat("There are no gene sets corresponding to these criteria \n </body> \n </html>", append = TRUE, file = htmlfile)
			}	
	}
#Produce the networks analysis results page
	htmlfile =  file.path(reportdir,"network.html")
	cat('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"http://www.w3.org/TR/html4/loose.dtd">', file = htmlfile)
	cat('\n <html> \n <link rel="stylesheet" href="./htsanalyzer.css" type="text/css"> ',append = TRUE, file = htmlfile)
	cat('\n <head> <title> HTSanalyzeR Experiment Report </title> </head>',append = TRUE, file = htmlfile)
	cat('\n <body> \n  <table class="border"> \n <tr class="border top"> \n <td class="border corner"> &nbsp&nbsp&nbsp&nbsp </td>',append = TRUE,file=htmlfile)
	cat(paste('\n <td class="border top"> <div class="header"> Report for Experiment <span class="header"> ',experimentName,'</span> </div>',sep=""),append = TRUE,file=htmlfile)
	cat(paste('<div class="timestamp">generated  ',date(),'(<small>version 1.1.7</small>) </div>',sep=""),append = TRUE,file=htmlfile)
	cat('\n <div class="HTSheader"> HTSanalyzeR </div> </td> </tr> <tr class="border middle"> <td class="border left"></td>',append = TRUE,file=htmlfile)
	cat('<td class="main"> <table> <tr> <div class="HTSlogos"> <img src="./Rlogo.png" width="50" height="40"/>',append = TRUE,file=htmlfile)
	cat('&nbsp <img src="./blue_cruklogo.gif" width="120" height="50"/>',append = TRUE,file=htmlfile)
	cat('&nbsp <img src="./goatcomputer.png" width="85" height="60"/> </div></tr><tr>',append = TRUE,file=htmlfile)
#Produce the tabs
	cat(paste('\n <table class="noframe"> <tr> <td class="tabs"> <h3><a href="./index.html" title="index">Index</a></h3> </td> <tr>',sep=""), append = TRUE, file = htmlfile)
	cat("\n <tr>", append = TRUE, file = htmlfile)
	l.HyperGeo.results<-length(enrichmentAnalysis$HyperGeo.results)
	l.GSEA.results<-length(enrichmentAnalysis$GSEA.results)
	l.Sig.adj.pvals.in.both<-length(enrichmentAnalysis$Sig.adj.pvals.in.both)
	cat('<td class="tabs"> <h3>  </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.HyperGeo.results){
		cat(paste('\n <td class="tabs"> <h3>',names(enrichmentAnalysis$HyperGeo.results)[i],'</h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)	
	cat('<td class="tabs"> <h3> Hypergeometric tests </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.HyperGeo.results){
		cat(paste('\n <td class="tabs"> <h3><a href="./hyperg',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
		cat(paste(names(enrichmentAnalysis$HyperGeo.results)[i],' Hyperg. Tests">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)
	cat('<td class="tabs"> <h3> GSEA </h3> </td>', append = TRUE, file = htmlfile)	
	for(i in 1:l.GSEA.results){
		cat(paste('\n <td class="tabs"> <h3><a href="./gsea',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
		cat(paste(names(enrichmentAnalysis$GSEA.results)[i],' GSEA">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)	
	cat('<td class="tabs"> <h3> Enrichment Summary </h3> </td>', append = TRUE, file = htmlfile)
	for(i in 1:l.Sig.adj.pvals.in.both){
			cat(paste('\n <td class="tabs"> <h3><a href="./enrichment',i,'.html" title="',sep=""), append = TRUE, file = htmlfile)
			cat(paste(names(enrichmentAnalysis$Sig.adj.pvals.in.both)[i],' Enrichment.summary">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
			}
	cat("\n </tr> <tr>", append = TRUE, file = htmlfile)
	cat(paste('\n  <td class="tabs"> <h3><a href="./network.html" title="network">Network Analysis</a></h3> </td> </tr> </table>',sep=""), append = TRUE, file = htmlfile)
#Produce a link to the graph file containing the module	
	cat(paste('\n <hr/> \n  <br>Click <a href="../',nwAnalysisGraphFile,'" target="_blank" title="Enriched Subnetwork">here</a> to get the enriched subnetwork ',sep=""), append = TRUE, file = htmlfile)
#Check that the nwAnalysisOutput has the right format
	if(is.list(nwAnalysisOutput) == FALSE) stop("The nwAnalysisOutput should be a list")
	if(length(intersect(names(nwAnalysisOutput),c("subnw","labels"))) != 2) stop("The nwAnalysisOutput should contain the following elements: 'subnw', 'labels'")
	if(class(nwAnalysisOutput$subnw) != "graphNEL") stop("The nwAnalysisOutput$subnw should be a graphNEL object")
	if(length(nwAnalysisOutput$labels) != length(nodes(nwAnalysisOutput$subnw))) stop("The nwAnalysisOutput$label should be a vector of the length nodes(nwAnalysisOutput$subnw)")
	EnrichSNnodes<-cbind(nodes(nwAnalysisOutput$subnw),nwAnalysisOutput$labels)
	colnames(EnrichSNnodes)<-c("Entrez Identifier","Symbol")
	cat(paste('\n <table class="noframe"> <tr> <td> <table class="results"> <tr class="head"> \n <th>',
		colnames(EnrichSNnodes)[1],"</th> \n <th>",
		colnames(EnrichSNnodes)[2],"</th> ",
		'<tr class="odd"> <td>',EnrichSNnodes[1,1],"</td> \n <td>",
		EnrichSNnodes[1,2],"</td> \n </tr>", sep=""), append = TRUE, file = htmlfile)
	for(i in 2:dim(EnrichSNnodes)[1]){
	if((i%%2) == 0){
		cat(paste('<tr class="even"> <td>',EnrichSNnodes[i,1],"</td> \n <td>",
		EnrichSNnodes[i,2],"</td> </tr>", sep=""),append = TRUE, file = htmlfile)
		}else{
		cat(paste('<tr class="odd"> <td>',EnrichSNnodes[i,1],"</td> \n <td>",
		EnrichSNnodes[i,2],"</td> </tr>", sep=""),append = TRUE, file = htmlfile)
		}	
	}
	cat(paste('</table> </td> <td> <img src="../',nwAnalysisGraphFile,'" align="top" width="800" height="800"> </td> </tr> </table>',sep=""),append = TRUE, file = htmlfile)	
	cat('\n </tr> \n </table> \n </td> \n </tr> \n </table> \n </body> \n </html>',append = TRUE, file = htmlfile)
}

