#This function writes an html report following a complete analysis of a dataset 
#with the two wrapper functions geneSetAnalyze and networkAnalysis. 
#Most of the parameters of this function can simply be copied and paste 
#from the parameters of the two above functions.
makeGSEAplots <- function(geneList, geneSet, exponent,
		filepath, filename,
		output='png'){
	test <- gseaScores(geneList=geneList, geneSet=geneSet,
			exponent=exponent, mode="graph")
	gseaPlots(Running.Score=test[['Running.Score']],
			Enrichment.Score=test[['Enrichment.Score']],
			Positions=test[['Positions']], geneList=geneList,
			output=output, filepath=filepath, filename=filename)
}

makeOverlapTable <- function(geneSet, hits, mapID, filepath, filename){
	overlap <- intersect(geneSet, hits)
	nOverlap <- length(overlap)
	overlapSymbols <- rep(0, nOverlap)
	if (nOverlap > 0)
		overlapSymbols <- sapply(mapID[overlap], function(x) ifelse(length(x) == 1, x, 0))
	filename <- file.path(filepath, paste(filename, ".txt", sep=""))
	overlap <- cbind(EntrezGene=overlap, Symbols=overlapSymbols)
	write.table(overlap, file=filename, row.names=FALSE, quote=FALSE)
}


writeReportHTSA <- function(experimentName, enrichmentAnalysis,
		cutoffHits, hits, listOfGeneSetCollections,
		geneList, geneListName, p.adj.method="BH",
		nPermutations=1000, min.gene.set.size=15,
		exponent=1, nwAnalysisOutput,
		controls="neg", alternative="two.sided",
		tests="T-test",
		columns=c("t.test.pvalues.two.samples","t.test.pvalues.one.sample"),
		species="Dm", fdr=0.001, genetic=FALSE,
		networkObject=NA, nGseaPlots=10,
		whichSetIsKEGGIds, whichSetIsGOIds,
		reportdir="HTSanalyzerReport"){
	##########################################
	#       		check input 		     #
	##########################################
	## whichSetIsKEGGIds: numeric vector
	## whichSetIsGOIds: numeric vector
	## enrichmentAnalysis: list and containing (HyperGeo.results,
	##                       GSEA.results, Sig.adj.pvals.in.both)
	## nGseaPlots: numeric (length 1, positive and <= than the smaller
	##             geneset)
	## listOfGeneSetCollections: named list
	
	## Get the indices of the gene set collections for which there will be
	## links to the relevant databases
	if (missing(whichSetIsKEGGIds))
		whichSetIsKEGGIds <- 0
	if (missing(whichSetIsGOIds))
		whichSetIsGOIds <- 0
	
	## Error checking: add helpers here for future S4
	species <- match.arg(species, c("Dm", "Hs", "Rn", "Mm", "Ce"))
	stopifnot(is.numeric(whichSetIsKEGGIds),
			is.numeric(whichSetIsGOIds),
			is.list(listOfGeneSetCollections),
			is.character(names(listOfGeneSetCollections)),
			is.list(enrichmentAnalysis),
			all(c("HyperGeo.results","GSEA.results","Sig.adj.pvals.in.both")
							%in% names(enrichmentAnalysis)),
			length(nGseaPlots) == 1, is.numeric(nGseaPlots),
			nGseaPlots > 0,
			nGseaPlots <= min(sapply(enrichmentAnalysis[['GSEA.results']], nrow)))
	##########################################
	#           create directories 		     #
	##########################################
	docdir<-file.path(reportdir,"doc")				#dir for documents or text files
	imgdir<-file.path(reportdir,"image")			#dir for image files
	htmdir<-file.path(reportdir,"html")				#dir for html files
	for(htsdir in c(reportdir,docdir,imgdir,htmdir)) {
		if(!file.exists(htsdir)) 
			dir.create(htsdir)
	}
	##########################################
	#       	produce GSEA plots		     #
	##########################################
	nGSC <- length(listOfGeneSetCollections)
	mapID <- as.list(get(paste("org", species, "egSYMBOL", sep=".")))
	
	for(s in 1:nGSC){
		gsc <- listOfGeneSetCollections[[s]]
		nmg <- names(listOfGeneSetCollections)[s]
		gseaR <- rownames(enrichmentAnalysis[['GSEA.results']][[nmg]])
		hgeoR <- rownames(enrichmentAnalysis[['HyperGeo.results']][[nmg]])
		
		## 1. for the non-GO and non-KEGG gene set collections
		if (!(s %in% union(whichSetIsKEGGIds, whichSetIsGOIds)))
			for(p in 1:nGseaPlots){
				makeGSEAplots(geneList=geneList, geneSet=gsc[[gseaR[p]]],
						exponent=exponent, filepath=dirs['image'],
						filename=gsub("/", "_", gseaR[p]))
				hgeo <- gsub("/", "_", hgeoR[p])
				makeOverlapTable(geneSet=gsc[[hgeoR[p]]], hits=hits,
						mapID=mapID, filepath=dirs['doc'],
						filename=hgeo)
			}
		
		## 2. for the GO gene set collections
		if (s %in% whichSetIsGOIds)
			for(p in 1:nGseaPlots){
				GOgsname <- gsub("[^0-9]*$", "", gseaR[p])
				makeGSEAplots(geneList=geneList, exponent=exponent,
						output='png', geneSet=gsc[[GOgsname]],
						filename=GOgsname, filepath=dirs['image'])
				GOgsname <- sub(pattern="(\\D*$)", replacement="",
						x=hgeoR[p], perl=TRUE)
				makeOverlapTable(geneSet=gsc[[GOgsname]], hits=hits,
						mapID=mapID, filepath=dirs['doc'],
						filename=GOgsname)
			}
		
		## 3.for the KEGG gene set collections
		if (s %in% whichSetIsKEGGIds)
			for(p in 1:nGseaPlots){
				makeGSEAplots(geneList=geneList, exponent=exponent, output='png',
						filepath=dirs['image'], filename=gseaR[p], geneSet=gsc[[gseaR[p]]])
				makeOverlapTable(geneSet=gsc[[hgeoR[p]]], hits=hits,
						mapID=mapID, filepath=dirs['doc'],
						filename=hgeoR[p])
			}
		rm(gsc, nmg, gseaR, hgeoR)
	}
	#The gseaScores function will check:
		#that the geneList has the right format
		#that the gene list has been named properly
		#that the exponent is a single integer
		#that the geneSet is a single vector	
	##########################################
	#       	produce html templates	     #
	##########################################
	#Copy css and logos in there
	cpfile<-dir(system.file("templates",package="HTSanalyzeR"),full=TRUE)
	file.copy(from=cpfile,to=imgdir,overwrite=TRUE)
	##########################################
	#       	produce index.html		     #
	##########################################
	#Produce the index html page
	htmlfile<-file.path(reportdir,"index.html")
	writeHTSAHtmlHead(experimentName=experimentName, htmlfile=htmlfile, rootdir=".")
	#Produce the tabs	
	writeHTSAHtmlTab(enrichmentAnalysis=enrichmentAnalysis,htmlfile=htmlfile,rootdir=".",index=FALSE)
	l.HyperGeo.results<-length(enrichmentAnalysis$HyperGeo.results)
	l.GSEA.results<-length(enrichmentAnalysis$GSEA.results)
	l.Sig.adj.pvals.in.both<-length(enrichmentAnalysis$Sig.adj.pvals.in.both)
	#check that the arguments experimentName, geneListName have length one
	if(length(experimentName) != 1  || length(geneListName) != 1 ) 
		warning("The experimentName and geneListName should be of length one")
	cat("\n <hr/> \n <br> The enrichment analysis was performed using the data: ", append = TRUE, file = htmlfile)
	cat(geneListName, append = TRUE, file = htmlfile)
	cat(paste(" ( ",length(geneList), " genes)",sep=""), append = TRUE, file = htmlfile)
	cat("\n <br> This analysis was performed using the gene set collection(s): ", append = TRUE, file = htmlfile)
	cat("\n \t <UL>", append = TRUE, file = htmlfile)
	for(i in 1:number.gsc) {
		cat(paste(" \n \t \t <LI>",names(listOfGeneSetCollections)[i],sep=""), append = TRUE, file = htmlfile)
		cat(paste(" ( ",length(listOfGeneSetCollections[[i]]), " gene sets, of which ",dim(enrichmentAnalysis$GSEA.results[[i]])[1] ," were above the minimum size )",sep=""), append = TRUE, file = htmlfile)
	}
	cat("\n \t </UL>", append = TRUE, file = htmlfile)
	cat("\n <br> The following methods were used: ", append = TRUE, file = htmlfile)
	cat("\n \t <UL> \n \t \t <LI>", append = TRUE, file = htmlfile)
	cat("Hypergeometric test", append = TRUE, file = htmlfile)
	cat("\n \t \t <UL> \n \t \t \t <LI>", append = TRUE, file = htmlfile)
	#check that the arguments cutoffHits, p.adj.method, min.gene.set.size, nPermutations have length one	
	if(length(cutoffHits) != 1  || length(p.adj.method) != 1 || length(min.gene.set.size) != 1 || length(nPermutations) != 1) 
		warning("Thearguments cutoffHits, p.adj.method, min.gene.set.size, nPermutations should be of length one")
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
	if(is.na(networkObject)) {
		cat(paste("\n \t \t \t <LI> Interaction dataset: The Biogrid organism: ",species), append = TRUE, file = htmlfile)
	} else {
		cat(paste("\n \t \t \t <LI> Interaction dataset: Biogrid object: ",deparse(substitute(networkObject))), append = TRUE, file = htmlfile)
	}
	if(!genetic) {
		cat(" (excluding genetic interactions)", append = TRUE, file = htmlfile)
	} else {
		cat(" (including genetic interactions)", append = TRUE, file = htmlfile)
	}	
	cat(paste("\n \t \t \t <LI> FDR for score calculation: ",fdr), append = TRUE, file = htmlfile)
	cat("\n \t \t </UL>", append = TRUE, file = htmlfile)
	cat("\n \t </UL>", append = TRUE, file = htmlfile)
	cat('\n </tr> \n </table> \n </td> \n </tr> \n </table> \n </body> \n </html>',append = TRUE, file = htmlfile)
	##########################################
	#       produce htmls for hyperGeo	     #
	##########################################
	#-Produce the HyperG results pages: one page is produced for each element of the enrichmentAnalysis$HyperGeo.results list (i.e. for each gene set collection)
	for(n in 1:l.HyperGeo.results) {
		htmlfile =  file.path(htmdir,paste("hyperg",n,".html",sep=""))
		##########################################
		#       non-GO and non-KEGG Gene sets    #
		##########################################
		#First, take care of the non-GO and non-KEGG gene set collections		
		if(length(intersect(n,keggGS)) == 0 && length(intersect(n,GOGS)) == 0) {
			writeHTSAHtmlHead(experimentName=experimentName, htmlfile=htmlfile, rootdir="..")
			#Produce the tabs
			writeHTSAHtmlTab(enrichmentAnalysis=enrichmentAnalysis,htmlfile=htmlfile,rootdir="..",index=TRUE)
			#Replace '/' in gene names with '_'
			t.names<-rownames(enrichmentAnalysis$HyperGeo.results[[n]])
			t.names<-sapply(t.names, 
					function(t.name) {
						while(length(grep(pattern="/",x=t.name)) != 0)
							t.name<-sub(pattern="/",replacement="_",x=t.name,perl=TRUE)
						t.name
					}
			)	
			#Produce table
			#data table
			dat.tab<-signif(enrichmentAnalysis$HyperGeo.results[[n]],digits=4)
			dat.tab<-cbind(rownames(dat.tab),dat.tab)
			colnames(dat.tab)[1]<-"Gene Set name"
			
			this.row<-nrow(enrichmentAnalysis$HyperGeo.results[[n]])
			this.col<-ncol(enrichmentAnalysis$HyperGeo.results[[n]])
			#hyperlink table 
			href.tab<-array(NA,dim=c(this.row,this.col+1,3))
			dimnames(href.tab)[[3]]<-c("href","target","title")
			href.tab[1:nGseaPlots,6,1]<-paste("../doc/",t.names[1:nGseaPlots],".txt",sep="")
			href.tab[1:nGseaPlots,6,2]<-"_blank"
			href.tab[1:nGseaPlots,6,3]<-"Observed.hits"
			#highlight table
			signif.tab<-matrix(NA,this.row,this.col+1)
			colnames(signif.tab)<-rep("class",this.col+1)
			signif.tab[which(enrichmentAnalysis$HyperGeo.results[[n]][,8] < 0.05), 9]<-"signif"
			#row attribute table
			row.attr.tab<-matrix("even",this.row,1)
			colnames(row.attr.tab)<-"class"
			row.attr.tab[which(1:this.row%%2==1),1]<-"odd"
			#Generate and write table 
			writeHTSAHtmlTable(
					dat.tab=dat.tab, 
					href.tab=href.tab, 
					signif.tab=signif.tab, 
					row.attr.tab=row.attr.tab,
					tab.class="result",
					tab.name=paste(names(enrichmentAnalysis$HyperGeo.results)[n],' Hyperg. Tests',sep=""),
					htmlfile=htmlfile
			)
			writeHTSAHtmlTail(htmlfile=htmlfile)
		}
		##########################################
		#       		KEGG Gene sets 			 #
		##########################################
		#Second, take care of the KEGG gene set collections	
		if(length(intersect(n,keggGS)) != 0){
			writeHTSAHtmlHead(experimentName=experimentName, htmlfile=htmlfile, rootdir="..")
			#Produce the tabs
			writeHTSAHtmlTab(enrichmentAnalysis=enrichmentAnalysis,htmlfile=htmlfile,rootdir="..",index=TRUE)

			kegggsnames<-rownames(enrichmentAnalysis$HyperGeo.results[[n]])
			#Produce table
			#data table
			dat.tab<-signif(enrichmentAnalysis$HyperGeo.results[[n]],digits=4)
			dat.tab<-cbind(rownames(dat.tab),dat.tab)
			colnames(dat.tab)[1]<-"Gene Set name"
			
			this.row<-nrow(enrichmentAnalysis$HyperGeo.results[[n]])
			this.col<-ncol(enrichmentAnalysis$HyperGeo.results[[n]])
			#hyperlink table 
			href.tab<-array(NA,dim=c(this.row,this.col+1,3))
			dimnames(href.tab)[[3]]<-c("href","target","title")
			href.tab[,1,1]<-paste("http://www.genome.jp/dbget-bin/www_bget?pathway:",kegggsnames,sep="")
			href.tab[,1,2]<-"_blank"
			href.tab[,1,3]<-rownames(enrichmentAnalysis$HyperGeo.results[[n]])
			href.tab[1:nGseaPlots,6,1]<-paste("../doc/",kegggsnames[1:nGseaPlots],".txt",sep="")
			href.tab[1:nGseaPlots,6,2]<-"_blank"
			href.tab[1:nGseaPlots,6,3]<-"Observed.hits"
			#highlight table
			signif.tab<-matrix(NA,this.row,this.col+1)
			colnames(signif.tab)<-rep("class",this.col+1)
			signif.tab[which(enrichmentAnalysis$HyperGeo.results[[n]][,8] < 0.05), 9]<-"signif"
			#row attribute table
			row.attr.tab<-matrix("even",this.row,1)
			colnames(row.attr.tab)<-"class"
			row.attr.tab[which(1:this.row%%2==1),1]<-"odd"
			#Generate and write table 
			writeHTSAHtmlTable(
					dat.tab=dat.tab, 
					href.tab=href.tab, 
					signif.tab=signif.tab, 
					row.attr.tab=row.attr.tab,
					tab.class="result",
					tab.name=paste(names(enrichmentAnalysis$HyperGeo.results)[n],' Hyperg. Tests',sep=""),
					htmlfile=htmlfile
			)
			writeHTSAHtmlTail(htmlfile=htmlfile)
		}
		##########################################
		#       		GO Gene sets		     #
		##########################################
		#Third, take care of the GO gene set collections	
		if(length(intersect(n,GOGS)) != 0) {
			writeHTSAHtmlHead(experimentName=experimentName, htmlfile=htmlfile, rootdir="..")
			#Produce the tabs
			writeHTSAHtmlTab(enrichmentAnalysis=enrichmentAnalysis,htmlfile=htmlfile,rootdir="..",index=TRUE)			
			gogsnames2web<-sapply(rownames(enrichmentAnalysis$HyperGeo.results[[n]]), 
					function(gogsname) {
						sub(pattern="(\\D*$)",replacement="",x=sub(pattern="(\\D*)",replacement="",x=gogsname,perl=TRUE),perl=TRUE)
					}
			)
			gogsnames2doc<-sapply(rownames(enrichmentAnalysis$HyperGeo.results[[n]]), 
					function(gogsname) {
						sub(pattern="(\\D*$)",replacement="",x=gogsname,perl=TRUE)
					}
			)
			#Produce table
			#data table
			dat.tab<-signif(enrichmentAnalysis$HyperGeo.results[[n]],digits=4)
			dat.tab<-cbind(rownames(dat.tab),dat.tab)
			colnames(dat.tab)[1]<-"Gene Set name"
			
			this.row<-nrow(enrichmentAnalysis$HyperGeo.results[[n]])
			this.col<-ncol(enrichmentAnalysis$HyperGeo.results[[n]])
			#hyperlink table 
			href.tab<-array(NA,dim=c(this.row,this.col+1,3))
			dimnames(href.tab)[[3]]<-c("href","target","title")
			href.tab[,1,1]<-paste("http://www.ebi.ac.uk/QuickGO/GTerm?id=GO:",gogsnames2web,sep="")
			href.tab[,1,2]<-"_blank"
			href.tab[,1,3]<-rownames(enrichmentAnalysis$HyperGeo.results[[n]])
			href.tab[1:nGseaPlots,6,1]<-paste("../doc/",gogsnames2doc[1:nGseaPlots],".txt",sep="")
			href.tab[1:nGseaPlots,6,2]<-"_blank"
			href.tab[1:nGseaPlots,6,3]<-"Observed.hits"
			#highlight table
			signif.tab<-matrix(NA,this.row,this.col+1)
			colnames(signif.tab)<-rep("class",this.col+1)
			signif.tab[which(enrichmentAnalysis$HyperGeo.results[[n]][,8] < 0.05), 9]<-"signif"
			#row attribute table
			row.attr.tab<-matrix("even",this.row,1)
			colnames(row.attr.tab)<-"class"
			row.attr.tab[which(1:this.row%%2==1),1]<-"odd"
			#Generate and write table 
			writeHTSAHtmlTable(
					dat.tab=dat.tab, 
					href.tab=href.tab, 
					signif.tab=signif.tab, 
					row.attr.tab=row.attr.tab,
					tab.class="result",
					tab.name=paste(names(enrichmentAnalysis$HyperGeo.results)[n],' Hyperg. Tests',sep=""),
					htmlfile=htmlfile
			)
			writeHTSAHtmlTail(htmlfile=htmlfile)
		}
	}
	##########################################
	#       produce htmls for GSEA		     #
	##########################################
	#-Produce the gsea results pages:  one page is produced for each element of the enrichmentAnalysis$GSEA.results list (i.e. for each gene set collection)
	for(n in 1:l.GSEA.results) {
		htmlfile =  file.path(htmdir,paste("gsea",n,".html",sep=""))
		##########################################
		#     non-GO and non-KEGG Gene sets	     #
		##########################################
		#First, take care of the non-GO and non-KEGG gene set collections		
		if(length(intersect(n,keggGS)) == 0 && length(intersect(n,GOGS)) == 0) {
			writeHTSAHtmlHead(experimentName=experimentName, htmlfile=htmlfile, rootdir="..")
			#Produce the tabs
			writeHTSAHtmlTab(enrichmentAnalysis=enrichmentAnalysis,htmlfile=htmlfile,rootdir="..",index=TRUE)
			#Replace '/' in gene names with '_'
			t.names<-rownames(enrichmentAnalysis$GSEA.results[[n]])
			t.names<-sapply(t.names, 
					function(t.name) {
						while(length(grep(pattern="/",x=t.name)) != 0)
							t.name<-sub(pattern="/",replacement="_",x=t.name,perl=TRUE)
						t.name
					}
			)	
			#Produce table
			#data table
			this.row<-nrow(enrichmentAnalysis$GSEA.results[[n]])
			this.col<-ncol(enrichmentAnalysis$GSEA.results[[n]])
			dat.tab<-signif(enrichmentAnalysis$GSEA.results[[n]],digits=4)
			dat.tab<-cbind(rownames(dat.tab),dat.tab,rep("",this.row))
			dat.tab[1:nGseaPlots,this.col+2]<-"plot"
			colnames(dat.tab)[1]<-"Gene Set name"
			colnames(dat.tab)[this.col+2]<-"Plots"
			
			#hyperlink table 
			href.tab<-array(NA,dim=c(this.row,this.col+2,3))
			dimnames(href.tab)[[3]]<-c("href","target","title")
			href.tab[1:nGseaPlots,6,1]<-paste("../image/gsea_plots",t.names[1:nGseaPlots],".png",sep="")
			href.tab[1:nGseaPlots,6,2]<-"_blank"
			href.tab[1:nGseaPlots,6,3]<-"gseaplots"
			#highlight table
			signif.tab<-matrix(NA,this.row,this.col+2)
			colnames(signif.tab)<-rep("class",this.col+2)
			signif.tab[which(enrichmentAnalysis$GSEA.results[[n]][,3] < 0.05), 4]<-"signif"
			#row attribute table
			row.attr.tab<-matrix("even",this.row,1)
			colnames(row.attr.tab)<-"class"
			row.attr.tab[which(1:this.row%%2==1),1]<-"odd"
			#Generate and write table 
			writeHTSAHtmlTable(
					dat.tab=dat.tab, 
					href.tab=href.tab, 
					signif.tab=signif.tab, 
					row.attr.tab=row.attr.tab,
					tab.class="result",
					tab.name=paste(names(enrichmentAnalysis$GSEA.results)[n],' GSEA',sep=""),
					htmlfile=htmlfile
			)
			writeHTSAHtmlTail(htmlfile=htmlfile)			
		}
		##########################################
		#   		  Kegg Gene sets		     #
		##########################################
		#Second:take care of the Kegg gene sets		
		if(length(intersect(n,keggGS)) != 0){
			writeHTSAHtmlHead(experimentName=experimentName, htmlfile=htmlfile, rootdir="..")
			#Produce the tabs
			writeHTSAHtmlTab(enrichmentAnalysis=enrichmentAnalysis,htmlfile=htmlfile,rootdir="..",index=TRUE)
			kegggsnames<-rownames(enrichmentAnalysis$GSEA.results[[n]])
			#Produce table
			#data table
			this.row<-nrow(enrichmentAnalysis$GSEA.results[[n]])
			this.col<-ncol(enrichmentAnalysis$GSEA.results[[n]])
			dat.tab<-signif(enrichmentAnalysis$GSEA.results[[n]],digits=4)
			dat.tab<-cbind(rownames(dat.tab),dat.tab,rep("",this.row))
			dat.tab[1:nGseaPlots,this.col+2]<-"plot"
			colnames(dat.tab)[1]<-"Gene Set name"
			colnames(dat.tab)[this.col+2]<-"Plots"
			
			#hyperlink table 
			href.tab<-array(NA,dim=c(this.row,this.col+2,3))
			dimnames(href.tab)[[3]]<-c("href","target","title")
			href.tab[,1,1]<-paste("http://www.genome.jp/dbget-bin/www_bget?pathway:",kegggsnames,sep="")
			href.tab[,1,2]<-"_blank"
			href.tab[,1,3]<-rownames(enrichmentAnalysis$GSEA.results[[n]])
			href.tab[1:nGseaPlots,6,1]<-paste("../image/gsea_plots",kegggsnames[1:nGseaPlots],".png",sep="")
			href.tab[1:nGseaPlots,6,2]<-"_blank"
			href.tab[1:nGseaPlots,6,3]<-"gseaplots"
			#highlight table
			signif.tab<-matrix(NA,this.row,this.col+2)
			colnames(signif.tab)<-rep("class",this.col+2)
			signif.tab[which(enrichmentAnalysis$GSEA.results[[n]][,3] < 0.05), 4]<-"signif"
			#row attribute table
			row.attr.tab<-matrix("even",this.row,1)
			colnames(row.attr.tab)<-"class"
			row.attr.tab[which(1:this.row%%2==1),1]<-"odd"
			#Generate and write table 
			writeHTSAHtmlTable(
					dat.tab=dat.tab, 
					href.tab=href.tab, 
					signif.tab=signif.tab, 
					row.attr.tab=row.attr.tab,
					tab.class="result",
					tab.name=paste(names(enrichmentAnalysis$GSEA.results)[n],' GSEA',sep=""),
					htmlfile=htmlfile
			)
			writeHTSAHtmlTail(htmlfile=htmlfile)
		}
		##########################################
		#   		  	GO Gene sets		     #
		##########################################
		#Third, take care of the GO gene sets		
		if(length(intersect(n,GOGS)) != 0){
			writeHTSAHtmlHead(experimentName=experimentName, htmlfile=htmlfile, rootdir="..")
			#Produce the tabs
			writeHTSAHtmlTab(enrichmentAnalysis=enrichmentAnalysis,htmlfile=htmlfile,rootdir="..",index=TRUE)
			gogsnames2web<-sapply(rownames(enrichmentAnalysis$GSEA.results[[n]]), 
					function(gogsname) {
						sub(pattern="(\\D*$)",replacement="",x=sub(pattern="(\\D*)",replacement="",x=gogsname,perl=TRUE),perl=TRUE)
					}
			)
			gogsnames2doc<-sapply(rownames(enrichmentAnalysis$GSEA.results[[n]]), 
					function(gogsname) {
						sub(pattern="(\\D*$)",replacement="",x=gogsname,perl=TRUE)
					}
			)
			#Produce table
			#data table
			this.row<-nrow(enrichmentAnalysis$GSEA.results[[n]])
			this.col<-ncol(enrichmentAnalysis$GSEA.results[[n]])
			dat.tab<-signif(enrichmentAnalysis$GSEA.results[[n]],digits=4)
			dat.tab<-cbind(rownames(dat.tab),dat.tab,rep("",this.row))
			dat.tab[1:nGseaPlots,this.col+2]<-"plot"
			colnames(dat.tab)[1]<-"Gene Set name"
			colnames(dat.tab)[this.col+2]<-"Plots"
			#hyperlink table 
			href.tab<-array(NA,dim=c(this.row,this.col+2,3))
			dimnames(href.tab)[[3]]<-c("href","target","title")
			href.tab[,1,1]<-paste("http://www.ebi.ac.uk/QuickGO/GTerm?id=GO:",gogsnames2web,sep="")
			href.tab[,1,2]<-"_blank"
			href.tab[,1,3]<-rownames(enrichmentAnalysis$GSEA.results[[n]])
			href.tab[1:nGseaPlots,6,1]<-paste("../image/gsea_plots",gogsnames2doc[1:nGseaPlots],".png",sep="")
			href.tab[1:nGseaPlots,6,2]<-"_blank"
			href.tab[1:nGseaPlots,6,3]<-"gseaplots"
			#highlight table
			signif.tab<-matrix(NA,this.row,this.col+2)
			colnames(signif.tab)<-rep("class",this.col+2)
			signif.tab[which(enrichmentAnalysis$GSEA.results[[n]][,3] < 0.05), 4]<-"signif"
			#row attribute table
			row.attr.tab<-matrix("even",this.row,1)
			colnames(row.attr.tab)<-"class"
			row.attr.tab[which(1:this.row%%2==1),1]<-"odd"
			#Generate and write table 
			writeHTSAHtmlTable(
					dat.tab=dat.tab, 
					href.tab=href.tab, 
					signif.tab=signif.tab, 
					row.attr.tab=row.attr.tab,
					tab.class="result",
					tab.name=paste(names(enrichmentAnalysis$GSEA.results)[n],' GSEA',sep=""),
					htmlfile=htmlfile
			)
			writeHTSAHtmlTail(htmlfile=htmlfile)
		}	
	}	
	##########################################
	#  produce htmls for enrichment summary  #
	##########################################	
	#Produce the enrichment summary page
	for(n in 1:l.Sig.adj.pvals.in.both) {
		htmlfile =  file.path(htmdir,paste("enrichment",n,".html",sep=""))
		writeHTSAHtmlHead(experimentName=experimentName, htmlfile=htmlfile, rootdir="..")
		#Produce the tabs
		writeHTSAHtmlTab(enrichmentAnalysis=enrichmentAnalysis,htmlfile=htmlfile,rootdir="..",index=TRUE)
		#data table
		this.row<-nrow(enrichmentAnalysis$Sig.adj.pvals.in.both[[n]])
		this.col<-ncol(enrichmentAnalysis$Sig.adj.pvals.in.both[[n]])
		dat.tab<-signif(enrichmentAnalysis$Sig.adj.pvals.in.both[[n]],digits=4)
		dat.tab<-cbind(rownames(dat.tab),dat.tab)
		colnames(dat.tab)[1]<-"Gene Set name"
		#row attribute table
		row.attr.tab<-matrix("even",this.row,1)
		colnames(row.attr.tab)<-"class"
		row.attr.tab[which(1:this.row%%2==1),1]<-"odd"
		#Generate and write table 
		writeHTSAHtmlTable(
				dat.tab=dat.tab, 
				href.tab=NULL, 
				signif.tab=NULL, 
				row.attr.tab=row.attr.tab,
				tab.class="result",
				tab.name=paste("Gene sets with significant adjusted p value with both hypergeometric test and GSEA: ",names(enrichmentAnalysis$GSEA.results)[n],sep=""),
				htmlfile=htmlfile
		)
		writeHTSAHtmlTail(htmlfile=htmlfile)
	}
	
	##########################################
	#  produce htmls for network analyses	 #
	##########################################	
	#Produce the networks analysis results page
	nwAnalysisGraphFile<-"EnrichedSubNw.png"
	networkPlot(nwAnalysisOutput=nwAnalysisOutput,phenotypeVector=geneList,filepath=imgdir,filename=nwAnalysisGraphFile)
	htmlfile =  file.path(htmdir,"network.html")
	writeHTSAHtmlHead(experimentName=experimentName, htmlfile=htmlfile, rootdir="..")
	#Produce the tabs
	writeHTSAHtmlTab(enrichmentAnalysis=enrichmentAnalysis,htmlfile=htmlfile,rootdir="..",index=TRUE)
	#Check that the nwAnalysisOutput has the right format
	if(!is.list(nwAnalysisOutput)) stop("The nwAnalysisOutput should be a list")
	if(!all(names(nwAnalysisOutput) %in% c("subnw","labels"))) stop("The nwAnalysisOutput should contain the following elements: 'subnw', 'labels'")
	if(!is(nwAnalysisOutput$subnw,"graphNEL")) stop("The nwAnalysisOutput$subnw should be a graphNEL object")
	if(length(nwAnalysisOutput$labels) != length(nodes(nwAnalysisOutput$subnw))) stop("The nwAnalysisOutput$label should be a vector of the length nodes(nwAnalysisOutput$subnw)")
	EnrichSNnodes<-cbind(nodes(nwAnalysisOutput$subnw),nwAnalysisOutput$labels)
	colnames(EnrichSNnodes)<-c("Entrez Identifier","Symbol")
	this.row<-nrow(EnrichSNnodes)
	#row attribute table
	row.attr.tab<-matrix("even",this.row,1)
	colnames(row.attr.tab)<-"class"
	row.attr.tab[which(1:this.row%%2==1),1]<-"odd"
	cat('\n <table class="noframe"> <tr> <td>',append=TRUE,file=htmlfile)
	#Generate and write table 
	writeHTSAHtmlTable(
			dat.tab=EnrichSNnodes, 
			href.tab=NULL, 
			signif.tab=NULL, 
			row.attr.tab=row.attr.tab,
			tab.class="result",
			tab.name=paste('Click <a href="../image/',nwAnalysisGraphFile,'" target="_blank" title="Enriched Subnetwork">here</a> to get the enriched subnetwork ',sep=""),
			htmlfile=htmlfile
	)
	cat('</td> <td> <img src="../image/',nwAnalysisGraphFile,'" align="top" width="800" height="800"> </td> </tr> </table>',sep="",append=TRUE,file=htmlfile)
	writeHTSAHtmlTail(htmlfile=htmlfile)
}

