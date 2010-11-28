###############################################################################
# Write html reports
# Xin Wang <xw264@cam.ac.uk>
# Advisor: Florian Markowetz <florian.markowetz@cancer.org.uk> 
# University of Cambridge Deparment of Oncology
# Cancer Research UK - Cambridge Research Institute
# At 11:13:32, on 21 Nov 2010
###############################################################################
makeGSEAplots <- function(geneList, geneSet, exponent,
		filepath, filename,
		output='png'){
	test <- gseaScores(geneList=geneList, geneSet=geneSet,
			exponent=exponent, mode="graph")
	gseaPlots(runningScore=test[['runningScore']],
			enrichmentScore=test[['enrichmentScore']],
			positions=test[['positions']], geneList=geneList,
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

#generate html report for both GSCA and NWA
setClassUnion("GSCA_Or_NULL",c("GSCA","NULL"))
setClassUnion("NWA_Or_NULL",c("NWA","NULL"))
setMethod(
		"reportAll",
		c("GSCA_Or_NULL", "NWA_Or_NULL"),
		function(gsca, nwa, experimentName="Unknown", species=NULL, ntop=NULL, allSig=FALSE, keggGSCs=NULL, goGSCs=NULL, reportDir="HTSanalyzerReport", ...) {
			#call writeReportHTSA
			if(missing(gsca)) gsca<-NULL
			if(missing(nwa)) nwa<-NULL
			writeReportHTSA(gsca=gsca, nwa=nwa, experimentName=experimentName, species=species, ntop=ntop, allSig=allSig, 
					keggGSCs=keggGSCs, goGSCs=goGSCs, reportDir=reportDir)
		}
)
writeReportHTSA <- function(gsca=NULL, nwa=NULL, experimentName="Unknown", species=NULL, ntop=NULL, allSig=FALSE, keggGSCs=NULL, goGSCs=NULL, reportDir="HTSanalyzerReport", ...) {
	#check input arguments
	##check gsca and nwa
	Rep.gsca<-FALSE
	Rep.nwa<-FALSE
	if(is(gsca, "GSCA"))
		Rep.gsca<-TRUE
	if(is(nwa, "NWA"))
		Rep.nwa<-TRUE
	if(!Rep.gsca && !Rep.nwa)
		stop("Please input a 'GSCA' or 'NWA' object, or both of them! \n")
	##check experimentName
	paraCheck(name="experimentName",para=experimentName)
	if(Rep.gsca) {
		##check gscs
		gscs<-names(gsca@listOfGeneSetCollections)
		rslt.gscs<-names(gsca@result$GSEA.results)
		#if(missing(gscs) || !is.character(gscs) || length(gscs)==0)
		#	stop("Please specify the name(s) of Gene Set Collections in 'gscs'! \n")
		#if(!all(gscs %in% gsc.names))
	#	stop("Wrong Gene Set Collection name(s) in 'gscs'! \n")
	}
	##check species
	if(!is.null(species)) {
		paraCheck(name="species",para=species)	
	}
	##check ntop and allSig
	if(!is.null(ntop))
		paraCheck(name="ntop",para=ntop)
	paraCheck(name="allSig",para=allSig)
	if((is.null(ntop) && !allSig)||(!is.null(ntop) && allSig))
		stop("Either specify 'ntop' or set 'allSig' to be TRUE!\n")
	if(Rep.gsca) {
		##check keggGSCs and goGSCs
		if(!is.null(keggGSCs)) {
			paraCheck(name="keggGSCs",para=keggGSCs)
			if(!all(keggGSCs %in% gscs))
				stop("Wrong gene set collection names specified in 'keggGSCs'!\n")
		}
		if(!is.null(goGSCs)) {
			paraCheck(name="goGSCs",para=goGSCs)
			if(!all(goGSCs %in% gscs))
				stop("Wrong gene set collection names specified in 'goGSCs'!\n")
		}	
	}

	##check reportDir
	paraCheck(name="reportDir",para=reportDir)
	##########################################
	#           create directories 		     #
	##########################################
	## Directories: alphabetically sorted (1 - doc; 2 - html; 3 - image)
	dirs <- file.path(reportDir, c('','doc', 'html', 'image'))
	names(dirs) <- c('root', 'doc', 'html', 'image')
	sapply(dirs, function(diri) if(!file.exists(diri)) dir.create(diri))
	##########################################
	#       	produce html templates	     #
	##########################################
	#Copy css and logos in there
	cpfile<-dir(system.file("templates",package="HTSanalyzeR"),full=TRUE)
	file.copy(from=cpfile,to=dirs['image'],overwrite=TRUE)
	##########################################
	#          		 index page 		     #
	##########################################
	#Produce the index html page
	htmlfile<-file.path(reportDir,"index.html")
	writeHTSAHtmlHead(experimentName=experimentName, htmlfile=htmlfile, rootdir=".")
	#Produce the tabs	
	writeHTSAHtmlTab(enrichmentAnalysis=gsca@result,tab=c("GSCA","NWA")[c(Rep.gsca,Rep.nwa)],htmlfile=htmlfile,rootdir=".",index=FALSE)
	writeHTSAHtmlSummary(gsca=gsca, nwa=nwa, htmlfile=htmlfile)
	writeHTSAHtmlTail(htmlfile=htmlfile)
	##########################################
	#           	GSCA pages	 		     #
	##########################################
	if(Rep.gsca) {
		for(i in 1:length(rslt.gscs)) {
			##########################################
			#           	HyperGeo pages 		     #
			##########################################
			#write hits
			if(rslt.gscs[i] %in% gscs)
				writeHits(object=gsca, gscs=rslt.gscs[i], species=species, ntop=ntop, allSig=allSig, filepath=dirs['doc'])
			hyper.filenames<-getTopGeneSets(object=gsca, resultName="HyperGeo.results", gscs=rslt.gscs[i], ntop=ntop, allSig=allSig)
			htmlfile <-  file.path(dirs['html'],paste("hyperg",i,".html",sep=""))
			#create htmls
			writeHTSAHtmlHead(experimentName=experimentName, htmlfile=htmlfile, rootdir="..")
			#Produce the tabs
			writeHTSAHtmlTab(enrichmentAnalysis=gsca@result,tab=c("GSCA","NWA")[c(Rep.gsca,Rep.nwa)],htmlfile=htmlfile,rootdir="..",index=TRUE)
			#Produce table
			gs.names<-rownames(gsca@result$HyperGeo.results[[rslt.gscs[i]]])
			#data table
			dat.tab<-signif(gsca@result$HyperGeo.results[[rslt.gscs[i]]],digits=4)
			dat.tab<-cbind(rownames(dat.tab),dat.tab)
			colnames(dat.tab)[1]<-"Gene Set name"
			
			this.row<-nrow(gsca@result$HyperGeo.results[[rslt.gscs[i]]])
			this.col<-ncol(gsca@result$HyperGeo.results[[rslt.gscs[i]]])
			if(this.row>0) {
				#hyperlink table 
				href.tab<-array(NA,dim=c(this.row,this.col+1,3))
				dimnames(href.tab)[[3]]<-c("href","target","title")
				dimnames(href.tab)[[1]]<-gs.names
				##hyperlinks for kegg gene sets
				if(rslt.gscs[i] %in% keggGSCs) {
					href.tab[,1,1]<-paste("http://www.genome.jp/dbget-bin/www_bget?pathway:",gs.names,sep="")
					href.tab[,1,2]<-"_blank"
					href.tab[,1,3]<-gs.names
				}
				##hyperlinks for go gene sets
				if(rslt.gscs[i] %in% goGSCs) {
					gogsnames2web<-sapply(gs.names, 
						function(gogsname) {
							sub(pattern="(\\D*$)",replacement="",x=sub(pattern="(\\D*)",replacement="",x=gogsname,perl=TRUE),perl=TRUE)
						}
					)
					gogsnames2doc<-sapply(gs.names, 
						function(gogsname) {
							sub(pattern="(\\D*$)",replacement="",x=gogsname,perl=TRUE)
						}
					)
					href.tab[,1,1]<-paste("http://www.ebi.ac.uk/QuickGO/GTerm?id=GO:",gogsnames2web,sep="")
					href.tab[,1,2]<-"_blank"
					href.tab[,1,3]<-gs.names
				}
				href.tab[names(hyper.filenames[[1]]),6,1]<-paste("../doc/",hyper.filenames[[1]],".txt",sep="")
				href.tab[names(hyper.filenames[[1]]),6,2]<-"_blank"
				href.tab[names(hyper.filenames[[1]]),6,3]<-"Observed.hits"
				#highlight table
				signif.tab<-matrix(NA,this.row,this.col+1)
				colnames(signif.tab)<-rep("class",this.col+1)
				#rownames(signif.tab)<-rownames(gsca@result$HyperGeo.results[[gscs[i]]])
				signif.tab[which(gsca@result$HyperGeo.results[[rslt.gscs[i]]][,7] < gsca@para$pValueCutoff), 8]<-"signif"
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
						tab.name=paste(rslt.gscs[i],' Hyperg. Tests',sep=""),
						htmlfile=htmlfile
				)
			}
			writeHTSAHtmlTail(htmlfile=htmlfile)
			##########################################
			#           	GSEA pages	 		     #
			##########################################
			if(rslt.gscs[i] %in% gscs)
				plotGSEA(object=gsca, gscs=rslt.gscs[i], ntop=ntop, allSig=allSig, filepath=dirs['image'])
			gsea.filenames<-getTopGeneSets(object=gsca, resultName="GSEA.results", gscs=rslt.gscs[i], ntop=ntop, allSig=allSig)
			htmlfile <-  file.path(dirs['html'],paste("gsea",i,".html",sep=""))
			#create htmls
			writeHTSAHtmlHead(experimentName=experimentName, htmlfile=htmlfile, rootdir="..")
			#Produce the tabs
			writeHTSAHtmlTab(enrichmentAnalysis=gsca@result,tab=c("GSCA","NWA")[c(Rep.gsca,Rep.nwa)],htmlfile=htmlfile,rootdir="..",index=TRUE)
			#Produce table
			gs.names<-rownames(gsca@result$GSEA.results[[rslt.gscs[i]]])
			top.gs.id<-match(names(gsea.filenames[[1]]),gs.names)
			#data table
			this.row<-nrow(gsca@result$GSEA.results[[rslt.gscs[i]]])
			this.col<-ncol(gsca@result$GSEA.results[[rslt.gscs[i]]])
			if(this.row>0) {
				dat.tab<-signif(gsca@result$GSEA.results[[rslt.gscs[i]]],digits=4)
				dat.tab<-cbind(rownames(dat.tab),dat.tab,rep("",this.row))
				dat.tab[top.gs.id,this.col+2]<-"plot"
				colnames(dat.tab)[1]<-"Gene Set name"
				colnames(dat.tab)[this.col+2]<-"Plots"
				#hyperlink table 
				href.tab<-array(NA,dim=c(this.row,this.col+2,3))
				dimnames(href.tab)[[3]]<-c("href","target","title")				
				##hyperlinks for kegg gene sets
				if(rslt.gscs[i] %in% keggGSCs) {
					href.tab[,1,1]<-paste("http://www.genome.jp/dbget-bin/www_bget?pathway:",gs.names,sep="")
					href.tab[,1,2]<-"_blank"
					href.tab[,1,3]<-gs.names
				}
				##hyperlinks for go gene sets
				if(rslt.gscs[i] %in% goGSCs) {
					gogsnames2web<-sapply(gs.names, 
							function(gogsname) {
								sub(pattern="(\\D*$)",replacement="",x=sub(pattern="(\\D*)",replacement="",x=gogsname,perl=TRUE),perl=TRUE)
							}
					)
					gogsnames2doc<-sapply(gs.names, 
							function(gogsname) {
								sub(pattern="(\\D*$)",replacement="",x=gogsname,perl=TRUE)
							}
					)
					href.tab[,1,1]<-paste("http://www.ebi.ac.uk/QuickGO/GTerm?id=GO:",gogsnames2web,sep="")
					href.tab[,1,2]<-"_blank"
					href.tab[,1,3]<-gs.names
				}
				href.tab[top.gs.id,6,1]<-paste("../image/gsea_plots",gsea.filenames[[1]],".png",sep="")
				href.tab[top.gs.id,6,2]<-"_blank"
				href.tab[top.gs.id,6,3]<-"gseaplots"
				#highlight table
				signif.tab<-matrix(NA,this.row,this.col+2)
				colnames(signif.tab)<-rep("class",this.col+2)
				signif.tab[which(gsca@result$GSEA.results[[rslt.gscs[i]]][,3] < gsca@para$pValueCutoff), 4]<-"signif"
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
						tab.name=paste(rslt.gscs[i],' GSEA',sep=""),
						htmlfile=htmlfile
				)
			}
			writeHTSAHtmlTail(htmlfile=htmlfile)	
			##########################################
			#	     enrichment summary pages	 	 #
			##########################################
			htmlfile =  file.path(dirs['html'],paste("enrichment",i,".html",sep=""))
			writeHTSAHtmlHead(experimentName=experimentName, htmlfile=htmlfile, rootdir="..")
			#Produce the tabs
			writeHTSAHtmlTab(enrichmentAnalysis=gsca@result,tab=c("GSCA","NWA")[c(Rep.gsca,Rep.nwa)],htmlfile=htmlfile,rootdir="..",index=TRUE)
			#data table
			this.row<-nrow(gsca@result$Sig.adj.pvals.in.both[[rslt.gscs[i]]])
			this.col<-ncol(gsca@result$Sig.adj.pvals.in.both[[rslt.gscs[i]]])
			if(this.row>0) {
				dat.tab<-signif(gsca@result$Sig.adj.pvals.in.both[[rslt.gscs[i]]],digits=4)
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
						tab.name=paste("Gene sets with significant adjusted p value with both hypergeometric test and GSEA: ",rslt.gscs[i],sep=""),
						htmlfile=htmlfile
				)
			}
			writeHTSAHtmlTail(htmlfile=htmlfile)
		}
	}
	##########################################
	#           	NWA pages	 		     #
	##########################################
	if(Rep.nwa) {
		nwAnalysisGraphFile<-"EnrichedSubNw.png"
		#plot subnetwork
		plotSubNet(nwa, filepath=dirs['image'], filename=paste(nwAnalysisGraphFile,sep=""))
		htmlfile =  file.path(dirs['html'],"network.html")
		writeHTSAHtmlHead(experimentName=experimentName, htmlfile=htmlfile, rootdir="..")
		#Produce the tabs
		writeHTSAHtmlTab(enrichmentAnalysis=gsca@result,tab=c("GSCA","NWA")[c(Rep.gsca,Rep.nwa)],htmlfile=htmlfile,rootdir="..",index=TRUE)	
		cat(paste('\n <hr/> \n <br>','Click <a href="../doc/subnetNodes.txt" target="_blank" title="Enriched Subnetwork">here</a> to get Entrez identifiers and symbols of nodes in identified subnetwork! <br> \n\n',sep=""),
				file=htmlfile,append=TRUE)
		cat('<table class="result"><tr><td> <img src="../image/',nwAnalysisGraphFile,'" align="top" width="800" height="800"> </td></tr></table>',sep="",append=TRUE,file=htmlfile)
		writeHTSAHtmlTail(htmlfile=htmlfile)	
		htmlfile =  file.path(dirs['html'],"subnetNodes.html")
		
		#Check that the nwAnalysisOutput has the right format
		if(!is.null(nwa@result$label)) {
			EnrichSNnodes<-cbind(nodes(nwa@result$subnw),nwa@result$labels)
			colnames(EnrichSNnodes)<-c("Entrez Identifier","Symbol")
		} else {
			EnrichSNnodes<-matrix(nodes(nwa@result$subnw),ncol=1)
			colnames(EnrichSNnodes)<-c("Entrez Identifier")			
		}
		this.row<-nrow(EnrichSNnodes)
		if(this.row>0) {
			write.table(EnrichSNnodes,file=file.path(dirs['doc'],"subnetNodes.txt"),row.names=FALSE,quote=FALSE,sep="\t")			
		}
	}
	
}



