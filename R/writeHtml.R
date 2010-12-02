##A set of functions for writing different parts of a HTML web page
##Write head (html routines, logos, etc.) of report htmls
writeHTSAHtmlHead<-function(experimentName, htmlfile, rootdir=".") {
	##logos and templates
	cat('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"http://www.w3.org/TR/html4/loose.dtd">', file = htmlfile)
	cat('\n <html> \n <link rel="stylesheet" href="',file.path(rootdir,"image","htsanalyzer.css"),'" type="text/css"> ',append = TRUE, file = htmlfile)
	cat('\n <head> <title> HTSanalyzeR Experiment Report </title> </head>',append = TRUE, file = htmlfile)
	cat('\n <body> \n  <table class="border"> \n <tr class="border top"> \n <td class="border corner"> &nbsp&nbsp&nbsp&nbsp </td>',append = TRUE,file=htmlfile)
	cat(paste('\n <td class="border top"> <div class="header"> Report for Experiment <span class="header"> ',experimentName,'</span> </div>',sep=""),append = TRUE,file=htmlfile)
	cat(paste('<div class="timestamp">generated  ',date(),'(<small>version 1.2.3</small>) </div>',sep=""),append = TRUE,file=htmlfile)
	cat('\n <div class="HTSheader"> HTSanalyzeR </div> </td> </tr> <tr class="border middle"> <td class="border left"></td>',append = TRUE,file=htmlfile)
	cat('<td class="main"> <table> <tr> <div class="HTSlogos"> <img src="',file.path(rootdir,"image","Rlogo.png"),'" width="50" height="40"/>',append = TRUE,file=htmlfile)
	cat('&nbsp <img src="',file.path(rootdir,"image","blue_cruklogo.gif"),'" width="120" height="50"/>',append = TRUE,file=htmlfile)
	cat('&nbsp <img src="',file.path(rootdir,"image","goatcomputer.png"),'" width="85" height="60"/> </div></tr><tr>',append = TRUE,file=htmlfile)
}
##Write tabs of report htmls
writeHTSAHtmlTab<-function(enrichmentAnalysis,tab=c("GSCA","NWA"),htmlfile,rootdir=".",index=TRUE) {
	if(index) {
		cat(paste('\n <table class="noframe"> <tr> <td class="tabs"> <h3><a href="',file.path(rootdir,"index.html"),'" title="index">Index</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)		
	}
	if("GSCA" %in% tab) {
		cat("\n <tr>", append = TRUE, file = htmlfile)
		cat('<td class="tabs"> <h3>  </h3> </td>', append = TRUE, file = htmlfile)
		l.HyperGeo.results<-length(enrichmentAnalysis$HyperGeo.results)
		l.GSEA.results<-length(enrichmentAnalysis$GSEA.results)
		l.Sig.adj.pvals.in.both<-length(enrichmentAnalysis$Sig.adj.pvals.in.both)
		for(i in 1:l.HyperGeo.results){
			cat(paste('\n <td class="tabs"> <h3>',names(enrichmentAnalysis$HyperGeo.results)[i],'</h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
		cat("\n </tr> <tr>", append = TRUE, file = htmlfile)	
		cat('<td class="tabs"> <h3> Hypergeometric tests </h3> </td>', append = TRUE, file = htmlfile)
		for(i in 1:l.HyperGeo.results){
			cat(paste('\n <td class="tabs"> <h3><a href="',file.path(rootdir,"html",paste("hyperg",i,".html",sep="")),'" title="',sep=""), append = TRUE, file = htmlfile)
			cat(paste(names(enrichmentAnalysis$HyperGeo.results)[i],' Hyperg. Tests">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
		cat("\n </tr> <tr>", append = TRUE, file = htmlfile)
		cat('<td class="tabs"> <h3> GSEA </h3> </td>', append = TRUE, file = htmlfile)	
		for(i in 1:l.GSEA.results){
			cat(paste('\n <td class="tabs"> <h3><a href="',file.path(rootdir,"html",paste("gsea",i,".html",sep="")),'" title="',sep=""), append = TRUE, file = htmlfile)
			cat(paste(names(enrichmentAnalysis$GSEA.results)[i],' GSEA">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
		cat("\n </tr> <tr>", append = TRUE, file = htmlfile)	
		cat('<td class="tabs"> <h3> Enrichment Summary </h3> </td>', append = TRUE, file = htmlfile)
		for(i in 1:l.Sig.adj.pvals.in.both){
			cat(paste('\n <td class="tabs"> <h3><a href="',file.path(rootdir,"html",paste("enrichment",i,".html",sep="")),'" title="',sep=""), append = TRUE, file = htmlfile)
			cat(paste(names(enrichmentAnalysis$Sig.adj.pvals.in.both)[i],' Enrichment.summary">','here</a></h3> </td>',sep=""), append = TRUE, file = htmlfile)
		}
		cat("\n </tr>", append = TRUE, file = htmlfile)
	}
	if("NWA" %in% tab) {
		cat("\n <tr>", append = TRUE, file = htmlfile)
		cat(paste('\n  <td class="tabs"> <h3><a href="',file.path(rootdir,"html","network.html"),'" title="network">Network Analysis</a></h3> </td> </tr> ',sep=""), append = TRUE, file = htmlfile)
	}
	cat("\n </table>", append = TRUE, file = htmlfile)
}
##Write tail part of report htmls
writeHTSAHtmlTail<-function(htmlfile) {
	cat('\n </tr> \n </table> \n </td> \n </tr> \n </table> \n </body> \n </html>',append = TRUE, file = htmlfile)
}
##Write the summary page
writeHTSAHtmlSummary<-function(gsca=NULL, nwa=NULL, htmlfile) {
	cat("\n <hr/> \n ",append=TRUE,file=htmlfile)
	if(!is.null(gsca)) {
		l.HyperGeo.results<-length(gsca@result$HyperGeo.results)
		l.GSEA.results<-length(gsca@result$GSEA.results)
		l.Sig.adj.pvals.in.both<-length(gsca@result$Sig.adj.pvals.in.both)
		cat("<br> The enrichment analysis was performed using the gene list ", append = TRUE, file = htmlfile)
		##cat(geneListName, append = TRUE, file = htmlfile)
		cat(paste(" including ",gsca@summary$gl[,"input"], " genes",sep=""), append = TRUE, file = htmlfile)
		cat("\n <br><br> This analysis was performed using the gene set collection(s): ", append = TRUE, file = htmlfile)
		cat("\n \t <UL>", append = TRUE, file = htmlfile)
		for(i in 1:length(gsca@listOfGeneSetCollections)) {
			cat(paste(" \n \t \t <LI>",names(gsca@listOfGeneSetCollections)[i],sep=""), append = TRUE, file = htmlfile)
			cat(paste(" ( ",length(gsca@listOfGeneSetCollections[[i]]), " gene sets, of which ",dim(gsca@result$GSEA.results[[i]])[1] ," were above the minimum size )",sep=""), append = TRUE, file = htmlfile)
		}
		cat("\n \t </UL>", append = TRUE, file = htmlfile)
		cat("\n The following methods were used: ", append = TRUE, file = htmlfile)
		cat("\n \t <UL> \n \t \t <LI>", append = TRUE, file = htmlfile)
		cat("Hypergeometric test", append = TRUE, file = htmlfile)
		cat("\n \t \t <UL>", append = TRUE, file = htmlfile)
		##cat(paste("\n \t \t \t <LI> Cutoff for hits: ",gsca@summary$para$hypergeo[,"pValueCutoff"],sep=""), append = TRUE, file = htmlfile)
		cat(paste("\n \t \t \t <LI> MHT correction method: ",gsca@summary$para$hypergeo[,"pAdjustMethod"]), append = TRUE, file = htmlfile)
		cat(paste("\n \t \t \t <LI> Minimum gene set size: ",gsca@summary$para$hypergeo[,"minGeneSetSize"]), append = TRUE, file = htmlfile)
		cat("\n \t \t </UL>", append = TRUE, file = htmlfile)
		cat("\n \t \t <LI>", append = TRUE, file = htmlfile)
		cat("Gene Set Enrichment Analysis", append = TRUE, file = htmlfile)
		cat("\n \t \t <UL> \n \t \t \t <LI>", append = TRUE, file = htmlfile)
		cat(paste("Minimum gene set size: ",gsca@summary$para$gsea[,"minGeneSetSize"],sep=""), append = TRUE, file = htmlfile)
		cat(paste("\n \t \t \t <LI> MHT correction method: ",gsca@summary$para$gsea[,"pAdjustMethod"]), append = TRUE, file = htmlfile)
		cat(paste("\n \t \t \t <LI> Number of permutations: ",gsca@summary$para$gsea[,"nPermutations"]), append = TRUE, file = htmlfile)
		cat(paste("\n \t \t \t <LI> Exponent: ",gsca@summary$para$gsea[,"exponent"]), append = TRUE, file = htmlfile)
		cat("\n \t \t </UL>", append = TRUE, file = htmlfile)
	}

	if(!is.null(nwa)) {
		cat("\n \t \t <LI>", append = TRUE, file = htmlfile)
		cat("Network Analysis", append = TRUE, file = htmlfile)
		cat("\n \t \t <UL>", append = TRUE, file = htmlfile)

		if(!is.na(nwa@summary$db[,"species"])) {
			cat(paste("\n \t \t \t <LI> Interaction dataset: The Biogrid organism: ",ifelse(!is.na(nwa@summary$db[,"species"]),nwa@summary$db,"Unknown")), append = TRUE, file = htmlfile)
		} else {
			cat(paste("\n \t \t \t <LI> Interaction dataset: User-input Biogrid object "), append = TRUE, file = htmlfile)
		}
		if(!is.na(nwa@summary$db[,"genetic"])) {
			if(nwa@summary$db[,"genetic"]=="FALSE") {
				cat(" (excluding genetic interactions)", append = TRUE, file = htmlfile)
			} else {
				cat(" (including genetic interactions)", append = TRUE, file = htmlfile)
			}
		}
		cat(paste("\n \t \t \t <LI> FDR for score calculation: ",nwa@fdr), append = TRUE, file = htmlfile)
		cat("\n \t \t </UL>", append = TRUE, file = htmlfile)
		cat("\n \t </UL>", append = TRUE, file = htmlfile)
	}

}
##Write a table
writeHTSAHtmlTable<-function(dat.tab, href.tab=NULL, signif.tab=NULL, row.attr.tab=NULL, tab.class,tab.name, htmlfile) {
	##produce the name part of the table and append it to file 
	cat(paste('\n <hr/> \n <br>', tab.name,' <br>','\n',sep=""),"\n",file=htmlfile,append=TRUE)
	if(!is.null(dat.tab)) {
		if(!is.matrix(dat.tab))
			stop("input dat.tab must be a matrix including all data to display in the table")
		if(!is.null(href.tab) && length(dim(href.tab))!=3) 
			stop("input href.tab must be a 3d array specifying href, target and title of hyperlinks")
		if(!is.null(signif.tab) && !is.matrix(signif.tab))
			stop("input signif.tab must be a matrix specifying highlighted units")
		if(!is.null(row.attr.tab) && !is.matrix(row.attr.tab))
			stop("input row.attr.tab must be a nx1 matrix specifying classes of units")
		n.col<-ncol(dat.tab)
		n.row<-nrow(dat.tab)
		##produce head of the table
		tab.header<-sapply(colnames(dat.tab),GenHTSAHtmlRowUnit,header=TRUE)
		tab.header<-paste('<tr class="head">',paste(tab.header,collapse="\n"), '</tr>',sep="\n")
		tab.rows<-t(
				sapply(1:n.row,
						function(i) {
							sapply(1:n.col, 
									function(j) {
										GenHTSAHtmlRowUnit(content=dat.tab[i,j], href.attr=href.tab[i,j,], td.attr=signif.tab[i,j],header=FALSE)
									}
							)
						}
				)
		)
		##produce rows
		tab.rows<-apply(tab.rows, 1, paste, collapse="\n")
		tab.rows<-paste(paste('<tr class="',row.attr.tab,'">',sep=""), tab.rows,'</tr>',sep="\n")
		tab.rows<-unlist(lapply(list(tab.rows),paste,collapse="\n"))
		##paste head and rows together
		tab<-paste(paste('<table class="',tab.class,'">',sep=""),tab.header,tab.rows,"</table>",sep="\n")
		##append this table to file
		cat(tab,append=TRUE,file=htmlfile)
	}	
}
##Generate a unit of a table according to specified attributes
GenHTSAHtmlRowUnit<-function(content, href.attr=NULL,td.attr=NULL,header=FALSE) {
	mark<-ifelse(header,"th","td")

	row.unit<-paste(
			"<",mark, 
			ifelse(is.null(td.attr) || all(is.na(td.attr)),"",htmlAttrVectorPaste(td.attr)), ">", 
			ifelse(is.null(href.attr) || all(is.na(href.attr)),content,paste("<a", htmlAttrVectorPaste(href.attr), ">", content, "</a>", sep="")),
			"</",mark,">", sep="")
	return(row.unit)
}
##Collapse an attribute vector for a table unit
htmlAttrVectorPaste<-function(vec) {

	if(all(!is.na(vec)) && !is.null(vec)) {
		if(!(all(!is.na(names(vec))) && !is.null(names(vec)) && is.character(vec) && length(vec)>0)) {
			stop("Wrong vector of attributes")
		} else {
			attr<-unlist(lapply(list(paste(names(vec),paste('"',vec,'"',sep=""),sep="=")),paste,collapse=" "))
			return(paste(" ",attr,sep=""))
		}
	} else {
		return("")
	}		
}




