#This function takes the output of gseaScores and the gene list, and plot figures for GSEA.
gseaPlots <- function(runningScore, enrichmentScore, positions, geneList, output="png", filepath=".", filename="test") {
	#check input arguments
	paraCheck("genelist",geneList)
	paraCheck("output",output)
	paraCheck("filepath",filepath)
	paraCheck("filename",filename)
	#check that the 'runningScore' is a vector of length=length of geneList
	if(!is.numeric(runningScore) || length(runningScore)==0) 
		stop("'runningScore' should be a numerical vector!\n")
	#check that the 'positions' vector contains only zeros and ones, and is of the right length and class
	if(!(is.numeric(positions) || is.integer(positions)) || length(positions)==0) 
		stop("'positions' should be a numerical vector!\n")
	if(!all(unique(positions) %in% c(0,1))) 
		stop("'positions' should contains only 0s and 1s, see help(gseaScores)!\n")
	if(length(runningScore) != length(geneList)) 
		stop("The length of 'runningScore' should be the same as the length of 'geneList'!\n")
	if(length(positions) != length(geneList)) 
		stop("The length of 'positions' should be the same as the length of 'geneList'!\n")
	#open a file	
	if(output == "pdf" ) 
		pdf(file.path(filepath,paste("gsea_plots",filename,".pdf",sep="")))
	if(output == "png" ) 
		png(file.path(filepath,paste("gsea_plots",filename,".png",sep="")))
	#set the graphical parameters	
	par(pin=c(5,1.5),mfrow=c(2,1),lwd=1,mai=c(0.2,1,1,1))
	#Plot the phenotypes along the geneList, and add a vertical line for each match between geneList and gene set
	#this is done using the 'positions' output of gseaScores, which stores a one for each match position and a zero otherwise	
	plot(x=seq(1,length(geneList)),type="l",y=geneList,ylab="Phenotypes",xlab=NA,col="magenta",lwd=2,xaxt="n")
	abline(v=which(positions == 1))
	abline(h=0)
	#Plot the running score and add a vertical line at the position of the enrichment score (maximal absolute value of the running score)	
	par(mai=c(1,1,0.1,1))
	plot(x=c(1:length(runningScore)),y=runningScore,type="l",xlab="Position in the ranked list of genes",ylab="Running enrichment score")
	abline(h=0)
	abline(v=which(runningScore == enrichmentScore),lty=3,col="magenta",lwd=3)	
	dev.off()
}

