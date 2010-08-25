#Function that takes the output of gseaScores and the gene list, and does the plots for GSEA.
gseaPlots <-
function(Running.Score,Enrichment.Score,Positions,geneList,output="png",name){
#check that the geneList is a single numerical vector
	if(class(geneList) != "numeric" && class(geneList) != "integer") stop("The geneList should be a single vector of numerical phenotypes")
#check that the 'Running.Score' is a vector of length=length of geneList
	if(class(Running.Score) != "numeric") stop("The Running.Score should be a numerical vector")
	if(length(Running.Score) != length(geneList)) warning("The length of the Running.Score vector should be the same as the length of the geneList")
#check that the 'Positions' vector contains only zeros and ones, and is of the right length and class
	if(length(Positions) != length(geneList)) warning("The length of the Positions vector should be the same as the length of the geneList")
	if(class(Positions) != "numeric") stop("The Positions vector should be a single numerical vector")
	if(sum(unique(Positions)) != 1) warning("The Positions vector should contains only 0s and 1s, see help(gseaScores)")
#check that the output argument is correct (if it is not then the graph should simply be plotted on the screen since no file will be opened)	
	if(class(output) != "character" | length(intersect(output, c("png","pdf"))) == 0) {
		warning("If the 'output' argument is not specified as one of the following character strings: 'png' or 'pdf', no graphical file will be saved")
		}
#check that the name argument has length 1 (it could be of any class)
	if(length(name) != 1) {
		warning("The name argument should have length one")
		}
#open a file	
	if(output == "pdf" ) {
		pdf(paste("gsea_plots",name,".pdf",sep=""))
		}
	if(output == "png" ) {
		png(paste("gsea_plots",name,".png",sep=""))
		}
#set the graphical parameters	
	par(pin=c(5,1.5),mfrow=c(2,1),lwd=1,mai=c(0.2,1,1,1))
#Plot the phenotypes along the geneList, and add a vertical line for each match between geneList and gene set
#this is done using the 'Positions' output of gseaScores, which stores a one for each match position and a zero otherwise	
	plot(x=seq(1,length(geneList)),type="l",y=geneList,ylab="Phenotypes",xlab=NA,
		col="magenta",lwd=2,xaxt="n")
	abline(v=which(Positions == 1))
	abline(h=0)
#Plot the running score and add a vertical line at the position of the enrichment score (maximal absolute value of the running score)	
	par(mai=c(1,1,0.1,1))
	plot(x=c(1:length(Running.Score)),y=Running.Score,type="l",
		xlab="Position in the ranked list of genes",ylab="Running enrichment score")
	abline(h=0)
	abline(v=which(Running.Score == Enrichment.Score),lty=3,col="magenta",lwd=3)	
	dev.off()
	}

