#This function produces a couple of GSEA plot, for two given phenotypes and one gene set 
#it takes as input:
	#-gl1:a named vector where names are gene identifiers of the same type as the ones in the 
#gene set collection, and values are the measurement on phenotype 2 corresponding to those genes
#this vector MUST be ordered
	#-gl2:a named vector where names are gene identifiers of the same type as the ones in the 
#gene set collection, and values are the measurement on phenotype 2 corresponding to those genes
#this vector MUST be ordered
	#-exponent: exponent for the GSEA, keep it at one (see the help on function gseaScores)
	#-geneSet: a vector of identifiers (see the help on function gseaScores)
	#-a name (single character vector) that characterizes the plot, e.g. the gene set name 
	#-output: either "png" or "pdf"
#it produces: a single plot	
pairwiseGseaPlot<-function(gl1,gl2,geneSet,exponent=1,output="png",filepath,filename){
	#check that the name argument has length 1 (it could be of any class)
	if(length(filename) != 1)
		warning("The name argument should have length one")
	#The other arguments will be tested by the gseaScores function	
	#Compute the score, running sum and positions
	gseaPh1<-gseaScores(geneList=gl1,geneSet=geneSet,exponent=exponent,mode="graph")
	gseaPh2<-gseaScores(geneList=gl2,geneSet=geneSet,exponent=exponent,mode="graph")
	#open a file	
	if(output == "pdf") 
		pdf(file=file.path(filepath,paste("pairwiseGsea_plots",filename,".pdf",sep="")))
	if(output == "png") 
		png(file=file.path(filepath,paste("pairwiseGsea_plots",filename,".png",sep="")))
	#set the graphical parameters	
	par(pin=c(5,1.5),mfrow=c(2,2),lwd=1,mai=c(0.2,0.8,0.5,0.1))
	#Plot the phenotypes along the geneList, and add a vertical line for each match between geneList and gene set
	#this is done using the 'Positions' output of gseaScores, which stores a one for each match position and a zero otherwise	
	plot(x=seq(1,length(gl1)),type="l",y=gl1,ylab="Phenotypes",xlab=NA,col="magenta",lwd=2,xaxt="n",main="Phenotype 1")
	abline(v=which(gseaPh1$Positions == 1))
	abline(h=0)
	par(mai=c(0.2,0.8,0.5,0.1))
	plot(x=seq(1,length(gl2)),type="l",y=gl2,ylab="Phenotypes",xlab=NA,col="magenta",lwd=2,xaxt="n",main="Phenotype 2")
	abline(v=which(gseaPh2$Positions == 1))
	abline(h=0)
	#Plot the running score and add a vertical line at the position of the enrichment score (maximal absolute value of the running score)	
	par(mai=c(1,0.8,0.1,0.1))
	plot(x=c(1:length(gseaPh1$Running.Score)),y=gseaPh1$Running.Score,type="l",xlab="Position in the ranked list of genes",ylab="Running enrichment score")
	abline(h=0)
	abline(v=which(gseaPh1$Running.Score == gseaPh1$Enrichment.Score),lty=3,col="magenta",lwd=3)	
	par(mai=c(1,0.8,0.1,0.1))
	plot(x=c(1:length(gseaPh2$Running.Score)),y=gseaPh2$Running.Score,type="l",xlab="Position in the ranked list of genes",ylab="Running enrichment score")
	abline(h=0)
	abline(v=which(gseaPh2$Running.Score == gseaPh2$Enrichment.Score),lty=3,col="magenta",lwd=3)
	dev.off()
}
