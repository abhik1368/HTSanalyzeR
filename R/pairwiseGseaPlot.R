##This function produces plots from the GSEA analysis, in parallel for 
##two phenotypes for one gene set.

pairwiseGseaPlot <- function(gl1, gl2, geneSet, exponent = 1, 
	output = "png", filepath, filename, ...){
	##check arguments
	paraCheck("genelist", gl1)
	paraCheck("genelist", gl2)
	paraCheck("gs", geneSet)
	paraCheck("exponent", exponent)
	paraCheck("output", output)
	paraCheck("filepath", filepath)
	paraCheck("filename", filename)
	##The other arguments will be tested by the gseaScores function	
	##Compute the score, running sum and positions
	gseaPh1 <- gseaScores(geneList = gl1, geneSet = geneSet, 
		exponent = exponent, mode = "graph")
	gseaPh2 <- gseaScores(geneList = gl2, geneSet = geneSet, 
		exponent = exponent, mode = "graph")
	##open a file	
	if(output == "pdf") 
		pdf(file = file.path(filepath, paste("pairwiseGsea_plots", 
			filename, sep = "")), ...=...)
	if(output == "png") 
		png(filename = file.path(filepath, paste("pairwiseGsea_plots", 
			filename, sep = "")), ...=...)
	##set graphical parameters	
	par(pin = c(5, 1.5), mfrow = c(2, 2), lwd = 1, mai = c(0.2, 0.8, 0.5, 0.1))
	##Plot the phenotypes along the geneList, and add a vertical line 
	##for each match between geneList and gene set. This is done using 
	##the 'positions' output of gseaScores, which stores a one for each 
	##match position and a zero otherwise	
	plot(x = seq(1, length(gl1)), type = "l", y = gl1, ylab = "Phenotypes", 
		xlab = NA, col = "magenta", lwd = 2, xaxt = "n", main = "Phenotype 1")
	abline(v = which(gseaPh1$positions == 1))
	abline(h = 0)
	par(mai = c(0.2, 0.8, 0.5, 0.1))
	plot(x = seq(1, length(gl2)), type = "l", y = gl2, ylab = "Phenotypes", 
		xlab = NA, col = "magenta", lwd = 2, xaxt = "n", main = "Phenotype 2")
	abline(v = which(gseaPh2$positions == 1))
	abline(h = 0)
	##Plot the running score and add a vertical line at the position of 
	##the enrichment score (maximal absolute value of the running score)	
	par(mai = c(1, 0.8, 0.1, 0.1))
	plot(x = c(1:length(gseaPh1$runningScore)), y = gseaPh1$runningScore, 
		type = "l", xlab = "Position in the ranked list of genes", 
		ylab = "Running enrichment score")
	abline(h = 0)
	abline(v = which(gseaPh1$runningScore == gseaPh1$enrichmentScore), 
		lty = 3, col = "magenta", lwd = 3)	
	par(mai = c(1, 0.8, 0.1, 0.1))
	plot(x = c(1:length(gseaPh2$runningScore)), y = gseaPh2$runningScore, 
		type = "l", xlab = "Position in the ranked list of genes", 
		ylab = "Running enrichment score")
	abline(h = 0)
	abline(v = which(gseaPh2$runningScore == gseaPh2$enrichmentScore), 
		lty = 3, col = "magenta", lwd = 3)
	dev.off()
}
