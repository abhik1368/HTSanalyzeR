##This function takes in a single gene set (GeneSet), a vector 
##(GeneList) of gene symbols for all tested genes, a vector of "hits" 
##(hits), and a p-value adjustment method. It outputs a vector 
##containing the size of the gene universe, the size of the gene set 
##within this universe (i.e. how many genes from the universe map to 
##this gene set), the total number of hits, the number of hits expected 
##to occur in the gene set, the actual hits observed in the gene set, 
##and the pvalue from a hypergeometric test.

hyperGeoTest <- function(geneSet, universe, hits) {
	##number of genes in universe
	N <- length(universe) 			
	##remove genes from gene set that are not in universe			
	geneSet <- intersect(geneSet[[1]], universe) 
	##size of gene set	
	m <- length(geneSet) 							
	Nm <- N-m	
	##hits in gene set
	overlap <- intersect(geneSet, hits) 	
	##number of hits in gene set		
	k <- length(overlap) 							
	n <- length(hits)	
	HGTresults <- phyper(k, m, Nm, n, lower.tail = F)
	ex <- (n/N)*m
	if(m == 0) HGTresults <- NA
	hyp.vec <- c(N, m, n, ex, k, HGTresults)
	names(hyp.vec) <- c("Universe Size", "Gene Set Size", "Total Hits", 
		"Expected Hits", "Observed Hits", "Pvalue")
	return(hyp.vec)
}

