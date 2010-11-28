hyperGeoTest <-
function(geneSet,universe,hits) {

		N<-length(universe) 						#number of genes in universe
		geneSet<-intersect(geneSet[[1]],universe) 	#remove genes from gene set that are not in universe
		m<-length(geneSet) 							#size of gene set
		Nm<-N-m	
		overlap<-intersect(geneSet,hits) 			#hits in gene set
		k<-length(overlap) 							#number of hits in gene set
		n<-length(hits)	
		HGTresults<-phyper(k,m,Nm,n,lower.tail=F)
		ex<-(n/N)*m
		if (m==0) HGTresults<-NA
		hyp.vec<-c(N,m,n,ex,k,HGTresults)
		names(hyp.vec)<-c("Universe Size","Gene Set Size","Total Hits","Expected Hits","Observed Hits","Pvalue")
		return(hyp.vec)
}

