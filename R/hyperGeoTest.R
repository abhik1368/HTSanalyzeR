hyperGeoTest <-
function(GeneSet,GeneList,hits) {
	
	N<-length(GeneList) ##number of genes in universe
	GeneSet<-intersect(GeneSet,GeneList) ##remove genes from gene set that are not in universe
	m<-length(GeneSet) ##size of gene set
	

	Nm<-N-m
	
	overlap<-intersect(GeneSet,hits) ##hits in gene set
	
	k<-length(overlap) ##number of hits in gene set
	n<-length(hits)
	
	HGTresults<-phyper(k,m,Nm,n,lower.tail=F)

	
	ex<-(n/N)*m
	
	if (m==0) HGTresults<-NA

	hyp.vec<-c(N,m,n,ex,k,HGTresults)
	names(hyp.vec)<-c("Universe Size","Gene Set Size","Total Hits","Expected Hits","Observed Hits","Pvalue")

	#return(list("Gene Set Size"=m,"Hits in Gene Set"=k,"p-value"=hyp.results))
	return(hyp.vec)

	
}

