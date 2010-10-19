#This function gets rid of the duplicates in a gene list. It also orders the gene list by phenotype.
duplicateRemover <- function(geneList, method=c("max","min","average","fc.avg"), absValue=FALSE){
	
	#Check that geneList is a single vector
	if(!is.vector(geneList)) 
		stop("Your geneList should be a vector")
	#Check that the data vector is named 	
	if(is.null(names(geneList))) 
		stop("Your data vector is not named")
	#Check that the method argument is correctly specified
	if(!(length(method) == 1 && (method %in% c("max","min","average","fold.change.average"))))
		stop("The 'method' argument should be only one of the following character strings: 'max', 'min', 'average', 'fc.avg(fold change average)'")
	#Get the unique names and create a vector that will store the processed values corresponding to those names
	geneList.names<-names(geneList)
	datanames<-unique(geneList.names)
	data.processed<-rep(0,length(datanames))
	names(data.processed)<-datanames
	data.processed2<-data.processed
	l.data.processed<-length(data.processed)	
	#If the absolute value of the min is bigger than the absolute value of the max, then it is the min that is kept		
	if(method=="max") {
		data.processed<-sapply(datanames, 
				function(i) {
					this.range<-range(geneList[which(geneList.names==i)])
					this.range[which.max(abs(this.range))]
				}
		)
	} else if(method=="min") {
		data.processed<-sapply(datanames, 
				function(i) {
					this.range<-range(geneList[which(geneList.names==i)])
					this.range[which.min(abs(this.range))]
				}
		)
	} else if(method=="average") {
		data.processed<-sapply(datanames, 
				function(i) {
					mean(geneList[which(geneList.names==i)])
				}
		)
	} else if(method=="fc.avg") {
		neg.fcs<-which(geneList < 1) 
		geneListRatios<-geneList
		#convert from fold change to ratio		
		geneListRatios[neg.fcs]<-abs(1/geneList[neg.fcs])
		#average the values across replicates		
		data.processed<-sapply(datanames,function(i) {mean(geneListRatios[which(geneList.names==i)])})	
		#convert back to fold change				
		neg.fcs<-which(data.processed<1)
		data.processed[neg.fcs]<-(-1/data.processed[neg.fcs]) 
	}
	#order the final vector: considering absolute values or signed values and return the apporpriate vector	
	if(!absValue) {
		data.processed<-data.processed[order(data.processed,decreasing=TRUE)]
		return(data.processed)
	} else {
		data.processed.abs<-abs(data.processed)[order(abs(data.processed),decreasing=TRUE)]	
		return(data.processed.abs)
	}	
}

