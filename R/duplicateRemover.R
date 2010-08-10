#This function gets rid of the duplicates in a gene list. It also orders the gene list by phenotype.
duplicateRemover <-function(geneList,method=c("max","min","average","fold.change.average"),absValue=FALSE){
#Check that geneList is a single vector
	if(is.null(dim(geneList)) != TRUE) stop("Your geneList should be a vector")
#Check that the data vector is named 	
	if(is.null(names(geneList))) stop("Your data vector is not named")
#Check that the method argument is correctly specified
	if(length(method) != 1 | length(match(method,c("max","min","average","fold.change.average"))) != 1 | is.na(match(method,c("max","min","average","fold.change.average")))){
		stop("The 'method' argument should be only one of the following character strings: 'max', 'min', 'average', 'fold.change.average'")
		}
#Get the unique names and create a vector that will store the processed values corresponding to those names
	datanames<-unique(names(geneList))
	data.processed<-rep(0,length(datanames))
	names(data.processed)<-datanames
	l.data.processed<-length(data.processed)	
#Fill the new data vector, by element of that vector (matching the element by names): if method is 'max'
	if(method=="max"){
		for(i in 1:l.data.processed){
#Compute the min and max (signed) values across observations with the same name		
			Max<-max(geneList[which(names(geneList)==names(data.processed)[i])])
			Min<-min(geneList[which(names(geneList)==names(data.processed)[i])])
#If the absolute value of the min is bigger than the absolute value of the max, then it is the min that is kept		
			if(abs(Min) > abs(Max)){
				data.processed[i]<-Min
				}else{
					data.processed[i]<-Max
					}
			}
	}
#Fill the new data vector, by element of that vector (matching the element by names): if method is 'min'
	if(method =="min"){
		for(i in 1:l.data.processed){
#Compute the min and max (signed) values across observations with the same name			
			Max<-max(geneList[which(names(geneList)==names(data.processed)[i])])
			Min<-min(geneList[which(names(geneList)==names(data.processed)[i])])
#If the absolute value of the min is bigger than the absolute value of the max, then it is the max that is kept				
			if(abs(Min) > abs(Max)){
				data.processed[i]<-Max
				}else{
					data.processed[i]<-Min
					}
			}
	}
#Fill the new data vector, by element of that vector (matching the element by names): if method is 'average'	
	if(method == "average"){
		for(i in 1:l.data.processed){
			data.processed[i]<-mean(geneList[which(names(geneList)==names(data.processed)[i])])
			}
	}	
#Fill the new data vector, by element of that vector (matching the element by names): if method is 'fold.change.average'			
	if(method == "fold.change.average"){
#get the indices of negative fold changes	
		neg.fcs<-which(geneList < 1) 
		geneListRatios<-geneList
#convert from fold change to ratio		
		geneListRatios[neg.fcs]<-abs(1/geneList[neg.fcs])
#average the values across replicates			
		for(i in 1:l.data.processed){	
			data.processed[i]<-mean(geneListRatios[which(names(geneList)==names(data.processed)[i])])
			}
#convert back to fold change				
		neg.fcs<-which(data.processed<1)
		data.processed[neg.fcs]<-(-1/data.processed[neg.fcs]) 
	}	
#order the final vector: considering absolute values or signed values			
	data.processed.abs<-abs(data.processed)[order(abs(data.processed),decreasing=TRUE)]	
	data.processed<-data.processed[order(data.processed,decreasing=TRUE)]
#return the apporpriate vector	
	if (absValue==FALSE) {
		return(data.processed)
		}else{
			return(data.processed.abs)
			}
	}

