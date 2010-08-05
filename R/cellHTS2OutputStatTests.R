cellHTS2OutputStatTests <-
function(cellHTSobject,annotationColumn="GeneID",controls="neg",alternative=c("two.sided","less","greater"),logged=FALSE,tests=c("T-test","MannWhitney","RankProduct")){
#check that 'cellHTSobject' is a cellHTS cellHTSobject
	if(class(cellHTSobject) != "cellHTS") stop("The argument 'cellHTSobject' is not a cellHTS object")
#check that the cellHTS2 cellHTSobject is in the right format 
#(configured, so that we know which rows are samples and which are controls
#normalized and annotated, so that we can match the names of the constructs)
	if(state(cellHTSobject)["configured"] == FALSE) stop("The cellHTS object should be configured to perform the statistical tests")
	if(state(cellHTSobject)["normalized"] == FALSE) warning("Your cellHTS object has not been normalized, this could impact the results of these tests",immediate.=TRUE)
	if(state(cellHTSobject)["scored"] == TRUE) warning("This cellHTS object has been scored and the statistical analysis will be perfomed on the score rather than the normalized signal intensities",immediate.=TRUE)
	if(state(cellHTSobject)["annotated"] == FALSE) warning("This cellHTS object has not been annotated",immediate.=TRUE)
#check that the annotationColumn has been specified as a single character string
	if(is.character(annotationColumn) != TRUE | length(annotationColumn) !=1 ) stop("The 'annotationColumn' parameter does not have the right format")
#check that the annotationColumn matches a column in the fData(cellHTSobject) dataframe	
	if(is.null(match(annotationColumn, colnames(fData(cellHTSobject))))) stop("The 'annotationColumn' parameter does not match to any column in your cellHTS object")
#check that the 'controls' parameter has been specified as a single character string
	if(is.character(controls) != TRUE | length(controls) !=1 ) stop("The 'controls' parameter does not have the right format")
#check that the  'controls' parameter matches a status in the 'controlStatus' column of the cellHTS cellHTSobject	
	if(is.null(match(controls, fData(cellHTSobject)[,"controlStatus"]))) stop("The 'controls' parameter does not match to any status in the 'controlStatus' column of your cellHTS object")	
#check that the 'alternative' parameter has been correctly set 
	if(length(match(alternative,c("two.sided","less","greater"))) != 1) stop("The 'alternative' parameter has not been set properly, please provide only one of the following charcater strings: 'two.sided','less','greater'")
#make a named data matrix (only samples) rows=features, columns=replicates, with row names = identifiers in the "annotationColumn" of the fData() data frame
	dataNw<-Data(cellHTSobject)[,1:dim(Data(cellHTSobject))[2],1]
	rownames(dataNw)<-fData(cellHTSobject)[,annotationColumn] 
	dataNw<-dataNw[-which(fData(cellHTSobject)[,"controlStatus"] != "sample"),]
	dataNw<-dataNw[-which(is.na(rownames(dataNw))==TRUE),]
#make a vector of data for the control 	population
	controlData<-Data(cellHTSobject)[which(fData(cellHTSobject)[,"controlStatus"]== controls),1:dim(Data(cellHTSobject))[2],1]
	controlData<-as.vector(controlData)
#compute the median of all samples, for the one sample tests
	mu=median(as.vector(dataNw),na.rm=TRUE)
#make a list of the data (one entry per unique ID): each entry in the list correspond to a unique name, 
#and the element under that entry is a vector of data of replicates for that unique construct
#formatting this as a list allows us to have different number of replicates for each construct 
	replicatesNames<-unique(rownames(dataNw))
	replicates<-as.list(rep(0,length(replicatesNames)))
	names(replicates)<-replicatesNames
	nreplicates<-length(replicates)
	for(i in 1:nreplicates){
		replicates[[i]]<-c(dataNw[which(rownames(dataNw)==names(replicates)[i]),])
		}
if(!is.na(match("T-test",tests))){		
#Compute the one sample t-test (only possible for those entries of the list that contain more than one replicate measurement
#otherwise the pvalue will be left at the default value of 1)		
	t.test.pvalues.one.sample<-rep(1,nreplicates)
	names(t.test.pvalues.one.sample)<-names(replicates)
	for(i in 1:nreplicates){
		if(sum(!is.na(replicates[[i]])) >= 2) t.test.pvalues.one.sample[i]<-t.test(x=replicates[[i]],mu=mu,alternative=alternative)$p.value
		}	
#Compute the two samples t-test (only possible for those entries of the list that contain more than one replicate measurement
#otherwise the pvalue will be left at the default value of 1)	
	t.test.pvalues.two.samples<-rep(1,nreplicates)
	names(t.test.pvalues.two.samples)<-names(replicates)
	for(i in 1:nreplicates){
		if(sum(!is.na(replicates[[i]])) >= 2) t.test.pvalues.two.samples[i]<-t.test(x=replicates[[i]],y=controlData,alternative=alternative)$p.value
		}	
	}
if(!is.na(match("MannWhitney",tests))){		
#Compute the one sample mann-whitney test(only possible for those entries of the list that contain more than one replicate measurement
#otherwise the pvalue will be left at the default value of 1)	
	mannW.test.pvalues.one.sample<-rep(1,nreplicates)
	names(mannW.test.pvalues.one.sample)<-names(replicates)
	for(i in 1:nreplicates){
		if(sum(!is.na(replicates[[i]])) >= 2) mannW.test.pvalues.one.sample[i]<-wilcox.test(x=replicates[[i]],mu=mu,alternative=alternative)$p.value
		}
#Compute the two samples mann-whitney test(only possible for those entries of the list that contain more than one replicate measurement
#otherwise the pvalue will be left at the default value of 1)	
	mannW.test.pvalues.two.samples<-rep(1,nreplicates)
	names(mannW.test.pvalues.two.samples)<-names(replicates)
	for(i in 1:nreplicates){
		if(sum(!is.na(replicates[[i]])) >= 2) mannW.test.pvalues.two.samples[i]<-wilcox.test(x=replicates[[i]],y=controlData,alternative=alternative)$p.value
		}
	}
if(!is.na(match("RankProduct",tests))){		
#Prepare the data for the Rank Product test: the function 'RP' requires as input a matrix with a row for each construct and a column for each replicate 
#(this function was built for microarrays, where each column could correspond to an array with a different treatment class
#this is not the case here, hence we set the class argument of the RP function to 1 for all columns
#Since our data might include varying number of replicates, a matrix of maximal dimensions (number of columns=max number of replicates) will be built
#with NAs when necessary		
	lengthreplicates<-sapply(replicates,length)
	maxlength<-max(lengthreplicates)
	replicatesmatrix<-c(replicates[[1]],rep(NA,(maxlength-length(replicates[[1]]))))
	for(i in 2:length(replicates)){
		if(length(replicates[[i]]) < maxlength)   replicatesmatrix<-rbind(replicatesmatrix,c(replicates[[i]],rep(NA,(maxlength-length(replicates[[i]])))))
		if(length(replicates[[i]]) == maxlength)   replicatesmatrix<-rbind(replicatesmatrix,c(replicates[[i]]))
		}
	rownames(replicatesmatrix)<-names(replicates)	
#Compute the Rank Product test
	rankptest<-RP(data=replicatesmatrix,cl=rep(1,ncol(replicatesmatrix)),logged=logged,gene.names=rownames(replicatesmatrix))	
	}
#Assemble the results as a column for each test and a row for each construct: if all 3 tests are performed: 
#the RP test produces one column for up-regulated genes and one for down-regulated ones (which constrasts with the
#other two types of tests that produce only one result per alternative, therefore, the RP column that will be in the output
#is different depending on the alternative chosen, which is why there are 3 parts to this assembly)
if(length(tests) == 3){
	stats<-cbind(t.test.pvalues.one.sample,t.test.pvalues.two.samples,mannW.test.pvalues.one.sample,mannW.test.pvalues.two.samples)
	rownames(stats)<-names(replicates)
	if(alternative == "two.sided"){
		upRP<-rankptest$pfp[,1]
		names(upRP)<-rownames(replicatesmatrix)
		downRP<-rankptest$pfp[,2]
		names(downRP)<-rownames(replicatesmatrix)
		stats<-cbind(stats,upRP,downRP)
		colnames(stats)<-c("t.test.pvalues.one.sample","t.test.pvalues.two.samples","mannW.test.pvalues.one.sample","mannW.test.pvalues.two.samples","rank.product.pfp.greater","rank.product.pfp.less")
		}
	if(alternative == "greater"){
		upRP<-rankptest$pfp[,1]
		names(upRP)<-rownames(replicatesmatrix)
		stats<-cbind(stats,upRP)
		colnames(stats)<-c("t.test.pvalues.one.sample","t.test.pvalues.two.samples","mannW.test.pvalues.one.sample","mannW.test.pvalues.two.samples","rank.product.pfp.greater")
		}
	if(alternative == "less"){
		downRP<-rankptest$pfp[,2]
		names(downRP)<-rownames(replicatesmatrix)
		stats<-cbind(stats,downRP)
		colnames(stats)<-c("t.test.pvalues.one.sample","t.test.pvalues.two.samples","mannW.test.pvalues.one.sample","mannW.test.pvalues.two.samples","rank.product.pfp.less")
		}
	}
#Assemble the results as a column for each test and a row for each construct: if only 1 test is performed: 	
if(length(tests) == 1){
	if(tests == "T-test") {
		stats<-cbind(t.test.pvalues.one.sample,t.test.pvalues.two.samples)
		rownames(stats)<-names(replicates)
		colnames(stats)<-c("t.test.pvalues.one.sample","t.test.pvalues.two.samples")
		}
	if(tests == "MannWhitney"){
		stats<-cbind(mannW.test.pvalues.one.sample,mannW.test.pvalues.two.samples)
		rownames(stats)<-names(replicates)
		colnames(stats)<-c("mannW.test.pvalues.one.sample","mannW.test.pvalues.two.samples")
		}
	if(tests == "RankProduct"){
		if(alternative == "two.sided"){
			upRP<-rankptest$pfp[,1]
			names(upRP)<-rownames(replicatesmatrix)
			downRP<-rankptest$pfp[,2]
			names(downRP)<-rownames(replicatesmatrix)
			stats<-cbind(upRP,downRP)
			rownames(stats)<-names(replicates)
			colnames(stats)<-c("rank.product.pfp.greater","rank.product.pfp.less")
			}
		if(alternative == "greater"){
			upRP<-rankptest$pfp[,1]
			names(upRP)<-rownames(replicatesmatrix)
			stats<-as.matrix(upRP,ncol=1)
			rownames(stats)<-names(replicates)
			colnames(stats)<-c("rank.product.pfp.greater")
			}
		if(alternative == "less"){
			downRP<-rankptest$pfp[,2]
			names(downRP)<-rownames(replicatesmatrix)
			stats<-as.matrix(downRP,ncol=1)
			rownames(stats)<-names(replicates)
			colnames(stats)<-c("rank.product.pfp.less")
			}
		}
	}
#Assemble the results as a column for each test and a row for each construct: if 2 test are performed (one block for each combination of tests):	
if(length(tests) == 2){	
		if(!is.na(match("T-test",tests)) & !is.na(match("MannWhitney",tests))){
			stats<-cbind(t.test.pvalues.one.sample,t.test.pvalues.two.samples,mannW.test.pvalues.one.sample,mannW.test.pvalues.two.samples)
			rownames(stats)<-names(replicates)
			colnames(stats)<-c("t.test.pvalues.one.sample","t.test.pvalues.two.samples","mannW.test.pvalues.one.sample","mannW.test.pvalues.two.samples")
			}
		if(!is.na(match("T-test",tests)) & !is.na(match("RankProduct",tests))){
			stats<-cbind(t.test.pvalues.one.sample,t.test.pvalues.two.samples)
			rownames(stats)<-names(replicates)
			if(alternative == "two.sided"){
				upRP<-rankptest$pfp[,1]
				names(upRP)<-rownames(replicatesmatrix)
				downRP<-rankptest$pfp[,2]
				names(downRP)<-rownames(replicatesmatrix)
				stats<-cbind(stats,upRP,downRP)
				colnames(stats)<-c("t.test.pvalues.one.sample","t.test.pvalues.two.samples","rank.product.pfp.greater","rank.product.pfp.less")
				}
			if(alternative == "greater"){
				upRP<-rankptest$pfp[,1]
				names(upRP)<-rownames(replicatesmatrix)
				stats<-cbind(stats,upRP)
				colnames(stats)<-c("t.test.pvalues.one.sample","t.test.pvalues.two.samples","rank.product.pfp.greater")
				}
			if(alternative == "less"){
				downRP<-rankptest$pfp[,2]
				names(downRP)<-rownames(replicatesmatrix)
				stats<-cbind(stats,downRP)
				colnames(stats)<-c("t.test.pvalues.one.sample","t.test.pvalues.two.samples","rank.product.pfp.less")
				}
			}
		if(!is.na(match("MannWhitney",tests)) & !is.na(match("RankProduct",tests))){
			stats<-cbind(mannW.test.pvalues.one.sample,mannW.test.pvalues.two.samples)
			rownames(stats)<-names(replicates)
			if(alternative == "two.sided"){
				upRP<-rankptest$pfp[,1]
				names(upRP)<-rownames(replicatesmatrix)
				downRP<-rankptest$pfp[,2]
				names(downRP)<-rownames(replicatesmatrix)
				stats<-cbind(stats,upRP,downRP)
				colnames(stats)<-c("mannW.test.pvalues.one.sample","mannW.test.pvalues.two.samples","rank.product.pfp.greater","rank.product.pfp.less")
				}
			if(alternative == "greater"){
				upRP<-rankptest$pfp[,1]
				names(upRP)<-rownames(replicatesmatrix)
				stats<-cbind(stats,upRP)
				colnames(stats)<-c("mannW.test.pvalues.one.sample","mannW.test.pvalues.two.samples","rank.product.pfp.greater")
				}
			if(alternative == "less"){
				downRP<-rankptest$pfp[,2]
				names(downRP)<-rownames(replicatesmatrix)
				stats<-cbind(stats,downRP)
				colnames(stats)<-c("mannW.test.pvalues.one.sample","mannW.test.pvalues.two.samples","rank.product.pfp.less")
				}				
			}
	}
return(stats)
	}

