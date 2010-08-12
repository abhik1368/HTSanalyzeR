#This function takes as input a matrix of p values for example obtained from a GSEA on multiple phenotypes,
#with a row for each gene set and a column for each phenotype
#and aggregates the p values by row (i.e. one aggregated p value for each gene set)
#according to Fisher or Stouffer's methods.
#It takes as input:
	#-pvalMatrix: a matrix of p values, with rows named according to the gene set 
#(rows= gene sets, columns=multiple p values to be aggregated for that gene set)
	#-method: "stouffers" or "fishers"
	#-pAdjustMethod: a method for adjustment for multiple hypothesis testing 
#of the aggregated p values (see p.adjust)
	#-order: TRUE or FALSE, if TRUE the results table is ordered 
#according to the aggregated p values
#It outputs:
	#a matrix with a row for each gene set and two columns: 
#"Aggregated.p.value","Adjusted.aggregated.p.value"
aggregatePvals<-function(pvalMatrix,method=c("stouffers","fishers"),pAdjustMethod="BH",order=TRUE){
#Check that the argument pvalMatrix has the right format
	if(class(pvalMatrix) != "matrix") stop("The argument pvalMatrix should be a matrix")
	if(is.null(rownames(pvalMatrix))) stop("The argument pvalMatrix should be a matrix with rownames")
#Check that the argument method has the right format
	if(is.na(match(method,c("stouffers","fishers")))) stop('The argument method should be one of "stouffers" or "fishers"')
	if(length(match(method,c("stouffers","fishers"))) != 1) stop('The argument method should be one of "stouffers" or "fishers"')
#Check that the argument method has the right format
	if(is.na(match(pAdjustMethod,c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")))) stop('The argument method should be one of "stouffers" or "fishers"')
	if(length(match(pAdjustMethod,c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"))) != 1) stop('The argument method should be one of "stouffers" or "fishers"')	
#Perform the aggregation
	aggr.pval<-rep(1,dim(pvalMatrix)[1])
	names(aggr.pval)<-rownames(pvalMatrix)
	if (method=="fishers") {
		pvalMatrixLogged<-log(pvalMatrix)
		aggr.pval<-(-2)*apply(pvalMatrixLogged,1,sum)		
		aggr.pval<-pchisq(aggr.pval,df=2*dim(pvalMatrix)[2],lower.tail=FALSE)
	}	
	if (method=="stouffers") {
		aggr.pval<-apply(pvalMatrix,1, function(x) {
			newPval<-pnorm(sum(qnorm(x))/sqrt(length(x)))
			newPval
			})
		}
	adj.aggr.pval<-p.adjust(aggr.pval,method=pAdjustMethod)	
	aggr.pval.table<-cbind(aggr.pval,adj.aggr.pval)
	rownames(aggr.pval.table)<-names(aggr.pval)
	colnames(aggr.pval.table)<-c("Aggregated.p.value","Adjusted.aggregated.p.value")
	if(order == TRUE){
		aggr.pval.table<-aggr.pval.table[order(aggr.pval.table[,"Aggregated.p.value"]),]
		}
}	