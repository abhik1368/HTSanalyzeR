#This function computes the overlap between a gene list and a gene set or a collection of gene sets.  
#This might be useful e.g. to check that the identifier mapping was done correctly.
overlapComputing <-
function(geneList,geneSet){
#Check that the geneList is a named vector
	if(is.null(dim(geneList)) != TRUE | is.null(names(geneList))) stop("Please provide a geneList as a single named vector")
#Check that the geneList does not have length zero
	if(length(geneList) == 0) stop("The geneList has length zero")
#Compute the overlap between the geneList names and the GeneSet/GeneSetCollection identifiers
	if(class(geneSet)[1]=="GeneSet"){
		overlap<-length(intersect(names(geneList),geneIds(geneSet)))
		uniqueIds<-length(geneIds(geneSet))
		}else{
			if(class(geneSet)[1]=="GeneSetCollection"){
				geneSetIncidence<-incidence(geneSet)
				overlap<-length(intersect(names(geneList),colnames(geneSetIncidence)))
				uniqueIds<-dim(geneSetIncidence)[2]
				}else{
#this computes the overlap	between the geneList names and any type of object that can be coerced by as.character				
					overlap<-length(intersect(names(geneList),as.character(geneSet)))
					uniqueIds<-length(geneSet)
					}
			}	
#this prints the information from this exploration on the screen			
	print(paste("Number of genes in your gene list: ", length(geneList)))
	print(paste("Number of genes in your gene set(s): ", uniqueIds))
	print(paste("Number of genes from your gene list matched to the gene set(s): ", overlap,
		"(",((overlap/length(geneList))*100),"%)"))
	}

