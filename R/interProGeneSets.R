##	Returns a gene set collection where gene sets are Interpro IDs/descriptions.
##	The gene sets are limited to those with interproIDs with genes in the GeneList	
##	Inputs:
##	-GeneList: a numeric vector named with Entrez Gene IDs
##	-species
##	
##	Output:
##	-a GeneSetCollection object
InterProGeneSets<-function(GeneList, species=c("Hs","Mm","Rn","Dm","Ce")) {		
	ensembl<-useMart("ensembl")
	if (species=="Hs") x="hsapiens_gene_ensembl"
	if (species=="Mm") x="mmusculus_gene_ensembl"
	if (species=="Rn") x="rnorvegicus_gene_ensembl"
	if (species=="Dm") x="dmelanogaster_gene_ensembl"
	if (species=="Ce") x="celegans_gene_ensembl"
	
	ensembl<-useDataset(x,ensembl)
	
	#Gets interpro signatures for genes in GeneList
	ipro=getBM(attributes=c("interpro","interpro_description","entrezgene"),
		filters="entrezgene",values = names(GeneList), mart = ensembl)
	int.pro<-ipro[order(ipro[,1]),]
	#Identify genes which do not map to Interpro
	noInterproInfo<-which(int.pro[,1]=="")
	cat(paste(length(noInterproInfo), "genes did not map to any Interpro signatures \n", sep=""))
	intPro<-int.pro[-noInterproInfo,]
	#Create list of all genes associated with each interpro ID
	intProNames<-unique(intPro[,1])
	intProList<-vector("list",length=length(intProNames))
	for(i in 1:length(intProNames)) {
		ind<-which(intPro[,1]==intProNames[i])
		geneIds<-as.character(intPro[ind,3])
		intProList[[i]]<-geneIds
	}
	#Add interpro description to ID in names
	fullIntProNames<-paste(intPro[,1],"::",intPro[,2])
	fullIntProNames<-unique(fullIntProNames)
	names(intProList)<-fullIntProNames
	#Create GeneSetCollection
	#gsc <- GeneSetCollection(mapply(function(geneIds, InterProId) {
	#	GeneSet(as.character(geneIds),geneIdType=EntrezIdentifier(),setName=InterProId)
	#	}, intProList, names(intProList)))
	return(intProList)
}
