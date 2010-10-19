###############################################################################
# TODO: 
# Xin Wang <xw264@cam.ac.uk>
# Advisor: Florian Markowetz <florian.markowetz@cancer.org.uk> 
# University of Cambridge Deparment of Oncology
# Cancer Research UK - Cambridge Research Institute
# At 15:22:08, on 9 Oct 2010
###############################################################################
networkPlot <-
function(
	nwAnalysisOutput,
	phenotypeVector=NA,
	filepath=".",
	filename="EnrichedSubNw.png"
) {
	subnw<-nwAnalysisOutput$subnw
	labels<-nwAnalysisOutput$labels
	#If no phenotype vector is specified, then we can just plot the module	
	if(length(phenotypeVector) == 1) {
		png("EnrichedSubNw.png",width = 900, height = 900)
		plotModule(subnw,labels=labels)
		dev.off()
	} else {
		#If a phenotype vector is specified, it will be used to color the nodes
		#check that the phenotypeVector has the right format (single numerical vector)
		if(!is.numeric(phenotypeVector) && !is.integer(phenotypeVector)) 
			stop("The phenotypeVector should be a single vector of numerical phenotypes")
		#Check that the phenotypeVector is a named vector
		if(!is.null(dim(phenotypeVector)) || is.null(names(phenotypeVector))) 
			stop("Please provide a phenotypeVector as a single named vector")	
		#this vector holds the phenotype for the nodes of the sub-network	
		diff.expr<-phenotypeVector[nodes(subnw)]
		names(diff.expr)<-nodes(subnw)
		#this vector contains the information of wether a node has an associated phenotype (1) or not (-1)
		#this information will be used to give a different shape to the nodes of the network
		present<-rep(1,length(nodes(subnw)))
		present[which(is.na(diff.expr))]<--1
		#this replaces all phenotypes of non-phenotyped nodes by a zero		
		diff.expr[which(is.na(diff.expr))]<-0
		names(present)<-nodes(subnw)
		#Plot the module	
		if(length(nodes(subnw)) == 1) {
			png(file.path(filepath,filename),width = 900, height = 900)
			plotModule(subnw,labels=labels,scores=present,diff.expr=diff.expr)
			dev.off()
		} else {
			Tcolors<-diff.expr
			Tcolors2<-diff.expr
			if(max(abs(Tcolors))<5)  
				Tcolors <- Tcolors*5
			# set red colors
			if(any(Tcolors>0)){
				maxRed <- max(ceiling(abs(Tcolors[which(Tcolors>0)])))
				redCols <- colorRampPalette(colors=c("white", "red"))
				redVec <- redCols(maxRed)
				Tcolors2[which(Tcolors>0)] <- redVec[ceiling(abs(Tcolors[which(Tcolors>0)]))]
			}
			#set the greens
			if(any(Tcolors<0)){
				maxGreen <- max(ceiling(abs(Tcolors[which(Tcolors<0)])))
				greenCols <- colorRampPalette(colors=c("white", "green"))
				greenVec <- greenCols(maxGreen)
				Tcolors2[which(Tcolors<0)] <- greenVec[ceiling(abs(Tcolors[which(Tcolors<0)]))]
			}	  	
			colScale<-unique(Tcolors2)
			colboundary<-rep(0,length(colScale))
			for(i in 1:length(colScale)) {
				values<-diff.expr[which(Tcolors2 == colScale[i])]
				colboundary[i]<-values[which(abs(values) == max(abs(values)))[1]]
			}
			colMatrix<-cbind(colboundary[order(colboundary)],colScale[order(colboundary)])
			png(file.path(filepath,filename), width = 900, height = 900)
			plotModule(subnw,labels=labels,scores=present,diff.expr=diff.expr)
			points(x=rep(-1.2,length(unique(Tcolors2))),y=seq(1.2,(1.2-(0.05*length(colMatrix[,2]))),length.out=length(colMatrix[,2])),pch=15,col=colMatrix[,2])
			text(x=rep(-1.1,length(unique(Tcolors2))),y=seq(1.2,(1.2-(0.05*length(colMatrix[,2]))),length.out=length(colMatrix[,2])),labels=signif(as.numeric(colMatrix[,1]),digits=2),cex=0.8)
			dev.off()
		}
	}
}