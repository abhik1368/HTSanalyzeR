###############################################################################
# generic functions
# Xin Wang <xw264@cam.ac.uk>
# Advisor: Florian Markowetz <florian.markowetz@cancer.org.uk> 
# University of Cambridge Deparment of Oncology
# Cancer Research UK - Cambridge Research Institute
# At 22:26:44, on 17 Nov 2010
###############################################################################
setGeneric("preprocess",function(object,...) standardGeneric("preprocess"), package="HTSanalyzeR")
setGeneric("score",function(object,...) standardGeneric("score"), package="HTSanalyzeR")
setGeneric("selectHits",function(object,...) standardGeneric("selectHits"), package="HTSanalyzeR")
setGeneric("calculateNWAPvals",function(object,...) standardGeneric("calculateNWAPvals"), package="HTSanalyzeR")
setGeneric("interactome",function(object,...) standardGeneric("interactome"), package="HTSanalyzeR")
setGeneric("analyze",function(object,...) standardGeneric("analyze"), package="HTSanalyzeR")
setGeneric("report",function(object,...) standardGeneric("report"), package="HTSanalyzeR")
setGeneric("reportAll",function(gsca,nwa,...) standardGeneric("reportAll"), package="HTSanalyzeR")
setGeneric("plotSubNet",function(object,...) standardGeneric("plotSubNet"), package="HTSanalyzeR")
setGeneric("writeHits",function(object,...) standardGeneric("writeHits"), package="HTSanalyzeR")
setGeneric("plotGSEA",function(object,...) standardGeneric("plotGSEA"), package="HTSanalyzeR")
setGeneric("getTopGeneSets",function(object,...) standardGeneric("getTopGeneSets"), package="HTSanalyzeR")
