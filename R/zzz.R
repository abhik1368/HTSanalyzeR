.onLoad<-function(libname, pkgname) {
	packageStartupMessage(
"/////////////////////////////////////////////////////////////////////////////\n
//------------------Thanks for using HTSanalyzeR v 2.9.3-------------------//\n 
//------------please use function changes() to see new changes-------------//\n
//------------please report any bug to xin.wang@cancer.org.uk--------------//\n
/////////////////////////////////////////////////////////////////////////////\n", 
appendLF=FALSE)
}
changes<-function() {
cat(
"//changes in v2.9.3\n
- The dependency of igraph was changed to igraph0 to adapt to the change in 
  BioNet.\n
//changes in v2.9.2\n
- Gene set overrepresentation and enrichment analysis can run independently 
  using the S4 method 'analyze' by specifying argument 'doGSOA' (for 
  hypergeometric test based overrepresentation analysis) and 'doGSEA' (for 
  GSEA). More details in ?analyze.
- One bug corrected for function 'analyzeGeneSetCollection' when calculating the 
  overlap of significant gene sets in both GSEA and HyperGeo.\n
//changes in v2.9.1\n
- dependent pakcages BioNet, cellHTS2, AnnotationDbi, biomaRt, RankProd were 
  moved from 'depend' to 'import' field in DESCRIPTION
- the problem of no global binding for 'org.Rn.egGO2EG' in function 'GOGeneSets' 
  and 'KeggGeneSets' was resolved
- default download link of the BioGRID database in function 'biogridDataDownload' 
  was updated to version 3.1.89 (tested on Jun 11, 2012)\n")
}
