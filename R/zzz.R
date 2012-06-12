.onLoad<-function(libname, pkgname) {
	packageStartupMessage(
"/////////////////////////////////////////////////////////////////////////////\n
//--------------------This is HTSanalyzeR version 2.9.1--------------------//\n 
//------------please use function changes() to see new changes-------------//\n
/////////////////////////////////////////////////////////////////////////////\n", 
appendLF=FALSE)
}
changes<-function() {
cat(
"- dependent pakcages BioNet, cellHTS2, AnnotationDbi, biomaRt, RankProd were 
  moved from 'depend' to 'import' field in DESCRIPTION\n
- the problem of no global binding for 'org.Rn.egGO2EG' in function 'GOGeneSets' 
  and 'KeggGeneSets' was resolved \n
- default download link of the BioGRID database in function 'biogridDataDownload' 
  was updated to version 3.1.89 (tested on Jun 11, 2012)\n")
}
