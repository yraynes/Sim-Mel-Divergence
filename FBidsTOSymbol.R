convertFBids <- function(FBids){
  require(drosophila2.db)
  require(AnnotationDbi)
  genes <- AnnotationDbi::select(drosophila2.db, keys = FBids, columns=c("SYMBOL"), keytype="FLYBASE")
  return (genes$SYMBOL)
}
