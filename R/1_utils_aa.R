lookupAA_CODON2symbol_functory <- function(codonTable){
  function(CODON){
    c <- data.table(CODON=CODON)
    codonTable[c,on="CODON"]$aaSymbol
  }
}
#' @export
lookupAA_CODON2symbol <- lookupAA_CODON2symbol_functory(codonTable)







lookupAA_codon2symbol_functory <- function(codonTable){
  function(codon){
    c <- data.table(codon=codon)
    codonTable[c,on="codon"]$aaSymbol
  }
}
#' @export
lookupAA_codon2symbol <- lookupAA_codon2symbol_functory(codonTable)







lookupAA_symbol2aaProperty_functory <- function(aaTable){
  function(aaSymbol){
    s <- data.table(aaSymbol=aaSymbol)
    aaTable[s,on="aaSymbol"]$aaProp
  }
}
#' @export
lookupAA_symbol2aaProperty <- lookupAA_symbol2aaProperty_functory(aaTable)


