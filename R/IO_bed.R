#' Read a bed file
#'
#' @param bedFname A character vector. The .bed filenames to read. Calls awk. [no default]
#' @param saveFnames Logical value. Include a bedFname column in the output? [TRUE]
#' @returns A data.table with .bed flavour.
#' @export
readBed<- function(
    bedFname,
    saveFnames=TRUE
){
  colNames <- c("seqId","start","end","name","score","strand","seqId", "start", "end", "name", "score", "strand","thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts", "bedFname")
  colTypes <- c("character", "numeric", "numeric", "character", "numeric", "character", "numeric", "numeric", "character", "integer", "character", "character", "character")
  out <- ldply(bedFname,function(fn){
    o <- fread(fn,header=F)
    colnames(o)<-colNames[1:ncol(o)]
    for(c in 1:ncol(o)){
      cat("Change to ",colTypes[c],"\n")
      o[,eval(colNames[c]):=as(get(colNames[c]),colTypes[c])]
      o[,eval(colNames[c]):=as(get(colNames[c]),colTypes[c])]
      print(str(o))
    }
    if(saveFnames){ o[,bedFname:=fn] }
    o
  }) %>% setDT
  if(!isCoordDT(out,error=F)){ stop("Resultant table does not pass validation. Check the input as at least three columns (seqId, start, end) with sensible entries.") }
  return(out)
}
# fn <- bedFname <- bed6Fname
# readBed(bed6Fname) %>% isBedDT


#' Write a bed file
#' @param bed A .fasta-flavoured data.table.
#' @param bedFname A character vector. The filename to write. If NULL will use the deparsed/substituted variable name given as bed [NULL]
#' @returns NULL
#' @export
writeBed <- function(
  bedDT,
  bedFname = NULL
){
  bedDT <- copy(bedDT)
  isCoordDT(bedDT)
  bedDT[,end:=end+1L] # To "trad" bed coords
  if(is.null(bedFname)){bedFname <- paste0(deparse(substitute(fasta)),".bed")}
  #CHOOSE ALL BED-VALID ROWS IN ORDER TO DO
  colNames <- c(  "seqId",     "start",   "end",     "name",      "score",   "strand",    "thickStart", "thickEnd", "itemRgb",   "blockCount", "blockSizes", "blockStarts", "bedFname")
  colClasses <- c("character", "numeric", "numeric", "character", "numeric", "character", "numeric",    "numeric",  "character", "integer",    "character",  "character",   "character")
  chooseCol_i <- 0
  for(i in seq_along(colNames)){
    if(!colNames[i]%in%colnames(bedDT)){ break }
    if(!colClasses[i]==class(bedDT[,get(colNames[i])])){ warning(paste0("Column ",colNames[i]," of bed must be type \"",colClasses[i],"\"; This and subsequent columns will not be included ...")); break; }
    chooseCol_i <- i
  }
  getColNames <- colNames[1:chooseCol_i]
  ce("Attempting to include columns: ")
  ce(getColNames)
  write.table(bedDT[,..getColNames],file=bedFname,col.names=F,row.names=F,sep="\t",quote=F)
}
