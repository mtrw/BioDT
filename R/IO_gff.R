#' Read annotation files in gff format into a data.table
#'
#' @param faiFname A character vector. The .gff filename(s). [no default]
#' @returns A data.table with coordDT flavour.
#' @export
readGff <- function(gffFname) {
  ldtply(gffFname,function(fn){
    fread(fn,select=c(1,3,4,5,7,9),col.names=c("seqId","class","start","end","strand","attributes"))[,c("start","end"):=.(pmin(end,start),pmax(end,start))][,gffFname:=fn][]
  })
}
