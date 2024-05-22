#' Read .fai files into a data.table
#'
#' @param faiFname A character vector. The .fai filename(s). [no default]
#' @returns A data.table with .fai flavour.
#' @export
readGff <- function( #read a gff file (untested)
  gffFname
){
  ldtply(gffFname,function(fn){
    fread(gffFnamesTable[genome==ldprCoords[i,genome],gffFname],select=c(1,3,4,5,7,9),col.names=c("seqId","class","start","end","strand","attributes"))[,c("start","end"):=.(end,start)][,gffFname:=fn][]
  })
}
