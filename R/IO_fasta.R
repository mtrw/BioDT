#' Read a fasta or fasta.gz file
#'
#' @param fastaFname A character vector. The .fasta filenames to read. Calls awk. [no default]
#' @param saveFnames Logical value. Include a fastaFname column in the output? [TRUE]
#' @param namesOnly Logical value. If true, will return only the seqName column. Calls samtools faidx.
#' @param awkBin Character value. Where to find awk [system("which awk", intern=TRUE)]
#' @param samtoolsBin Character value. In case namesOnly==TRUE, where to find samtools. [system("which samtools", intern=TRUE)]
#' @param gzipBin Character value. Where to find gzip. [system("which gzip", intern=TRUE)]
#' @returns A data.table with .fasta flavour.
#' @description
#' Calls awk normally, and samtools faidx is namesOnly==TRUE. Calls gzip if the filename ends in ".gz".
#' @export
readFasta <- function(
    fastaFname,  #fastaFname=c("/data/gpfs/projects/punim1869/shared_data/misc_sequence/btr_btrLike_queries_Morex_GP.fasta","/data/gpfs/projects/punim1869/shared_data/misc_sequence/btr_btrLike_queries_Morex_GP.fasta"); fn <- fastaFname[1]
    saveFnames=TRUE,
    namesOnly=FALSE,
    awkBin=system("which awk", intern=TRUE),
    samtoolsBin=system("which samtools", intern=TRUE)
){
  if(namesOnly){
    getFai(fastaFname)[,.(seqName)]
  } else {
    tf <- tempfile(fileext=".tsv")
    out <- ldply(fastaFname,function(fn){
      catBashArg <- if(grepl("\\.gz^",fn)){ "gunzip -c" } else { "cat" }

      awkScript <- '{sub(/\\r/,"")} /^>/ { if(NR>1){printf("\\n")}; printf substr($1,2,length($0))" "; next; } 1 {printf $1}'
      command <- paste0(catBashArg," ",fn," | ", awkBin , " '" , awkScript , "' ")
      #system(command,intern = T)
      o <- fread(cmd=command,header=F,col.names=c("seqId","seq"))
      if(saveFnames){ o[,fastaFname:=fn] }
      o
    }) %>% setDT
    unlink(tf)
    return(out)
  }
}




#' Write a fasta or fasta.gz file
#' @param fasta A .fasta-flavoured data.table.
#' @param fastaFname A character vector. The .fasta filenames to write. If NULL will use the deparsed/substituted variable name given as the arg to `fasta` [NULL]
#' @param gzip Logical value. Create gzipped outpit? Calls gzip via R's native gzfile() interface. [FALSE]
#' @param nameCol Character value. The column to take fasta header names from. ["seqId"]
#' @returns NULL
#' @description
#' Calls awk normally, and samtools faidx is namesOnly==TRUE. Calls gzip if the filename ends in ".gz". Note this writes all the given fasta-dts into ONE file, i.e., if you writeFasta(readFasta(someFastaFileNames)) you will concatenate all the someFastaFileNames into one output.
#' @export
writeFasta <- function(
  fasta,
  fastaFname = NULL,
  gzip = F,
  nameCol = "seqId"
){
  is_seqDT(fasta)
  gzExt <-   if(gzip==TRUE){".gz"}else{""}
  if(is.null(fastaFname)){fastaFname <- paste0(deparse(substitute(fasta)),".fasta",gzExt)}
  fileCon <- if(gzip==TRUE){gzfile(fastaFname,open="w")}else{file(fastaFname,open="w")}
  l_ply(1:nrow(fasta),function(i){
    writeLines(paste0(">",fasta[i,get(nameCol)]),fileCon)
    writeLines(paste0(fasta[i,seq]),fileCon)
  })
  close.connection(fileCon)
}
