#' Call blast* aligner family programs
#'
#' @param param1_name What it is. Defaults?
#' @param param2_name What it is. Defaults?
#' @returns What the function returns.
#' @description
#' Describe what the function does in detail.
#' @examples
#' Example 1 code
#'
#' Example 2 code
#' @export
blast <- function(
  subjectFname=NULL,
  queryFname=NULL,
  subjectFastaDT=NULL,
  queryFastaDT=NULL,
  outFastaFname=NULL,
  program="blastn",
  blastFmtArgs="",
  outFmtArg="6       saccver  qaccver slen qlen length sstart send qstart qend pident evalue bitscore",
  outputColNames=c( "sSeqId",    "qSeqId" , "sLength" ,  "qLength" , "matchLength" , "sStart" , "sEnd" ,   "qStart" ,  "qEnd" ,   "pctId" , "eValue" , "bitscore" ),
  outputColClasses=c("character", "character", "integer", "integer", "integer",      "numeric", "numeric", "numeric" , "numeric", "numeric","numeric", "numeric" ),
  addAlignments=FALSE,
  moreBlastArgs="",
  blastBinaryDir=sub("/blastn$","",system("which blastn",intern=T)),
  makeBlastDbBinary=system("which makeblastdb",intern=T)
){

  if(length(subjectFname)>1){
    stop("Woah there cowgirl! One subject file at a time please.")
  }

  cleanSfile <- F
  cleanQfile <- F

  blTypeChar <- "n"
  blTypeShort <- "nucl"
  if(program %in% c("blastp","blastx")){
    blTypeChar <- "p"
    blTypeShort <- "prot"
  }

  if(is.null(subjectFname) & !is.null(subjectFastaDT)){
    cleanSfile <- T
    writeFasta(subjectFastaDT,subjectFname<-tempfile(fileext = ".fasta"))
  }
  if(is.null(queryFname) & !is.null(queryFastaDT)){
    cleanQfile <- T
    writeFasta(queryFastaDT,queryFname<-tempfile(fileext = ".fasta"))
  }

  if(!blastDBexists(subjectFname,typePrefix=blTypeChar)){
    ce("Making nucleotide blastDB for ",subjectFname)
    mcmd <- paste0(makeBlastDbBinary," -dbtype ",blTypeShort," -in ",subjectFname)
    ce("\tRunning command: ",mcmd)
    system(mcmd)
  }

  if(addAlignments==TRUE){
    outFmtArg %<>% paste(outFmtArg,c("sseq qseq"))
    outputColNames %<>% c("sAlnSeq","qAlnSeq")
    outputColClasses %<>% c("character","character")
  }

  bcmd <- paste0(blastBinaryDir,"/",program," -query ",queryFname," -db ",subjectFname," -outfmt '",outFmtArg,"' ",moreBlastArgs)
  ce("Running command: ",bcmd)
  bl <- fread(cmd = bcmd,col.names=outputColNames,colClasses=outputColClasses)
  if(!is.null(outFastaFname)){
    write.table(bl,outFastaFname,row.names=F,sep="\t",quote=F)
  }

  if(cleanSfile){unlink(Sys.glob(paste0(subjectFname,".*")))}
  if(cleanQfile){unlink(subjectFname)}

  bl
}





