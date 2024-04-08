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
  subjectSeqDT=NULL,
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

  makeSfile <- F
  makeQfile <- F

  blTypeChar <- "n"
  blTypeShort <- "nucl"
  if(program %in% c("blastp","blastx")){
    blTypeChar <- "p"
    blTypeShort <- "prot"
  }

  if(is.null(subjectFname) & !is.null(subjectSeqDT)){
    is_seqDT(subjectSeqDT,objName=deparse(substitute(subjectSeqDT)),croak=T)
    makeSfile <- T
    writeFasta(subjectSeqDT,subjectFname<-tempfile(fileext = ".fasta"))
  }
  if(is.null(queryFname) & !is.null(queryFastaDT)){
    is_seqDT(queryFastaDT,objName=deparse(substitute(queryFastaDT)),croak=T)
    makeQfile <- T
    writeFasta(queryFastaDT,queryFname<-tempfile(fileext = ".fasta"))
  }

  l_ply(subjectFname,function(sFn){
    if(!blastDBexists(sFn,typePrefix=blTypeChar)){
      ce("Making nucleotide blastDB for ",sFn)
      mcmd <- paste0(makeBlastDbBinary," -dbtype ",blTypeShort," -in ",sFn)
      ce("\tRunning command: ",mcmd)
      system(mcmd)
    }
  })

  if(addAlignments==TRUE){
    outFmtArg %<>% paste(outFmtArg,c("sseq qseq"))
    outputColNames %<>% c("sAlnSeq","qAlnSeq")
    outputColClasses %<>% c("character","character")
  }



  bl <- ldtply(subjectFname,function(sFn){
    tmp <- ldtply(queryFname,function(qFn){
      bcmd <- paste0(blastBinaryDir,"/",program," -query ",qFn," -db ",sFn," -outfmt '",outFmtArg,"' ",moreBlastArgs)
      ce("Running command: ",bcmd)
      tmp <- fread(cmd = bcmd,col.names=outputColNames,colClasses=outputColClasses)
      if(!makeQfile){ tmp[,queryFname:=qFn][] }
      tmp
    })
    if(!makeSfile){ tmp[,subjectFname:=sFn][] }
    tmp
  })

  if(!is.null(outFastaFname)){
    write.table(bl,outFastaFname,row.names=F,sep="\t",quote=F)
  }

  if(makeSfile){unlink(Sys.glob(paste0(subjectFname,".*")))}
  if(makeQfile){unlink(subjectFname)}

  bl
}





