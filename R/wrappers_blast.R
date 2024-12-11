
#' A convenience function to read in tabular BLAST output in the format requested (under the bonnet) by this wrapper
#'
#' @export
readBlastTable <- function(fileName,outputColNames=c( "sSeqId",    "qSeqId" ,    "matchLength" , "sStart" , "qStart", "sEnd"  ,     "qEnd" ,   "pctId_noGaps" , "eValue" , "bitscore" ),outputColClasses=c("character", "character", "integer",      "numeric",  "numeric", "numeric" , "numeric", "numeric",       "numeric", "numeric" )){
  ldtply(fileName,function(fn){ fread(fn,col.names = outputColNames, colClasses = outputColClasses)[,blastTableFname:=fn] })
}



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
  querySeqDT=NULL,
  outFastaFname=NULL,
  program="blastn",
  addAlignments=FALSE,
  outFmtArg="6       saccver  qaccver length sstart qstart send qend pident evalue bitscore",
  outputColNames=c( "sSeqId",    "qSeqId" ,    "matchLength" , "sStart" , "qStart", "sEnd"  ,     "qEnd" ,   "pctId_noGaps" , "eValue" , "bitscore" ),
  outputColClasses=c("character", "character", "integer",      "numeric",  "numeric", "numeric" , "numeric", "numeric",       "numeric", "numeric" ),
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

  if(argNotGiven(subjectFname) & argGiven(subjectSeqDT)){
    is_seqDT(subjectSeqDT,objName=deparse(substitute(subjectSeqDT)),croak=T)
    makeSfile <- T
    writeFasta(subjectSeqDT,subjectFname<-tempfile(fileext = ".fasta"))
  }
  if(argNotGiven(queryFname) & argGiven(querySeqDT)){
    is_seqDT(querySeqDT,objName=deparse(substitute(querySeqDT)),croak=T)
    makeQfile <- T
    writeFasta(querySeqDT,queryFname<-tempfile(fileext = ".fasta"))
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
    outFmtArg %<>% paste(c("sseq qseq"))
    outputColNames %<>% c("sAlnSeq","qAlnSeq")
    outputColClasses %<>% c("character","character")
  }



  bl <- ldtply(subjectFname,function(sFn){
    tmp <- ldtply(queryFname,function(qFn){
      #dev sFn<-subjectFname[1]; qFn<-queryFname[1]
      oFile <- tempfile()
      bcmd <- paste0(blastBinaryDir,"/",program," -query ",qFn," -db ",sFn," -outfmt '",outFmtArg,"' ",moreBlastArgs," > ",oFile)
      ce("Running command: ",bcmd)
      system(bcmd) -> exitCode
      if(exitCode!=0){ stop(paste0("Blast command returned exit code ",exitCode," meaning something went wrong. Could be memory-killed, in which case, try breaking up the query (or subject) into fewer or smaller sequences.")) }
      tmp <- fread(oFile,col.names=outputColNames,colClasses=outputColClasses)
      unlink(oFile)
      if(!makeQfile){ tmp[,qFastaFname:=qFn][] }
      tmp
    })
    if(!makeSfile){ tmp[,sFastaFname:=sFn][] }
    tmp
  })

  if(argGiven(outFastaFname)){
    write.table(bl,outFastaFname,row.names=F,sep="\t",quote=F)
  }

  if(makeSfile){unlink(Sys.glob(paste0(subjectFname,".*")))}
  if(makeQfile){unlink(queryFname)}

  bl
}





