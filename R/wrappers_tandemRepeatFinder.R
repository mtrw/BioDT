#' @export
tandemRepeatFinder <- function(
    seqDT=NULL,
    fastaFname=NULL,
    trf_matchWeight=2L,
    trf_mismatchPenalty=7L,
    trf_indelPenalty=7L,
    trf_matchProb=c(80L,75L),
    trf_indelProb=c(10L,20L),
    trf_minScore=20L,
    trf_maxPeriod=2000L,
    moreCommandArgs="",
    trfBinary=system("which trf",intern=T),
    awkBinary=system("which awk",intern=T)
){
  # DEV #################
  # seqDT=NULL
  # fastaFname=NULL
  # trf_matchWeight=2L
  # trf_mismatchPenalty=7L
  # trf_indelPenalty=7L
  # trf_matchProb=c(80L,75L)
  # trf_indelProb=c(10L,20L)
  # trf_minScore=20L
  # trf_maxPeriod=2000L
  # moreCommandArgs=""
  # trfBinary=system("which trf",intern=T)
  # awkBinary=system("which awk",intern=T)
  # \DEV ###############

  seqDT_ <- copy(seqDT)

  if((argGiven(fastaFname) & argGiven(seqDT)) | (argNotGiven(fastaFname) & argNotGiven(seqDT))){
    stop("Sequence must be provided via one (and only one) of either `fastaFname` (a .fasta or fasta.gz file) or `seqDT` (a data.table meeting seqDT specs).")
  }

  if(any(duplicated(seqDT_$seqId))){
    warning("For clustering, the `seqId` column must contain unique values. seqIds will therefore be numbered, and a new column `seqId_original` will be added to contain the original seqIds.")
    seqDT_[,seqId_original:=seqId]
    seqDT_[,seqId:=paste0(seqId,"_",1:.N)]
  }

  if(argGiven(fastaFname)){
    seqDT_ <- readFasta(fastaFname)
  }

  seqDT_[,seqId_noSpace:=gsub("[[:space:]]","",seqId)]

  tfFastaFname <- tempfile(fileext = ".fasta")
  tfOutFastaFname <- tempfile(fileext = ".out.fasta")
  tfOutClstrFname <- paste0(tfOutFastaFname,".clstr")

  writeFasta(seqDT_[,.(seqId=seqId_noSpace,seq)],outFastaFname=tfFastaFname)

  command <- paste0(
    trfBinary,
    " ",tfFastaFname,
    " ", trf_matchWeight,
    " ", trf_mismatchPenalty,
    " ", trf_indelPenalty,
    " ", trf_matchProb[1],
    " ", trf_indelProb[1],
    " ", trf_minScore,
    " ", trf_maxPeriod,
    " -ngs ",
    " | ",awkBinary," '",aScriptTrf2Tbl,"'"
  )

  rpts <- fread(cmd=command,header=F,select=c(1:5,16),col.names = c("seqId_noSpace","start","end","nCopies","meanLength","repeatConsensusSeq"),colClasses = list(character=c(1,16),numeric=c(2,3,5),integer=4))[,repeatArrayLength:=abs((end-start)+1)]

  unlink(tfFastaFname)

  return(rpts[seqDT_,on=.(seqId_noSpace)][,seqId_noSpace:=NULL][])
  rpts[seqDT_,on=.(seqId_noSpace)][,seqId_noSpace:=NULL][]
}





