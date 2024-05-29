#' @export
clusterSeqs <- function(
    seqDT=NULL,
    fastaFname=NULL,
    program=c("cd-hit","uclust"),
    seqIdThres=0.8,
    cdhit_wordLength=5,
    cdhit_global=FALSE,
    cdhit_minCovPropLongerSeq=0.2,
    cdhit_minCovPropShorterSeq=0.8,
    cdhit_minAlnLength=10,
    cdhit_fast=FALSE,
    moreCommandArgs="",
    programBinary=system(paste0("which ",program),intern=T),
    awkBinary=system("which awk",intern=T),
    returnOnlyClusters=FALSE
){
  # DEV #################
  # seqDT=NULL
  # fastaFname="/data/gpfs/projects/punim1869/shared_data/misc_sequence/btr_btrLike_queries_Morex_GP.fasta"
  # program="cd-hit"
  # cdhit_seqIdThres=0.8
  # cdhit_global=FALSE
  # cdhit_minCovPropLongerSeq=0.2
  # cdhit_minCovPropShorterSeq=0.8
  # cdhit_minAlnLength=10
  # cdhit_fast=FALSE
  # moreCommandArgs=""
  # programBinary=system(paste0("which ",program),intern=T)
  # awkBinary=system(paste0("which awk"),intern=T)
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
  writeFasta(seqDT_[,.(seqId=seqId_noSpace,seq)],outFastaFname=tfFastaFname)

  if(program[1]=="cd-hit"){
    tfOutFastaFname <- tempfile(fileext = ".out.fasta")
    tfOutClstrFname <- paste0(tfOutFastaFname,".clstr")

    commandCdhit <- paste0(
      programBinary,
      " -B 1 -d 0 -i ",tfFastaFname,
      " -o ",  tfOutFastaFname,
      " -c ",  seqIdThres,
      " -n ",  cdhit_wordLength,
      " -G ",  as.integer(cdhit_global),
      " -aS ", cdhit_minCovPropShorterSeq,
      " -aL ", cdhit_minCovPropLongerSeq,
      " -A ",  cdhit_minAlnLength,
      " -g ",  as.integer(cdhit_fast)
    )
    commandParse <- paste0("cat ",tfOutClstrFname," | ",awkBinary," '",aScriptCdhitClstr2tbl,"'")

    ce("Running command: ",commandCdhit)
    system(commandCdhit)
    #system(paste0("cat ",tfOutFastaFname))
    #system(paste0("cat ",tfOutClstrFname))
    ce("Running command: ",commandParse)
    clAssign <- fread(cmd=commandParse,header=F,col.names = c("seqId_noSpace","cluster"))
    unlink(c(tfFastaFname,tfOutFastaFname,tfOutClstrFname))
    if(returnOnlyClusters==TRUE){return(clAssign)}
    return(clAssign[seqDT_,on=.(seqId_noSpace)][,seqId_noSpace:=NULL][])
  } else if (program=="uclust"){
    commandUclust <- paste0(
        programBinary,
        " -cluster_fast ",tfFastaFname,
        " -id ",seqIdThres,
        " -uc /dev/stdout | awk '",aScriptUclust2Tbl,"'"
    )
    ce("Running command: ",commandUclust)
    clAssign <- fread(cmd=commandUclust,header=F,col.names = c("seqId_noSpace","cluster"))
    unlink(tfFastaFname)
    if(returnOnlyClusters==TRUE){return(clAssign)}
    return(clAssign[seqDT_,on=.(seqId_noSpace)][,seqId_noSpace:=NULL][])
  } else {
    stop(paste0("Program ",program," not recognised or implemented, sorry."))
  }
}





