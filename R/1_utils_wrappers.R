
#' @export
callLastz <- function(
    binPath,
    sFile,
    sSubset,
    sRange,qFile,
    qSubset,
    qRange,
    outFmtArgs,
    outFmtColNames=c("sSeqId","sStart","sEnd","sLength","sStrand","qSeqId","qStart","qEnd","qLength","qStrand","pctId","pctId_noGaps"),
    otherArgs="",
    showCmd=FALSE
){
  # binPath=system("which lastz",intern=T)
  # sFile="/data/gpfs/projects/punim1869/shared_data/misc_sequence/btr_btrLike_queries_Morex_GP.fasta"
  # sSubset="Btr1_like_b_5_Morex_Btr1_like_b1"
  # sRange=c(1,450)
  # qFile="/data/gpfs/projects/punim1869/shared_data/misc_sequence/btr_btrLike_queries_Morex_GP.fasta"
  # qSubset="Btr1_functional_Morex"
  # qRange=c(5,500)
  # outFmtArgs = "general:name1,start1,end1,length1,strand1,name2,start2,end2,length2,strand2,id%,blastid%"
  # outFmtColNames=c("sSeqId","sStart","sEnd","sLength","sStrand","qSeqId","qStart","qEnd","qLength","qStrand","pctId","pctId_noGaps")
  # otherArgs=""

  sSubsetArg <- if(!is.null(sSubset)){paste0("\\[subset=<(echo \"",sSubset,"\")\\]")}else{""}
  qSubsetArg <- if(!is.null(qSubset)){paste0("\\[subset=<(echo \"",qSubset,"\")\\]")}else{""}
  sRangeArg  <- if(!is.null(sRange)){paste0("\\[",sRange[1],"..",sRange[2],"\\]")}else{""}
  qRangeArg  <- if(!is.null(qRange)){paste0("\\[",qRange[1],"..",qRange[2],"\\]")}else{""}

  cmd <- paste0(binPath," ",sFile,sSubsetArg,sRangeArg," ",qFile,qSubsetArg,qRangeArg," --format=",outFmtArgs," ",otherArgs)#,"\["sSubset,sRange,qFile,qSubset,qRange,outFmtArgs,otherArgs)
  if(showCmd==TRUE){ ce(cmd) }
  fread(cmd=cmd,col.names=outFmtColNames)
}
