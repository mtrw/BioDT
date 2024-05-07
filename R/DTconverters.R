#' @export
#'
alignmentDT2coordDT <- function(alignmentDT,coordsOf="subject",numberSeqIds=FALSE,keepCols=c("noAlnSeqs","minimal","all"),keepColList=NULL){
  is_alignmentDT(alignmentDT,croak=TRUE)
  out <- if(coordsOf=="subject"){
    alignmentDT[,
      seqId:=if(numberSeqIds==TRUE){paste0(sSeqId,"_",1:.N)}else{sSeqId}][,
      start:=pmin(sStart,sEnd)][,
      end  :=pmax(sStart,sEnd)
    ]
    if(has_sFastaFname(alignmentDT)){ alignmentDT[,fastaFname:=sFastaFname] }
  } else if (coordsOf=="query"){
    alignmentDT[,
      seqId := if(numberSeqIds==TRUE){ paste0(qSeqId,"_",1:.N) }else{ qSeqId }][,
      start := qStart ][,
      end   := qEnd
    ]
    if(has_qFastaFname(alignmentDT)){ alignmentDT[,fastaFname:=qFastaFname] }
  }
  selCols <- if( argGiven(keepColList) ){
    if(any(!keepColList %in% colnames(alignmentDT))){ stop(paste0("Requested columns are not in the provided alignmentDT: ",setdiff(keepColList,colnames(alignmentDT)))) }
    union(keepColList,c("seqId","start","end"))
  } else if (keepCols[1]=="minimal"){
    c("seqId","start","end",union(colnames(alignmentDT),c("fastaFname")))
  } else if (keepCols[1]=="noAlnSeqs"){
    setdiff(colnames(alignmentDT),c("sAlnSeq","qAlnSeq") )
  } else if (keepCols[1]=="all"){
    colnames(alignmentDT)
  }
  is_coordDT(out,croak=T,message = "After conversion, the result was not a legitimate coordDT. This probably shouldn't happen. Please let Tim know")
  out[,..selCols]
}






#' @export
alignmentDT2seqDT <- function(alignmentDT,seqOf="subject",numberSeqIds=FALSE,keepCols=c("noAlnSeqs","minimal","all"),keepColList=NULL){
  is_alignmentDT(alignmentDT,croak=TRUE)
  has_sqAlnSeqs(alignmentDT,croak=TRUE)
  out <- if(seqOf=="subject"){
    alignmentDT[,
      seqId:=if(numberSeqIds==TRUE){paste0(sSeqId,"_",1:.N)}else{sSeqId}][,
      seq:=fifelse(sEnd<sStart,gsub("-","",sAlnSeq),gsub("-","",rc(sAlnSeq)))
    ]
  } else if (seqOf=="query"){
    alignmentDT[,
      seqId:=if(numberSeqIds==TRUE){paste0(sSeqId,"_",1:.N)}else{sSeqId}][,
      seq:=gsub("-","",qAlnSeq)
    ]
  }
  selCols <- if( argGiven(keepColList) ){
    if(any(!keepColList %in% colnames(alignmentDT))){ stop(paste0("Requested columns are not in the provided alignmentDT: ",setdiff(keepColList,colnames(alignmentDT)))) }
    union(keepColList,c("seqId","seq"))
  } else if (keepCols[1]=="minimal"){
    c("seqId","seq")
  } else if (keepCols[1]=="noAlnSeqs") {
    setdiff(colnames(alignmentDT),c("sAlnSeq","qAlnSeq"))
  } else if (keepCols[1]=="all") {
    colnames(alignmentDT)
  }
  is_seqDT(out,croak=T,message = "After conversion, the result was not a legitimate coordDT. This probably shouldn't happen. Please let Tim know")
  out[,..selCols]
}


