#' @export
has_SeqId <- makeValidator( "has_seqId" , c("seqId") )
#' @export
has_strand <-makeValidator( "has_strand" , c("strand") )
#' @export
has_name <- makeValidator( "has_name" , c("name") )
#' @export
has_score <- makeValidator( "has_score" , c("score") )
#' @export
has_start_end <- makeValidator( "has_start_end" , c("start","end") )
#' @export
has_start <- makeValidator( "has_start" , c("start") )
#' @export
has_end <- makeValidator( "has_end" , c("end") )
#' @export
has_sqAlnSeqs <- makeValidator( "has_sqAlnSeqs" , c("sAlnSeq","qAlnSeq") )
#' @export
has_sAlnSeq <-makeValidator( "has_sAlnSeq" , c("sAlnSeq") )
#' @export
has_qAlnSeq <-makeValidator( "has_qAlnSeq" , c("qAlnSeq") )
#' @export
has_fastaFname <- makeValidator( "fastaFname" , c("fastaFname") )
#' @export
has_sFastaFname <- makeValidator( "has_sFastaFname" , c("sFastaFname") )
#' @export
has_qFastaFname <- makeValidator( "has_qFastaFname" , c("qFastaFname") )
#' @export
is_seqDT <- makeValidator( "is_seqDT" , c("seqId","seq") )
#' @export
is_seqDT_with_seqFname <- makeValidator( "is_seqDT_with_seqFname" , c("seqId","seq","seqFname") )
#' @export
is_coordDT <- makeValidator( "is_coordDT" , c("seqId","start","end") )
#' @export
is_seqListDT <- makeValidator( "is_seqListDT" , c("seqId") )
#' @export
has_minSize_maxSize <- makeValidator( "has_minSize_maxSize" , c("minSize","maxSize") )
#' @export
has_optimalSize <- makeValidator( "has_optimalSize" , c("optimalSize") )
#' @export
has_internal <- makeValidator( "has_internal" , c("internal") )
#' @export
has_minProductSize <- makeValidator( "has_minProductSize" , c("productMinSize") )
#' @export
has_maxProductSize <- makeValidator( "has_maxProductSize" , c("productMaxSize") )
#' @export
is_alignmentDT <- makeValidator( "is_alignmentDT" , c("sSeqId","qSeqId","sStart","sEnd","qStart","qEnd") )

