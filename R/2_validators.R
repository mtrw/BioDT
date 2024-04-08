#' @export
has_seqFname <- makeValidator( "has_seqFname" , c("seqFname") )
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
is_seqDT <- makeValidator( "is_seqDT" , c("seqId","seq") )
#' @export
is_seqDT_with_seqFname <- makeValidator( "is_seqDT_with_seqFname" , c("seqId","seq","seqFname") )
#' @export
is_coordDT <- makeValidator( "is_coordDT" , c("seqId","start","end") )
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




