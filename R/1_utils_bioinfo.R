
# Reverse complement a text DNA sequence.
#rc(c("AAAAAC","AgCtaGcTyYyrSdNbxxxx----agcT"))
#' @export
rc <- function(x){
  stringi::stri_reverse(x) %>% stringi::stri_trans_char("acgtACGTRYSWKMBDHVNryswkmbdhvn","tgcaTGCAYRSWMKVHDBNyrswmkvhdbn")
}


# Print seqs as if an alignment
# Give it strings in (e.g.) AAG-AT / A-GCAT format. s1[1] will align with s2[1], etc
# Coords are simple alignment positions
# If s2 is NULL, will "multiple align" all s1s
#' @export
printAln <- function(s1,s2=NULL){
  w <- options()$width
  l <- stri_length(s1[1])
  for(i in seq(1,l,by=w)){
    cat("[",i,"] \n")
    for(i_Aln in 1:(length(s1)-1)){
      cat(substr(s1[i_Aln],i,min(i+w,l)),"\n")
      for(j in i:min(i+w,l)){
        if((substr(s1[i_Aln],j,j)==substr(s1[i_Aln+1],j,j)) & substr(s1[i_Aln],j,j)!="-"){
          cat("|")
        } else {
          cat(" ")
        }
      }
      cat("\n")
    }
    cat(substr(s1[length(s1)],i,min(i+w,l)),"\n")
    cat("\n\n")
  }
}
#printAln(c("aagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgta","aagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgta","aagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgta","aagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgtaaagtcgta"))



# CHANGE NAME TO SOMETHING MORE SENSIBLE SOON
#' @export
alnGetClumps <- function(alnDT,distCutoff=1e3L,hclustAlgorithm="single"){ # CHANGE TO COORDS GET CLUMPS AND BUILD IN aln2coords CONVERSION!!!
  #alnDT <- alnRph12toPgF; distCutoff <- 1e4; nameColx = "hitId"; hclustAlgorithm="single"
  alnDT <- copy(alnDT)
  clumps <- alnDT[,{
    #debugonce(gridApplyDT)
    dm <- gridApplyDT(dtx=.SD,FUN=function(x,y){ min( abs(y$sStart-x$sEnd) , abs(x$sStart-y$sEnd) ) },symmetrical=TRUE) %>% as.dist()
    .(
      clump = hclust(dm,method=hclustAlgorithm) %>% cutree(h=distCutoff),
      sStart,
      sEnd
    )
  },by=.(sSeqId)]
  clumps[,.(span=diff(range(c(sStart,sEnd))),nAlns=.N,start=min(c(sStart,sEnd)),end=max(c(sStart,sEnd))),by=.(sSeqId,clump)]
}

#' @export
He <- function(x, inclNAs = F){
  if (inclNAs==F) {
    x <- x[!is.na(x)]
  }
  t <- table(x,useNA = "ifany")
  1-sum((t/sum(t))**2)
}

#' #' @export
# addStrand < function(alignmentDT){
# Should go on coordDT, alignDT
# }
