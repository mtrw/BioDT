
#' @export
msa <- function(seqDT,method="ClustalOmega",seqType=c("guess","dna","aa"),...){
  require(msa)
  require(Biostrings)
  is_seqDT(seqDT,objName = deparse(substitute(seqDT)))
  stringSetFn <- if(seqType[1]=="dna"){
    Biostrings::DNAStringSet
  } else if (seqType[1]=="aa"){
    Biostrings::AAStringSet
  } else if(seqType[1]=="guess"){
    testMe <- substr(seqDT$seq[1],1,1000)
    testRegex <- paste0("[^",paste(ntCodesLegal,collapse=""),"]")
    if(!grepl(testRegex,testMe)){
      ce(
        "First 100bp of first sequence contained only entries in the set ",
        paste0("[",paste(ntCodesLegal,collapse=""),"]."),
        " Guessing these are nucleotide sequences. To set it manually use 'seqType=\"dna\" or seqType=\"aa\".'"
      )
      Biostrings::DNAStringSet
    } else {
      ce(
        "First 100bp of first sequence contained entries not in the set ",
        paste0("[",paste(ntCodesLegal,collapse=""),"]."),
        " Guessing these are amino acid sequences. To set it manually use 'seqType=\"dna\" or seqType=\"aa\".'"
      )
      Biostrings::AAStringSet
    }
  }

  seqDT[,alnSeq:={
    seq %>%
      stringSetFn() %>%
      msa::msa(method=method,order="input",...) %>%
      msaConvert(type = "seqinr::alignment") %>%
      `[[`("seq")
  }][]
}


#' @export
localAlign <- function(seqDT1,seqDT2,seqType=c("guess","dna","aa"),subMat=pwalign::nucleotideSubstitutionMatrix(),...){
  require(pwalign)
  is_seqDT(seqDT1,objName = deparse(substitute(seqDT)))
  is_seqDT(seqDT2,objName = deparse(substitute(seqDT)))
  if(nrow(seqDT1)!=1){ stop("seqDT1 must have one row providing a sequence to be aligned") }
  if(nrow(seqDT2)!=1){ stop("seqDT2 must have one row providing a sequence to be aligned") }
  stringSetFn <- if(seqType[1]=="dna"){
    Biostrings::DNAString
  } else if (seqType[1]=="aa"){
    Biostrings::AAString
  } else if(seqType[1]=="guess"){
    testMe <- substr(seqDT1$seq[1],1,1000)
    testRegex <- paste0("[^",paste(ntCodesLegal,collapse=""),"]")
    if(!grepl(testRegex,testMe)){
      ce(
        "First 100bp of first sequence of seqDT1 contained only entries in the set ",
        paste0("[",paste(ntCodesLegal,collapse=""),"]."),
        " Guessing these are nucleotide sequences. To set it manually use 'seqType=\"dna\" or seqType=\"aa\".'"
      )
      Biostrings::DNAStringSet
    } else {
      ce(
        "First 100bp of first sequence of seqDT1 contained entries not in the set ",
        paste0("[",paste(ntCodesLegal,collapse=""),"]."),
        " Guessing these are amino acid sequences. To set it manually use 'seqType=\"dna\" or seqType=\"aa\".'"
      )
      Biostrings::AAStringSet
    }
  }

  a <- pairwiseAlignment(seqDT1$seq,seqDT2$seq,type="local",substitutionMatrix=subMat)

  data.table(
    alnSeq=c(
      alignedPattern(a) %>% as.character,
      alignedSubject(a) %>% as.character
    ),
    score=score(a)
  )
}

