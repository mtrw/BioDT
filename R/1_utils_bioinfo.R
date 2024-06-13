#' Reverse Complement
#'
#' @param x A character vector. With DNA sequences if the results are to make sense. [no default]
#' @returns A character vector.
#' @description
#' Calculate the reverse complements of DNA sequences.
#' @details
#' Accepts (and preserves) uppercase and lowercase nucleotide codes, including IUPAC codes. Other characters will be left as-is.
#' @examples
#' rc("acgt---ACGT---RYSWKMBDHVN---ryswkmbdhvn")
#' rc(rc(c("AaA-CcC","GgG-TtT")))
#' @export
rc <- function(x){
  require(stringi)
  stringi::stri_reverse(x) %>% stringi::stri_trans_char("acgtACGTRYSWKMBDHVNryswkmbdhvn","tgcaTGCAYRSWMKVHDBNyrswmkvhdbn")
}


#' Delete '-' characters from strings
#'
#' @param x A character vector.  [no default]
#' @returns A character vector the same length of x, with all '-' characters omitted.
#' @details
#' A common use is to delete gaps from aligned sequences to they can become raw sequences again, and thus aligned to yet other sequences etc.
#' @examples
#' killGaps("AGCT--No-More-Gaps--TCGA")
#' @export
killGaps <- function(x){
  gsub("-","",x)
}


#' Print sequences as an alignment
#'
#' @param x A character vector. With equal-length aligned DNA sequences if the results are to make sense. [no default]
#' @returns NULL
#' @description
#' A BLAST-esque text alignment with bars ('|') showing matches between neighbouring sequences.
#' @details
#' Wraps to the console width. Coordinates are given relative to the sequence starts.
#' @examples
#' printAln(c(
#'   "AAAGGCTGAGTCGAGTC",
#'   "AAAxxCT--GxCGxG-C"
#' ))
#' @export
printAln <- function(s){
  w <- options()$width
  l <- nchar(s[1])
  for(i in seq(1,l,by=w)){
    cat("[",i,"] \n")
    for(i_Aln in 1:(length(s)-1)){
      cat(substr(s[i_Aln],i,min(i+w,l)),"\n")
      for(j in i:min(i+w,l)){
        if((substr(s[i_Aln],j,j)==substr(s[i_Aln+1],j,j)) & substr(s[i_Aln],j,j)!="-"){
          cat("|")
        } else {
          cat(" ")
        }
      }
      cat("\n")
    }
    cat(substr(s[length(s)],i,min(i+w,l)),"\n")
    cat("\n\n")
  }
  return(NULL) %>% invisible
}


#' Effective heterozygosity of a vector of alleles
#'
#' @param x A  vector. Canonically, listing alleles in a population at some site. [no default]
#' @param inclNAs A logical value. Delete NA entries before making He calculation? [FALSE]
#' @returns A numeric value
#' @description
#' Calculates the probability that two randomly drawn alleles will be different.
#' @details
#' Does not account for ploidy in any way. Is agnostic to what the given alleles are. c("A","G","C","T","-") will give the c("same","result","as","this","would"). No characters including gaps ("-") or NAs are recognised as special (by default).
#' @examples
#' He(c("A","A","G"))
#' @export
He <- function(x, inclNAs = F){
  if (inclNAs==F) {
    x <- x[!is.na(x)]
  }
  t <- table(x,useNA = "ifany")
  1-sum((t/sum(t))**2)
}

#' #' @export
# addStrand < function(alignmentDT){ # implement me
  # Should go on coordDT, alignDT
# }


#' Convert a vector of equal-length character strings into a matrix.
#'
#' @param x A (most usefully, character) vector.
#' @returns A character matrix
#' @description
#' One string per row, one char per entry.
#' @examples
#' charVec2Matrix(c("Hello","dolly","sings","Satch"))
#' @export
charVec2Matrix <- function(x){
  matrix(strsplit(paste0(x,collapse=""),"")[[1]],ncol=nchar(x[1]),nrow=length(x),byrow=T)
}


#' Convert an alignment from a vector of character strings to a matrix
#'
#' @param seqDT A data.table, meeting the specs for a seqDT.  [no default]
#' @returns A numeric value
#' @description
#' Useful for column-wise calculations (like consensus sequences, minor allele frequencies, etc.)
#' @details
#' One row per sequence, one entry per alignment column (i.e. per nucleotide/amino, or gap). Row names are taken from the column 'seqId'. The method is sequence-agnostic; any characters trings of equal length will work.
#' @examples
#' He(c("A","A","G"))
#' @export
alnSeq2seqMatrix <- function(seqDT){
  is_seqDT(seqDT,croak = T)
  if(!all(nchar(seqDT$seq)==nchar(seqDT$seq[1]))){ stop("Sequences in 'seq' column must all be the same length.") }
  m <- charVec2Matrix(seqDT$seq)
  colnames(m) <- seqDT$seqId
  m
}



#' Plots alignments with coloured nucleotides
#'
#' @param seqDT A data.table, meeting the specs for a seqDT. [no default]
#' @param rowGap The size of the gap between alignment rows [0.05]
#' @returns A numeric value
#' @description
#' Alignments are plotted horizontally, with the first given sequence on top.
#' @details
#' To plot from an alignmentDT, do a manual conversion. E.g., to 'melt' the alignments:
#' ```alignmentDT[,alnId:=paste0("Alignment_",1:.N)][,seqId=c(sSeqId,qSeqId),seq=c(sAlnSeq,qAlnSeq),by=.(alnId)]```
#' or call `alignmentDT2seqDT()` (be aware this function preserves only one of the subject or query sequences). Also, be aware most alignmentDTs have different length alignments, so you will typically be selecting a limited number of rows.
#' @examples
#' seqDT <- data.table(
#'   seq = c("GC-AT",
#'           "TC-AA"),
#'   seq2 = c("ABCD","EFGH"),
#'   seqId = LETTERS[1:2]
#' )
#' plotNucleotideAlignment(seqDT)
#' plotNucleotideAlignment(seqDT,invariantColour = "black")
#' plotNucleotideAlignment(seqDT,alnItemCols = applyPalette(1:4,palettePresets$williams),gapColour = "white")
#' plotNucleotideAlignment(seqDT[,.(seqId,seq=seq2)])
#' @export
plotNucleotideAlignment <- function(seqDT,rowGap=0.05,alnItemCols=c("#990000","#000099","#009900","#BB9911"),invariantColour=NULL,gapColour="#BBBBBB",unknownNtColour="#FF55FF",labelsInAxis=TRUE,cex.axis=0.7,...){
  is_seqDT(seqDT,croak = T)
  sm <- alnSeq2seqMatrix(seqDT)
  plotTop <- (1+rowGap)*nrow(sm)
  seqNames <- seqDT$seqId

  if(any(!sm %in% c("A","a","T","t","G","g","C","c","-"))){
    warning("Unrecognised nucleotides (i.e. [^AGCTacgt-]) in alignments! These will be coloured according to `unknownNtColour`.")
  }

  colsm <- apply(sm,2,function(c){
    if(!is.null(invariantColour) & all(c==c[1]) & c[1]!="-"){
      rep(invariantColour,length(c))
    } else {
      return(swap(c,c("A","a","T","t","G","g","C","c","-"),c(alnItemCols[rep(1:4,each=2)],gapColour)))
    }
  })

  colsm[!sm %in% c("A","a","T","t","G","g","C","c","-")] <- unknownNtColour

  #par(mar=c(.5,15,.5,.5))
  null_plot(c(0,ncol(colsm)),c(0,plotTop),xaxt="n",yaxt="n",cex.axis=.3)
  for(i in 1:nrow(colsm)){
    for(j in 1:ncol(colsm)){
      rect(
        xleft = j-1,
        xright = j,
        ytop = plotTop-((i-1)*(1+rowGap)),
        ybottom = plotTop-((i-1)*(1+rowGap))-1,
        border=NA,
        col = colsm[i,j]
      )
    }
  }
  yPts <- plotTop-(((1:nrow(colsm))-1)*(1+rowGap))-0.5
  if(labelsInAxis==TRUE){ axis(2,at=yPts,labels = seqNames,las=2,cex.axis=cex.axis,...) }
  yPts
}




#' @export
pctIdMat <- function(seqDT=NULL,seqMatrix=NULL,inclGaps=T){
  if(all(is.null(c(seqDT,seqMatrix)))){ stop("At least one of `seqId` and `seqMatrix` must be given") }
  if(!is.null(seqDT)){
    is_seqDT(seqDT,croak = T)
    seqMatrix <- alnSeq2seqMatrix(seqDT)
  }

  out <- matrix(NA,nrow=nrow(seqMatrix),ncol=nrow(seqMatrix))
  colnames(out) <- rownames(out) <- seqDT$seqId


  if(inclGaps==TRUE){
    for(i in 1:(nrow(seqMatrix)-1)){
      for(ii in (i+1):nrow(seqMatrix)){
        #dev i<-1; ii<-2
        selCol <- which(!(seqMatrix[i,]=="-" & seqMatrix[ii,]=="-"))
        out[ii,i] <- sum(seqMatrix[i,selCol]==seqMatrix[ii,selCol])/length(selCol)
      }
    }
  } else {
    for(i in 1:(nrow(seqMatrix)-1)){
      for(ii in (i+1):nrow(seqMatrix)){
        selCol <- which(seqMatrix[i,]!="-" | seqMatrix[ii,]!="-")
        out[ii,i] <- sum(seqMatrix[i,selCol]==seqMatrix[ii,selCol])/length(selCol)
      }
    }
  }
  d <- as.dist(out)
}




#' @export
clumpCoordDT <- function(coordDT,distCutoff=0){
  is_coordDT(coordDT,croak=T)
  c <- copy(coordDT)
  c[end<start,c("start","end"):=.(end,start)]
  setkey(c,start,end)
  c[,pmxesf:=c(0,maxSoFar(end)[1:(.N-1)]),by=.(seqId)] # previous max end so far
  #c[,mnssf:=minSoFar(start),by=.(seqId)]
  c[,newBlockStart:=((start-distCutoff)>pmxesf)+0L,by=.(seqId)]
  c[,block:=cumsum(newBlockStart),by=.(seqId)]
  c[]
}

#' @export
unionCoordDT <- function(coordDT,distCutoff=0){
  c <- clumpCoordDT(coordDT,distCutoff)
  c[,.(start=min(start),end=(max(end))),by=.(seqId,block)]
}



