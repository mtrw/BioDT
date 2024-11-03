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
alignedSeqDT2seqMatrix <- function(seqDT){
  is_alignedSeqDT(seqDT)
  if(!all(nchar(seqDT$alnSeq)==nchar(seqDT$alnSeq[1]))){ stop("Sequences in 'seq' column must all be the same length.") }
  m <- charVec2Matrix(seqDT$alnSeq)
  rownames(m) <- seqDT$seqId
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
#' ```alignmentDT[,alnId:=paste0("Alignment_",1:.N)][,seqId=c(sSeqId,qSeqId),alnSeq=c(sAlnSeq,qAlnSeq),by=.(alnId)]```
#' or call `alignedSeqDT2seqMatrix()` (be aware this function preserves only one of the subject or query sequences). Also, be aware most alignmentDTs have different length alignments, so you will typically be selecting a limited number of rows.
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
plotSequenceAlignment <- function(seqDT,rowGap=0.05,alnPreset=c("guess","dna","aa"),alnItems=NULL,alnItemCols=NULL,invariantColour=NULL,gapColour="#BBBBBB",unknownNtColour="#FF55FF",labelsInAxis=TRUE,cex.axis=0.7,...){
  is_alignedSeqDT(seqDT,croak = T)

  if( sum(is.null(alnItemCols),is.null(alnItems))==1 | length(alnItemCols) != length(alnItems) ){ stop("With custom colours, both 'alnItems' and 'alnItemCols' must be given, and these must be the same length ...") }

  if(argNotGiven(alnItems)){

    if (alnPreset[1]=="guess"){
      testMe <- substr(seqDT$seq[1],1,1000)
      testRegex <- paste0("[^",paste(ntCodesLegal,collapse=""),"]")
      if(!grepl(testRegex,testMe)){
        ce(
          "First 100bp of first sequence contained only entries in the set ",
          paste0("[",paste(ntCodesLegal,collapse=""),"]."),
          " Guessing these are nucleotide sequences. To set it manually use 'alnPreset=\"dna\" or alnPreset=\"aa\".'"
        )
        alnPreset[1] <- "dna"
      } else {
        ce(
          "First 100bp of first sequence contained entries not in the set ",
          paste0("[",paste(ntCodesLegal,collapse=""),"]."),
          " Guessing these are amino acid sequences. To set it manually use 'alnPreset=\"dna\" or alnPreset=\"aa\".'"
        )
        alnPreset[1] <- "aa"
      }
    }

    if(alnPreset[1]=="dna"){
      alnItems <- c("A","C","G","T")
      alnItemCols <- c("#990000","#000099","#009900","#BB9911")
    } else if (alnPreset[1]=="aa"){
      alnItems <- aaCodesLegalUcOnlyNogap
      alnItemCols <- c(BioDT::applyPalette(1:(length(alnItems)-1),colChain=palettePresets$wheel),"#000000")
    }
  }

  sm <- alignedSeqDT2seqMatrix(seqDT)
  plotTop <- (1+rowGap)*nrow(sm)
  seqNames <- seqDT$seqId

  if(any(!sm %in% c(alnItems,alnItems %>% tolower,"-"))){
    warning(paste0("Some of the items (nucleotide or amino codes) in the alignments are not in 'alnItems'. These will be coloured according to `unknownNtColour`: ",sm[!sm %in% c(alnItems,alnItems %>% tolower,"-")] %>% unique %>% paste(collapse="")))
  }

  colsm <- apply(sm,2,function(c){
    if(!is.null(invariantColour) & all(c==c[1]) & c[1]!="-"){
      rep(invariantColour,length(c))
    } else {
      return( swap( c, c(alnItems %>% toupper,alnItems %>% tolower,"-"), c(rep(alnItemCols,2),gapColour)) )
    }
  })

  colsm[!sm %in% c(alnItems %>% toupper,alnItems %>% tolower,"-")] <- unknownNtColour

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


#' Cartesian alignment plots
#'
#' @param alignmentDT A data.table, meeting the specs for an alignmentDT. [no default]
#' @param newPlot Add to an existing plot, or create a new one? [FALSE]
#' @param bgColour ["#000011", aka a very dark navy blue]
#' @param colourAlignments If the alignments are to be coloured, then name the column from which to draw them. See examples. [NULL]
#' @param alignmentColourDefault If 'colourAlignments' isn't used, then what colour should they be? ["#FFFFEE", aka light cream]
#' @param lwd The line width [1.0]
#' @param xlim Use to limit the range along x (the 'subject sequence' axis) [NULL]
#' @param ylim Use to limit the range along x (the 'subject sequence' axis) [NULL]
#' @returns Nothing, it draws a plot.
#' @description
#' Often called a dotplot for historical (read "invalid") reasons. That's a stupid name since the word dot just doesn't relate to anything going on here in any way. Hence, the shiny new name "cartesian alignment plot" (c.f., say, Sankey alignment plots).
#' @details
#' This assumes you've included just one subject and query seq (i.e. `nu(c(sSeqId,qSeqId))==2`). In future this check may be made explicit. If you do include more, the alignments will all be overlapped on one plot, which is not likely to be what you intended.
#' @examples
#' Coming ...
#' @export
plotAlignmentCartesian <- function(alignmentDT,newPlot=TRUE,bgColour="#000011",colourAlignments=NULL,alignmentColourDefault="#FFFFEE",lwd=1.0,xlim=NULL,ylim=NULL,...){
  is_alignmentDT(alignmentDT,croak=T)
  if(alignmentDT[,nu(sSeqId)>1 | nu(qSeqId)>1]){
    ce("WARNING: alignmentDT has more than one unique 'sSeqId's and/or 'qSeqId's. All alignments will be plotted on the same axes; This probably not what you want.")
    wait("Please press <enter> to continue and bear responsiblity for the consequences, or <esc> to stop and try again.")
  }
  a <- copy(alignmentDT)
  a[,idx:=1:.N]
  a[,col:=if(argGiven(colourAlignments)){get(colourAlignments)}else{alignmentColourDefault}]

  is_alignmentDT(a,croak=T)
  xLims <- if(is.null(xlim)) {a[,range(sStart,sEnd,qStart,qEnd)]} else {xlim}
  yLims <- if(is.null(ylim)) {a[,range(sStart,sEnd,qStart,qEnd)]} else {ylim}

  if(newPlot==TRUE){
    null_plot(xLims,yLims,...)
    #null_plot(linesLims,linesLims)
  }
  rect(xLims[1],xLims[1],yLims[2],yLims[2],border = NA,col = bgColour)

  a[.N:1,{
    # wait()
    # print(.SD)
    lines(
      x=c(sStart,sEnd),
      y=c(qStart,qEnd),
      col=c(col,col),
      lwd=lwd
    )
  },by=.(idx)]

  return(NULL) %>% invisible
}




#' @export
pctIdMat <- function(seqDT=NULL,seqMatrix=NULL,inclGaps=T){
  if(all(is.null(c(seqDT,seqMatrix)))){ stop("At least one of `seqId` and `seqMatrix` must be given") }
  if(!is.null(seqDT)){
    is_seqDT(seqDT,croak = T)
    seqMatrix <- alignedSeqDT2seqMatrix(seqDT)
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
  if(any(c("clumpId","newBlockStart","pmxsf") %in% colnames(coordDT))){
    warning("Columns named 'clumpId', 'newBlockStart' and/or 'pmxsf' in input will be overwritten or deleted. If you wish to keep them, please rename them.")
  }
  c <- copy(coordDT)
  c[end<start,c("start","end"):=.(end,start)]
  setkey(c,start,end)
  c[,pmxesf:=if(.N==1){end}else{c(0,maxSoFar(end)[1:(.N-1)])},by=.(seqId)] # previous max end so far
  #c[,mnssf:=minSoFar(start),by=.(seqId)]
  c[,newBlockStart:=((start-distCutoff)>pmxesf)+0L,by=.(seqId)]
  c[,clumpId:=cumsum(newBlockStart),by=.(seqId)]

  c[,pmxesf:=NULL][,newBlockStart:=NULL][]
}

#' @export
unionCoordDT <- function(coordDT,distCutoff=0){
  is_coordDT(coordDT,croak = TRUE)
  c <- clumpCoordDT(coordDT,distCutoff)
  c[,.(start=min(start),end=(max(end))),by=.(seqId,clumpId)][,clumpId:=NULL][]
}


#' @export
noNestAlignmentDT <- function(alignmentDT,markOnly=F){
  #alignmentDT=a
  is_alignmentDT(alignmentDT,croak=TRUE)
  if(any(c("sStart_unstranded","sEnd_unstranded","qStart_unstranded","qEnd_unstranded","sSpanEnd","qSpanStart","qSpanEnd") %in% colnames(alignmentDT))){
    warning("Any columns in alignmentDT named \"sStart_unstranded\", \"sEnd_unstranded\", \"qStart_unstranded\", \"qEnd_unstranded\", \"sSpanEnd\", \"qSpanStart\", and  \"qSpanEnd\" will be overwritten or deleted in output. If you wish to preserve them, rename them.")
  }
  alignmentDT <- copy(alignmentDT)
  alignmentDT[,sStart_unstranded:=pmin(sStart,sEnd)]
  alignmentDT[,sEnd_unstranded:=pmax(sStart,sEnd)]
  alignmentDT[,qStart_unstranded:=pmin(qStart,qEnd)]
  alignmentDT[,qEnd_unstranded:=pmax(qStart,qEnd)]
  setorder(alignmentDT,sStart,-sEnd)
  alignmentDT[,sSpanEnd:=maxSoFar(sEnd_unstranded)]
  setkey(alignmentDT,sSpanEnd,sStart)
  alignmentDT[,qSpanStart:=qStart[1],by=.(sSpanEnd)]
  alignmentDT[,qSpanEnd:=qEnd[1],by=.(sSpanEnd)]
  alignmentDT[1:25]
  alignmentDT[,nested:=((1:.N)>1) & (qStart>=qSpanStart) & (qEnd<=qSpanEnd),by=.(sSpanEnd)]
  keepCols <- setdiff(colnames(alignmentDT),c("sSpanEnd","qSpanStart","qSpanEnd"))
  if(markOnly==T){
    keepCols <- setdiff(colnames(alignmentDT),c("sSpanEnd","qSpanStart","qSpanEnd"))
    return(alignmentDT[,..keepCols])
  } else {
    keepCols <- setdiff(colnames(alignmentDT),c("sSpanEnd","qSpanStart","qSpanEnd","nested"))
    return(alignmentDT[nested==FALSE][,nested:=NULL][,..keepCols])
  }
}


#' @export
translateCoordsClumpGapkill <- function(coordDT=NULL,coordVec=NULL,clumpCoordDT,interClumpGap=1e3){
  if(argGiven(coordDT) & argNotGiven(coordVec)){
    coordDT <- copy(coordDT)
  } else if (argGiven(coordVec) & argNotGiven(coordDT)){
    warning("'seqId' column of clumpCoordDT will not be used and all entries will be set to \"\"")
    clumpCoordDT$seqId <- ""
    coordDT <- data.table(
      start=coordVec,
      end=coordVec,
      seqId=""
    )
  } else {
    stop("One and only one of coordDT or coordVec must be provided ...")
  }
  clumpCoordDT <- copy(clumpCoordDT)
  setkey(clumpCoordDT,seqId,start,end)
  setkey(coordDT,seqId,start,end)

  clumpCoordDT[,distToPrev:=start-c(0,end[1:(.N-1)])]
  clumpCoordDT[,adjust:= cumsum(-distToPrev) + (1:.N) + interClumpGap*(0:(.N-1))]

  coordDT <- foverlaps(coordDT,clumpCoordDT,type="any")
  if(any(is.na(coordDT$clumpId))){
    warning("Some data points provided in coordDT fall into gaps! These will be omitted from the output:")
    print(coordDT[is.na(clumpId),.(seqId,start=i.start,end=i.end)])
    #coordDT <- coordDT[!is.na(clumpId),]
  }

  setnames(coordDT,c("i.start","i.end"),c("start_original","end_original"))
  coordDT[,start:=start_original+adjust]
  if(argGiven(coordVec)){ return(coordDT$start) }
  coordDT[,end:=end_original+adjust]
  coordDT[]
}
# d <- data.table(
#   position = as.numeric(sort(sample(1:1e5,30)))
# )[,geneId:=paste0("seq_",1:.N)][]
# d <- d[,.SD[c(1,1)],by=.(geneId)]
# d[,posType:=rep(c("start","end"),30)]
# d[posType=="end",position:=position+1000]
# plot(d$position,y=rep(1,nrow(d)),col=applyPalette(d$posType,palettePresets$mclaren2),pch=20)
# geneCoords <- dcast(d,formula=geneId~posType,value.var = "position")
# geneCoords[,seqId:="chr3H"]
#
# # Use unionCoordDT() to get the clump boundaries (this is different to how we started, listing gaps)
# clumps <- unionCoordDT(geneCoords,distCutoff = 5e3)
#
# # translateCoordsGapkill() will translate any coordDT to catenate the specified clumps together. The gap placed between them is given in 'interClumpGap='
# geneCoordsTransformed <- translateCoordsClumpGapkill(coordDT=geneCoords,clumpCoordDT=clumps,interClumpGap = 1e3)
# # use the same transformation to transform some axes
# axisRaw <- geneCoordsTransformed[,seq(from=min(start_original),to=max(end_original),by=1e3)]
# #debugonce(translateCoordsGapkill)
# axisTransformed <- translateCoordsClumpGapkill(coordVec = axisRaw,clumpCoordDT=clumps,interClumpGap = 1e3)
# clumpsTransformed <- translateCoordsClumpGapkill(coordDT=clumps,clumpCoordDT=clumps,interClumpGap = 1e3)
#
# null_plot(geneCoordsTransformed[,range(c(start,end))],0:2,xaxt="n",yaxt="n")
# axis( #bit shit ...
#   side = 1,
#   at = axisTransformed,
#   labels=axisRaw,
#   cex.axis=.7,
#   las=2
# )
# clumpsTransformed[,{
#   rect(xleft = start,xright=end,ybottom=0,ytop=2,col="#BB000066",border=F)
# }]
# geneCoordsTransformed[,{
#   lines(
#     x=c(start,end),
#     y=c(1,1),
#     lwd=2
#   )
# },by=.(geneId)]

#' @export
ambigRegex <- function(seq,conversionTable=ambig2ntAmbigOnly){
  conversionTable[,{
    seq <<- sub(nt,paste0("[",ambig,"]"),seq)
    NULL
  },by=.I]
  seq
}

#' @export
coordDTNameSeqIdChrSwitch <- function(coordDT){ coordDT[,chr:=seqId][,seqId:=name][,name:=NULL][] }


