#' @export
d.t <- function(...){data.table(...)}

#' @export
ce <- function(...){
  cat(paste0(...,"\n"), sep='', file=stderr()) %>% eval(envir = globalenv() ) %>% invisible()
}

#' @export
argGiven <- function(x){
  !is.null(x)
}
#' @export
argNotGiven <- function(x){
  is.null(x)
}

#' @export
printVecRows <- function(x){cat(paste(x,collapse="\n"))}


#' @export
isBehaved <- function(x){
  !(is.na(x) | is.null(x) | is.infinite(x) | is.nan(x))
}

#' @export
pd <- function(x,add=F,bw="nrd0",plotMain=NA,...){
  if(!add){
    x %>% density(na.rm=TRUE,bw=bw) %>% plot(main=plotMain,...)
  } else {
    x %>% density(na.rm=TRUE,bw=bw) %>% lines(main=plotMain,...)
  }
}



#' @export
filenameNopathNoext <- function(f,newPath="",newExt=""){
  sub("\\..*","",sub(".*[\\/]","",f)) %>% paste0(newPath,.,newExt)
}

#create an empty plot with ranges x=c(low,high) and y=ditto
#' @export
null_plot <- function(x,y,xlab=NA,ylab=NA,revx=F,revy=F,...){
  xl<-range(x,na.rm=T)
  yl<-range(y,na.rm=T)
  if(revx==T){ xl <- rev(xl) }
  if(revy==T){ yl <- rev(yl) }
  plot(NULL,xlim=xl,ylim=yl,xlab=xlab,ylab=ylab,...)
}

# For a list of lists, get item sliceIdx (name or number) from each, return as vector
#' @export
sliceList <- function(list,sliceIdx){
  sapply(list,function(li){li[[sliceIdx]]})
}

#' @export
wait <- function(message="Press [enter] to continue"){
  invisible(readline(prompt=message))
}

#simultaneously change the names of things to a particular thing if they match (EXACTLY!) a particular string.
#swap(vec=iris$Species,matches = c("virginica"),names = c("XXXXX"))
#' @export
swap <- function(vec,matches,names,na.replacement=NA){
  orig_vec <- copy(vec)
  #if(sum(! matches %in% names ) > 0 ) { stop("Couldn't find all matches in names") }
  if(length(matches) != length(names)) { stop("Lengths of `matches` and `names` vectors don't match, you old bison!") }
  if(is.factor(vec)) { levels(vec) <- c(levels(vec),names,na.replacement) }
  vec[is.na(orig_vec)] <- na.replacement
  l_ply( 1:length(matches) , function(n){
    vec[orig_vec==matches[n]] <<- names[n]
  })
  vec
}


# Test a blast DB exists
# `root` is typically (e.g.) "my_subject.fasta"
# `typePrefixes` are ("n")ucleotide and ("p")rotein
#' @export
blastDBexists <- function(root,typePrefix="n"){
  expectedFiles <- paste0(root,".",typePrefix,c("db")) # This could add files but I can't yet find a good ref to describe what the minimum necessary set of files to constitute a blastdb is
  if(all(file.exists(expectedFiles))){
    return(TRUE)
  } else {
    ce("Database not detected:")
    print(data.table(
      expectedFile=expectedFiles,
      exists=file.exists(expectedFiles)
    ))
    return(FALSE)
  }
}
#blastDBexists("/data/gpfs/projects/punim1869/shared_data/refs/morex_v1/assembly/morex_v1_pseudomolecules.fasta")

# Number of unique things
#' @export
nu <- function(x){
  sum(!duplicated(x))
}


# ldply but a data.table returns
#' @export
ldtply <- function(...){
  (ldply(...) %>% setDT)[]
}


#' @export
gridApplyDT <- function(dtx,dty=NULL,FUN,nameColx=NULL,nameColy=NULL,symmetrical=FALSE){
  namesx <- if(!is.null(nameColx)){
    dtx[,get(nameColx)]
  } else {
    NULL
  }
  if(symmetrical==TRUE){
    if(!is.null(dty) | !is.null(nameColy)){
      stop("With `symmetrical`==TRUE, only `dtx` and `nameColx` arguments should be provided")
    }
    dty <- dtx
    namesy <- namesx
  } else { # Asymmetrical
    namesy <- if(!is.null(nameColy)){
      dty[,get(nameColy)]
    } else {
      NULL
    }
  }
  xlen <- nrow(dtx)
  ylen <- nrow(dty)
  om <- matrix(as(NA,class(FUN(dtx[1,],dty[2,]))),nrow=xlen,ncol=ylen,dimnames=list(namesx,namesy))
  l_ply(1:xlen,function(i_x){
    l_ply(1:if(symmetrical==TRUE){i_x}else{ylen},function(i_y){
      #ce("ix:",i_x,"\tiy:",i_y)
      om[i_x,i_y] <<- FUN(dtx[i_x,],dty[i_y])
    })
  })
  if(symmetrical==TRUE){
    om[upper.tri(om)] <- om[lower.tri(om)]
  }
  om
}


#' @export
pmean2 <- function(x,y){
  x+y/2
}

#' @export
mostCommonThing <- function(x,threshold_prop=0,na.rm=T,draw_out=NULL,na_wins_out=NA,threshold_notmet_out=NA){
  #browser()
  if(all(is.na(x))){
    return(as(NA,class(x)))
  }
  if(na.rm){
    x <- x[!is.na(x)]
  }
  if(length(x)==0){
    stop("Length of x is zero")
  }
  x[is.na(x)] <- "NA"
  tbl <- table(x)
  Ma <- names(tbl[order(-tbl)])[1]

  if (length(tbl)==1){
    as(Ma,class(x))
  }
  else if (tbl[order(-tbl)][1]==tbl[order(-tbl)][2]){
    if(is.null(draw_out)){
      sort(c(names(tbl[order(-tbl)][1]),names(tbl[order(-tbl)][2])))[1]
    } else {
      as(draw_out,class(x))
    }
  } else if (Ma=="NA"){
    as(na_wins_out,class(x))
  } else if (tbl[order(-tbl)][1]/sum(tbl) < threshold_prop){
    as(threshold_notmet_out,class(x))
  } else {
    as(Ma,class(x))
  }
}
#c(2,1,1,2) %>% most_common_thing()
#most_common_thing(x=c("G",NA,"T"),draw_out = NA)
#most_common_thing(x=c(1,1,1,2,2,2,3,3,NA,NA,NA,NA,NA),draw_out="DRAW!",na_wins_out="na was the most common",na.rm=T)


#' @export
MSA <- function(seqDT,method="ClustalOmega",...){
  require(msa)
  require(Biostrings)
  is_seqDT(seqDT,objName = deparse(substitute(seqDT)))


  seqDT[,alnSeq:={
    seq %>%
      Biostrings::DNAStringSet() %>%
      msa::msa(method=method,order="input",...) %>%
      msaConvert(type = "seqinr::alignment") %>%
      `[[`("seq")
  }][]
}

#' @export
drawConsoleLine <- function(){
  require(stringi)
  w <- options()$width
  ce(paste0(rep("-",w),collapse=""))
}


#scale a list of values to between two points, proportionally spaced as they were originally
#rnorm(100) %>% scale_between(20,29) %>% pd
#' @export
scale_between <- function(x,lower,upper){
  if(all(x==mean(x,na.rm=T))) return(rep(mean(c(lower,upper),na.rm=T),length(x)))
  ( x - min(x,na.rm=T) ) / (max(x,na.rm=T)-min(x,na.rm=T)) * (upper-lower) + lower
}


#An infix wrapper for the above
#' @export
`%scale_between%` <- function(x,y){
  x %>% scale_between(y[1],y[2])
}


#' @export
interpolate <- function(in_x,xs,ys){
  sapply(in_x,function(ix){
    #browser()
    righti <- which(cumsum(xs>=ix)==1)
    if(ix==xs[righti]){return(ys[righti])}
    slope <- (ys[righti]-ys[righti-1])/(xs[righti]-xs[righti-1])
    ys[righti-1] + slope*(ix-xs[righti-1])
  })
}



#give it "sv" as AG or CT or AGT etc, must be sorted alphabetically.
#' @export
IUPAC <- function(sv){
  iupac <- c( "A", "G", "C", "T", "R" ,  "Y" ,  "S" ,  "W" ,  "K" ,  "M" ,  "B" ,   "D" ,   "H" ,   "V" ,   "N" ,   "a", "g", "c", "t", "r" ,  "y" ,  "s" ,  "w" ,  "k" ,  "m" ,  "b" ,   "d" ,   "h" ,   "v" ,   "n"    , "-")
  nucs  <- c( "A", "G", "C", "T", "AG" , "CT" , "CG" , "AT" , "GT" , "AC" , "CGT" , "AGT" , "ACT" , "ACG" , "ACGT", "a", "g", "c", "t", "ag" , "ct" , "cg" , "at" , "gt" , "ac" , "cgt" , "agt" , "act" , "acg" , "acgt" , "-")
  swap(sv,nucs,iupac)
}
#IUPAC(c("A","AT","CT","-","ACGT"))

# DELETE SOON
#' @export
getAlnCol <- function(seqList,pos){
  o <- vector(length(seqList),mode="character")
  plyr::l_ply(seq_along(seqList),function(i){
    o[i] <<- substr(seqList[[i]],pos,pos)
  })
  o
}


# REINSTATE TO ggetAlnCol SOON
#' @export
getAlnvCol <- function(seqList,pos){
  o <- vector(length(seqList),mode="character")
  plyr::l_ply(seq_along(seqList),function(i){
    o[i] <<- substr(seqList[i],pos,pos)
  })
  o
}



#' @export
consensus <- function(seqList){
  l <- stri_length(seqList[[1]])
  s <- paste0(rep(" ",l),collapse = "")
  for(i in 1:l){ #applyify
    #dev i <- 7
    col <- getAlnCol(seqList,i)
    #ce(i)
    if(mostCommonThing(col)=="-"){
      substr(s,i,i) <- "-"
    } else {
      #browser()
      substr(s,i,i) <- IUPAC(col[col!="-"] %>% unique %>% sort %>% paste(collapse = ""))
    }
    if(is.na(s)){browser()}
  }
  s
}


#' @export
consensusSeq <- function(alnSeqDT,seqId="unnamed_consensus_sequence"){
  if(!all(stringi::stri_length(alnSeqDT$alnSeq)==stringi::stri_length(alnSeqDT$alnSeq[1]))){stop("All entries in `alnSeq` column must be the same length.")}
  l <- stringi::stri_length(alnSeqDT$alnSeq[1])
  s <- paste0(rep(" ",l),collapse = "")
  for(i in 1:l){ #applyify
    #dev i <- 7
    col <- getAlnvCol(alnSeqDT$alnSeq,i)
    #ce(i)
    if(mostCommonThing(col)=="-"){
      substr(s,i,i) <- "-"
    } else {
      #browser()
      substr(s,i,i) <- IUPAC(col[col!="-"] %>% unique %>% sort %>% paste(collapse = ""))
    }
    if(is.na(s)){browser()}
  }
  data.table(
    seqId = seqId,
    seq=s
  )
}

#' @export
franktd <- function(...){ frank(...,ties.method="dense") }


#' @export
postpadChar <- function(x,len=2,pad="0"){
  paste0(x,substr(rep(paste0(rep(pad,len),collapse=""),length(x)),rep(0,length(x)),(len-stringi::stri_length(x))))
}

#' @export
prepadChar <- function(x,len=2,pad="0"){
  paste0(substr(rep(paste0(rep(pad,len),collapse=""),length(x)),rep(0,length(x)),(len-stringi::stri_length(x))),x)
}
