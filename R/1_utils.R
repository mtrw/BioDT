#' @export
ce <- function(...){
  cat(paste0(...,"\n"), sep='', file=stderr()) %>% eval(envir = globalenv() ) %>% invisible()
}

#' @export
isBehaved <- function(x){
  !is.na(x) & !is.nan(x) & !is.infinite(x)
}

#' @export
isBehavedPositive <- function(x){
  isBehaved(x) & x>0
}

#' @export
pstop <- function(...){
  stop(paste0(...))
}

#' @export
trueFun <- function(x){TRUE}

#' @export
makeValidator <- function(validatorSpecs,wholeObjectValidatorFun=NULL){
  vS <- copy(validatorSpecs)
  #valspecs has three cols: colName(character), colClasses(character), colValidatorFuns(list(function(x),function(x),...))
  function(validateMe,error=T){
    if(!is.data.table(validateMe))                                                { if(error==TRUE){ pstop("object must be a data.table"); return(FALSE) }                                                                        else {return(FALSE)} }
    sapply(1:nrow(vS),function(i){
      if(! vS[i]$colName %in% colnames(validateMe))                               { if(error==TRUE){ pstop("Column ",vS[i]$colName," not found"); return(FALSE) }                                                                 else {return(FALSE)} }
      if(! class(validateMe[,get(vS[i]$colName)])==vS[i]$colClass)                { if(error==TRUE){ pstop("Column ",vS[i]$colName," must be of class ",vS[i]$colClass); return(FALSE) }                                          else {return(FALSE)} }
      if(! all(vS[i]$colValidatorFun[[1]](validateMe[,get(vS[i]$colName)])))      { if(error==TRUE){ pstop("Column ",vS[i]$colName," failed validation function ",deparse(substitute(vS[i]$colValidatorFun))); return(FALSE) }    else {return(FALSE)} }
      if(!is.null(wholeObjectValidatorFun)){
        if(!wholeObjectValidatorFun(validateMe))                                  { if(error==TRUE){ pstop("Column ",vS[i]$colName," failed validation function ",deparse(substitute(wholeObjectValidatorFun))); return(FALSE) }    else {return(FALSE)} }
      }
      TRUE
    }) %>% all
  }
}

#' @export
containsWhitespace <- function(x){
  all(grepl("[[:space:]]",x))
}

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

#' @export
filename_nopath_noext <- function(f,newPath="",newExt=""){
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

# Reverse complement a text DNA sequence.
# Does not yet do all IUPAC codes.
#rc(c("AAAAAC","AgCtaGcTxxxx----agcT"))
#' @export
rc <- function(x){
  stringi::stri_reverse(x) %>% stringi::stri_trans_char("acgtACGT","tgcaTGCA")
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

# ldply but a data.table returns
#' @export
ldtply <- function(...){
  ldply(...) %>% setDT
}




#' @export
gridApplyDT <- function(dtx,dty=NULL,distFun,nameColx=NULL,nameColy=NULL,symmetrical=FALSE){
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
  om <- matrix(as(NA,class(distFun(dtx[1,],dty[2,]))),nrow=xlen,ncol=ylen,dimnames=list(namesx,namesy))
  l_ply(1:xlen,function(i_x){
    l_ply(1:if(symmetrical==TRUE){i_x}else{ylen},function(i_y){
      #ce("ix:",i_x,"\tiy:",i_y)
      om[i_x,i_y] <<- distFun(dtx[i_x,],dty[i_y])
    })
  })
  if(symmetrical==TRUE){
    om[upper.tri(om)] <- om[lower.tri(om)]
  }
  om
}

#' @export
alnGetClumps <- function(alnDT,dist,distCutoff=1e3L,hclustAlgorithm="single"){ # CHANGE TO COORDS GET CLUMPS AND BUILD IN aln2coords CONVERSION!!!
#alnDT <- alnRph12toPgF; distCutoff <- 1e4; nameColx = "hitId"; hclustAlgorithm="single"
  alnDT <- copy(alnDT)
  clumps <- alnDT[,{
    #debugonce(gridApplyDT)
    dm <- gridApplyDT(dtx=.SD,distFun=function(x,y){ min( abs(y$sStart-x$sEnd) , abs(x$sStart-y$sEnd) ) },symmetrical=TRUE) %>% as.dist()
    .(
      clump = hclust(dm,method=hclustAlgorithm) %>% cutree(h=distCutoff),
      sStart,
      sEnd
    )
  },by=.(sSeqId)]
  clumps[,.(span=diff(range(c(sStart,sEnd))),nAlns=.N,start=min(c(sStart,sEnd)),end=max(c(sStart,sEnd))),by=.(sSeqId,clump)]
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
msa <- function(seqs,method="ClustalOmega",...){
  require(msa)
  require(Biostrings)

  seqs %>%
  Biostrings::DNAStringSet() %>%
  msa::msa(method=method,order="input",...) %>%
  msaConvert(type = "seqinr::alignment") %>%
  `[[`("seq")
}

#' @export
drawConsoleLine <- function(){
  require(stringi)
  w <- options()$width
  ce(paste0(rep("-",w),collapse=""))
}


#' @export
colHexToDec <- function(col){
  sapply(col,function(c){
    strtoi(
      paste0("0x",c(
        substr(c,2,3),
        substr(c,4,5),
        substr(c,6,7),
        substr(c,8,9)
      )
      )
    )
  }) %>% t
}


#' @export
prepadChar <- function(x,len=2,pad="0"){
  paste0(substr(rep("00",length(x)),rep(0,length(x)),(2-stringi::stri_length(x))),x)
}



#' @export
colDecToHex <- function(col){
  apply(col,1,function(c){
    paste0("#",paste0(as.hexmode(c) %>% prepadChar(2,"0"),collapse=""))
  })
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





#' @export
timPalette <- function(colChain=c("#EEDD00FF","#009933FF"),n=100L){
  colHexToDec(colChain) %>%
    apply( . , 2 , function(c) interpolate( (1:n) %scale_between% c( 1 , length(colChain)) , 1:length(c) , c ) ) %>%
    round %>%
    colDecToHex()
}








#give it "sv" as AG or CT or AGT etc, must be sorted alphabetically.
#' @export
IUPAC <- function(sv){
  iupac <- c( "A", "G", "C", "T", "R" ,  "Y" ,  "S" ,  "W" ,  "K" ,  "M" ,  "B" ,   "D" ,   "H" ,   "V" ,   "N" ,   "a", "g", "c", "t", "r" ,  "y" ,  "s" ,  "w" ,  "k" ,  "m" ,  "b" ,   "d" ,   "h" ,   "v" ,   "n"    , "-")
  nucs  <- c( "A", "G", "C", "T", "AG" , "CT" , "CG" , "AT" , "GT" , "AC" , "CGT" , "AGT" , "ACT" , "ACG" , "ACGT", "a", "g", "c", "t", "ag" , "ct" , "cg" , "at" , "gt" , "ac" , "cgt" , "agt" , "act" , "acg" , "acgt" , "-")
  swap(sv,nucs,iupac)
}
#IUPAC(c("A","AT","CT","-","ACGT"))


#' @export
getAlnCol <- function(seqList,pos){
  o <- vector(length(seqList),mode="character")
  plyr::l_ply(seq_along(seqList),function(i){
    o[i] <<- substr(seqList[[i]],pos,pos)
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
