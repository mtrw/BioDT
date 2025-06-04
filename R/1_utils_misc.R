
#' @export
combinedVar <- function(n,mu,var){
  if(length(n)<2 | length(mu)<2 | length(var)<2 | length(n)!=length(mu) | length(mu)!=length(var)){ stop("`n`, `mu`, and `var` must be equal length vectors of length > 1") }
  cVar <- 0.0
  for(i in 1:(length(n)-1)){
    cVar <- ( (n[i]-1)*var[i]+(n[i+1]-1)*var[i+1] )/( n[i]+n[i+1]-1 )  + ( n[i]*n[i+1]*(mu[i]-mu[i+1])**2 )/( (n[i]+n[i+1])*(n[i]+n[i+1]-1) )
  }
  cVar
}
# You know the means, vars, and sample sizes of some datasets. What is the var of the dataset you would make by combining them?

#' @export
combinedSd <- function(n,mu,sd){
  combinedVar(n,mu,sd**2) %>% sqrt
}

#' @export
ce <- function(...){
  cat(paste0(...,"\n"), sep='', file=stderr()) %>% eval(envir = globalenv() ) %>% invisible()
}

#' @export
`%!in%` <- function(x,y){
  !(x %in% y)
}

#' @export
argGiven <- function(x){
  !is.null(x)
}
#' @export
argNotGiven <- function(x){
  is.null(x)
}

#' Collapse a vector
#' @export
pastec <- function(x,sep=""){
  paste0(x,collapse=sep)
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

#' Use it in place of a dirname to ensure it really exists.
#' @export
existentDir <- function(dirName){
  if(!dir.exists(dirName)){
    dir.create(dirName,recursive = T)
  }
  dirName
}

#'
#' @export
fileExistsNonZeroSize <- function(fName){
  file.exists(fName) & file.size(fName)>0
}

#'
#' @export
filenameNopathNoext <- function(f,newPath="",newExt=""){
  sub("\\..*","",sub(".*[\\/]","",f)) %>% paste0(newPath,.,newExt)
}

# ... to a list with no sublists
#' @export
flattenList <- function(x) {
  do.call(c, lapply(x, function(y) if(is.list(y)) flattenList(y) else list(y)))
}
#flattenList(palettePresets)

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
    print(d.t(
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


#' @export
ldtply <- function(...){
  (ldply(...) %>% setDT)[]
}

#' Create a grid of values, calculated from every combination of two rows, one drawn from each of two data.tables.
#'
#' @param dtx A data table. [no default]
#' @param dty [optional] A data.table. If omitted, [NULL]
#' @param nameColx [optional] Name of the column in `dtx` from which the output row names are to be taken. [NULL]
#' @param nameColy [optional] Name of the column in `dty` from which the output column names are to be taken. [NULL]
#' @param FUN  [no default]
#' @returns A matrix of class whatever FUN outputs.
#' @description
#' The i,j-th entry of the output matrix is the result of `FUN(dtx[i],dty[j])`.
#' @details
#' Wraps to the console width. Coordinates are given relative to the sequence start.
#' @examples
#' ldply but a data.table returns
#' d1 <- d.t(value=1:26,label1=LETTERS,label2=letters)
#' d2 <- d.t(value=26:1,label1=LETTERS,label2=letters)
#' gridApplyDT(d1,nameColx = "label1",nameColy = "label2",FUN = function(x,y){ abs(x$value - y$value) })
#' gridApplyDT(d1,d2,nameColx = "label1",nameColy = "label2",FUN = function(x,y){ abs(x$value - y$value) })
#' @export
gridApplyDT <- function(dtx,dty=NULL,FUN,nameColx=NULL,nameColy=NULL){
  namesx <- if(argGiven(nameColx)){
    dtx[,get(nameColx)]
  } else {
    NULL
  }
  symmetrical <- if(argNotGiven(dty)){ # Symmetrical
    dty <- dtx
    namesy <- namesx
    TRUE
  } else {
    FALSE
  }
  if(argGiven(nameColy)){
    namesy <- dty[,get(nameColy)]
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



#' Collect input from a user interactively.
#'
#' @param question [optional] A question to prompt the user with [NULL]
#' @returns Whatever the user inputs before pressing <return>, as a character string.
#' @examples
#' while(ask("Guess the number I'm thinking of!")!=sample(1:1e9L,1)){ ce("WRONG!") }
#' @export
ask <- function(question="NULL"){
  if(argGiven(question)){
    ce(question)
  }
  readline()
}


#' @export
pmean2 <- function(x,y){
  (x+y)/2
}

#' @export
pmean <- function(...,na.rm=FALSE){
  m <- do.call(cbind,list(...))
  d <- ncol(m)-apply(m,1,function(r){ sum(is.na(r)) })
  rowSums(m,na.rm=na.rm)/d
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
qq <- function(x,y=1:length(x)){
  stopifnot(length(x)==length(y))
  d.t(
    qx=sort(x)/length(x),
    qy=sort(y)/length(y)
  )
}


#' @export
qqPlot <- function(x,y=1:length(x),...){
  d <- qq(x,y)
  plot(d$qx,d$qy,...)
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
  out <- rep(NA_real_,length(in_x))
  in_x[!in_x %between% range(xs)] <- NA
  for(i in which(!is.na(in_x))){
    righti <- which(cumsum(xs>=in_x[i])==1)
    if(in_x[i]==xs[righti]){ out[i] <- ys[righti]; next; }
    slope <- (ys[righti]-ys[righti-1])/(xs[righti]-xs[righti-1])
    out[i] <- ys[righti-1] + slope*(in_x[i]-xs[righti-1])
  }
  out
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
  is_seqListDT(seqList,croak=TRUE,message="Object must be a valid seqListDT.")
  l <- nchar(seqList[[1]])
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
  d.t(
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


#' @export
dcast2matrix <- function(...){
  tmp <- dcast(...)
  rnames <- tmp[[1]]
  #cnames <- colnames(tmp)[-1]
  m <- as.matrix(tmp[,2:ncol(tmp)])
  rownames(m) <- rnames
  m
}


#' @export
groupify <- function(x,nClasses=10){
  x %scale_between% c(1,nClasses) %>% round %>% as.integer
}

#' @export
groupifyByBoundaries <- function(x,boundaries,values=NULL){
  if(is.null(values)){values=0:length(boundaries)+1}
  if(length(values)!=(length(boundaries)+1)){ stop("Values must be provided for each boundary gap, and above/below the limits") }
  if(is.unsorted(boundaries)){ stop("Boundaries in ascending order please!") }
  o <- x
  o[x<first(boundaries)] <- first(values)
  o[x>=last(boundaries)] <- last(values)
  if(length(boundaries)>1){
    for(i in 1:(length(boundaries)-1)){
      o[x >= boundaries[i] & x < boundaries[i+1]] <- values[i+1]
    }
  }
  return(o)
}

#' Test all vector members are the same thing
#' @export
same <- function(x){
  nu(x)==1
}

#' Draw an arch
#' @export
arch <- function(start,end,bottom,top,n=100L,...){
  d <- d.t(
    x=cos(seq(0,pi,length.out=n)) %>% scale_between(start,end),
    y=sin(seq(0,pi,length.out=n)) %>% scale_between(bottom,top)
  )
  lines(d$x,d$y,...)
}
# null_plot(-10:10,-10:10)
# arch(-5,2,-4,4,col="red",lwd=4,lty=2)

#' The largest value found in a vector so far
#' @export
maxSoFar <- function(x){
  #rewrite in Rcpp
  mxsf <- x[1]
  sapply(x,function(xi) {if(xi>mxsf){ mxsf<<-xi }; mxsf } )
}

#' The smallest value found in a vector so far
#' @export
minSoFar <- function(x){
  #rewrite in Rcpp
  mnsf <- x[1]
  sapply(x,function(xi) {if(xi<mnsf){ mnsf<<-xi }; mnsf } )
}

#' Return the argument (a directory name) but also make sure a directory called that exists.
#' @export
existentDir <- function(dirName){
  sapply(dirName,function(dn){
    if(!dir.exists(dn)){
      dir.create(dn,recursive = T)
    }
    dn
  })
}

#'
#' @export
mcldtply <- function(...){
  require(parallel)
  mclapply(...) %>% rbindlist
}

#' Expand grid wrapper that returns a data.table
#' @export
expandGridDt <- function(...){
  expand.grid(...,stringsAsFactors = FALSE) %>% setDT
}

#' Add a tiny bit to each number.
#' @description
#' Useful when a collection of numbers on [0,Inf) need to be logged
#' @param x The vector to add a bit to [no default]
#' @param inc How much to add [0.001 times the smallest non-zero gap between two consecutive data points (after sorting) ]
#' @param shiftToZero Shift all values to sit along \[0,...\] before adding `inc` [FALSE]
#' @returns A vector
#' @export
plusAtinyBit <- function(x,inc={t<-abs(diff(sort(x[!is.na(x)]))); min(t[t>0])}*0.001,shiftToZero=FALSE){
  if(shiftToZero){ x <- x-min(x,na.rm=T) }
  x+inc
}


#' Return the range of x, extended by a bit
#' @param x The vector to take the range of
#' @param m Subtract me from the low end
#' @param p Add me to the high end [m]
#' @returns A length two vector
#' @export
`%plusMinus%` <- function(x,m,p=m){
  range(x) + c(-m,p)
}



# Check must work for from to ranges going in different directions.
# #' Create a function to project values from one range to another
# #' @param fromRange The input values of the created function will be assumed to be given on this range.
# #' @param toRange The created function
# #' @returns A function of one argument ("values"), which given any value will transform it linearly from `fromRange` to `toRange`.
# #' @example
# #' tfm <- makeRangeTransformer(c(0,10),c(100,0))
# #' tfm(1:9) # returns (90, 80, 70, 60, 50, 40, 30, 20, 10)
# #' @export
# makeRangeTransformer <- function(fromRange,toRange){
#   fromRange <- range(fromRange)
#   toRange <- range(toRange)
#   function(values){
#     (((values - fromRange[1]) / abs(diff(fromRange))) * abs(diff(toRange))) + toRange[1]
#   }
# }






