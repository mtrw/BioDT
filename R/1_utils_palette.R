

#' @export
showPalettes <- function(colPalettes=palettePresets,gradientN=NULL){
  colPalettes <- flattenList(colPalettes)
  if(is.character(colPalettes)){
    colPalettes <- list(colPalettes)
  }
  span <- sapply(colPalettes,length)%>%max
  bdr <- if(argGiven(gradientN)){NA}else{"#000000FF"}
  null_plot(0:1,c(0,length(colPalettes)),yaxt="n",xaxt="n")
  l_ply(seq_along(colPalettes), function(j) {
    #dev j <- 1; i <- 1
    if(argNotGiven(gradientN)){
      colPalette <- colPalettes[[j]]
      n <- length(colPalette)
    } else {
      colPalette <- makePalette(colPalettes[[j]],gradientN)
      n <- gradientN
    }
    width <- 1/n
    l_ply(1:n,function(i){
      rect(
        (i-1)*width,(j-0.7),
        i*width    ,(j-0.3),
        border=bdr,
        col=colPalette[i]
      )
    })
    text(0,j,labels=names(colPalettes)[j],pos=4)
  })
}

parseColChain <- function(colChain,setAlpha=NULL){
  # Parse and wrangle chain to rgb
  chain <- col2rgb(colChain,a=T)/255

  # Set alpha as requested
  if(argGiven(setAlpha)){
    if(!all(setAlpha %between% 0:1)){ stop("'setAlpha' values should be scaled between 0 and 1") }
    chain[4,] <- if( length(setAlpha)==1 | length(setAlpha)==length(colChain) ){
      setAlpha
    } else {
      approx((1:length(setAlpha)) %scale_between% c(1,length(colChain)),setAlpha,1:length(colChain))$y
    }
  }

  return(chain)
}

#' Makes a simple linear palette interpolating the colChain across n evenly spaced points. Mostly useful for calling from other funs.
#' @export
makePalette <- function(colChain=palettePresets$wheel$wheel,n=NULL,at=NULL,setAlpha=NULL){ # colours in the palette defined by the colChain, matched to [`type=`] "discrete" or "continuous" data.
  if(sum(argGiven(n),argGiven(at))!=1){ stop("One of 'n' or 'at' must be supplied, not both.") }
  if(argGiven(n)){
    at <- 1:n
  }
  parseColChain(colChain, setAlpha) %>%
    apply( 1 , function(r) { approx(seq_along(r) %scale_between% range(at),r,at)$y } ) %>%
    apply( 1 , function(r) { rgb(r[1],r[2],r[3],r[4]) } )
}

#' @export
alpha <- function(colChain,setAlpha=1L){
  parseColChain(colChain,setAlpha) %>%
    t %>%
    apply( 1 , function(r) {rgb(r[1],r[2],r[3],r[4])} )
}



#' @export
applyPalette <- function(x,colChain=palettePresets$wheel$wheel,setAlpha=NULL,discreteOrContinuous=c("guess","discrete","continuous"),show=FALSE){ # An evenly spaced palette interpolating the colChain

  # Ascertain discreteness
  bi <- isBehaved(x)
  #n <- length(x)
  discrete <- if(discreteOrContinuous[1]=="guess"){
    !(is.numeric(x[bi]) | is.integer(x[bi]))
  } else if (discreteOrContinuous[1]=="discrete"){
    TRUE
  } else if (discreteOrContinuous[1]=="continuous"){
    FALSE
  } else {
    stop("Value for `discreteOrContinuous=` argument must be \"discrete\" or \"continuous\" or \"guess\" ...")
  }



  if(discrete==TRUE){
    # Build table per value
    tbl <- data.table(joiner=unique(x[bi]))[,col:=makePalette(colChain,n=.N,setAlpha=setAlpha)][]
    setkey(tbl,joiner)
    return( tbl[data.table(joiner=x[bi]),on=.(joiner)]$col )
  } else {
    # Check low number of unique vals -- if so, make a table using interpolation and merge
    if(nu(x[bi])/length(x[bi]) < 0.20){
      #x <- sample(5,10,r=T)
      tbl <- data.table(joiner=unique(x[bi]))[,col:=makePalette(colChain,at=joiner,setAlpha=setAlpha)][]
      setkey(tbl,joiner)
      return( tbl[data.table(joiner=x[bi]),on=.(joiner)]$col )
    }
    # ... else ...
    # Interpolate
    x[bi] <- makePalette(colChain=colChain,at=x[bi],setAlpha=setAlpha)
    x[!bi] <- NA_character_
    return( x )#makePalette(colChain=colChain,at=x,setAlpha=setAlpha) )
  }
  # if(returnLegend==TRUE){ ... } # wherever appropriate ...
}
