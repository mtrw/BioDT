

#' @export
showPalettes <- function(colPalettes=palettePresets,gradientN=NULL){
  if(is.list(colPalettes)){ colPalettes <- flattenList(colPalettes) }
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
applyPalette <- function(x,colChain=palettePresets$wheel$wheel,setAlpha=NULL,discreteOrContinuous=c("guess","discrete","continuous"),show=FALSE,returnLegend=FALSE,assignLegend=NULL){ # An evenly spaced palette interpolating the colChain

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
    tbl <- d.t(joiner=unique(x[bi]))[,col:=makePalette(colChain,n=.N,setAlpha=setAlpha)][]
    setkey(tbl,joiner)
    out <- tbl[d.t(joiner=x[bi]),on=.(joiner)]$col
    if(returnLegend==TRUE){ return( list(out,tbl[,.(label=joiner,col)]) ) }
    if(argGiven(assignLegend)){ assign(assignLegend,tbl[,.(label=joiner,col)],envir=globalenv()) }
    return( out )
  } else {
    # Check low number of unique vals -- if so, make a table using interpolation and merge
    if(nu(x[bi])/length(x[bi]) < 0.20){
      tbl <- d.t(joiner=unique(x[bi]))[,col:=makePalette(colChain,at=joiner,setAlpha=setAlpha)][]
      if( returnLegend==TRUE || argGiven(assignLegend) ) {
        leg <- tbl[,.(label=if(.N>1){c(first(joiner),rep(NA,.N-2),last(joiner))}else{joiner},col)]
        if( returnLegend==TRUE ){ return( list(out,leg) ) }
        if( argGiven(assignLegend) ) { assign(assignLegend,leg,envir=globalenv()) }
      }
      setkey(tbl,joiner)
      return( tbl[d.t(joiner=x[bi]),on=.(joiner)]$col )
    } else {
      # Interpolate
      p <- character(length(x))
      p[bi] <- makePalette(colChain=colChain,at=x[bi],setAlpha=setAlpha)
      p[!bi] <- NA_character_
      if( returnLegend==TRUE || argGiven(assignLegend) ) {
        leg <- d.t( col = p[bi] )[order(x[bi]),][,label:=if(.N>1){c(min(x[bi]),rep(NA,.N-2),max(x[bi]))}else{x[bi]}][]
        if( returnLegend==TRUE ){ return( list(out,leg) ) }
        if( argGiven(assignLegend) ) { assign(assignLegend,leg,envir=globalenv()) }
      }
      return( p )
    }
  }
}



#' @export
plotLegend <- function( xRange , yRange , legendDT , col_bg="#FFFFFFFF" , col_border="#000000FF" , col_text=legendDT$col , pos_middle=0.5 , gap=NULL , ...){

  # applyPalette(strsplit("Of the world’s 20 largest economies, Australia is the only one not using nuclear energy, or moving towards using it."," ")[[1]],assignLegend = "l")
  # applyPalette(sample(22),assignLegend = "l")
  # legendDT = l
  # col_text=legendDT$col
  # col_border="#000000FF"
  # pos_middle=0.5
  # xRange <- c(1,5)
  # yRange <- c(1,5)
  # gap=NULL
  # #plot(1:nrow(legendDT),col=legendDT$col,cex=3,pch=20)

  if(argNotGiven(gap)){
    gap <- if(any(is.na(legendDT))){ 0.0 } else { 0.2 }
  }
  legendDT$col <- col_text
  nCols <- nrow(legendDT)
  xInt <- diff(xRange)
  yInt <- diff(yRange)
  rect(xRange[1],yRange[1],xRange[2],yRange[2],border=col_border,col=col_bg)
  colBoxLeft  <-  xRange[1]+(xInt*pos_middle)
  colBoxRight <-  xRange[1]+(xInt*(29/30))
  colBoxTop <-    yRange[2]-(yInt*0.05)
  colBoxBottom <- yRange[1]+(yInt*0.05)
  textEndRight <- xRange[1]+(xInt*pos_middle)
  #colBoxWidth <- diff(c(colBoxRight,colBoxLeft))
  colBoxHeight <-  diff(c(colBoxBottom,colBoxTop))
  eachColHeight <-colBoxHeight/(nCols+(nCols-1)*gap)
  eachGapHeight <- gap * eachColHeight
  eachColAndGapHeight <- eachColHeight + eachGapHeight
  colTops <- colBoxTop - ((1:nCols)-1)*eachColAndGapHeight
  l_ply(1:nCols,function(i_col){
    #i_col <- 1
    rect(colBoxLeft,colTops[i_col]-eachColHeight,colBoxRight,colTops[i_col],border=NA,col=legendDT$col[i_col])
    text(textEndRight,colTops[i_col]-eachColHeight/2,labels=legendDT$label[i_col],pos=2)#,...)
  })


}


