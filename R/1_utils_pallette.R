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
colDecToHex <- function(col){
  apply(col,1,function(c){
    paste0("#",paste0(as.hexmode(c) %>% prepadChar(2,"0"),collapse=""))
  })
}


#' @export
makePalette <- function(colChain=palettePresets$wheel,n=100L,setAlpha="ff",show=FALSE){ # An evenly spaced palette interpolating the colChain
  colChain <- postpadChar(colChain,9,setAlpha)
  c <- colHexToDec(colChain) %>%
    apply( . , 2 , function(c) interpolate( (1:n) %scale_between% c( 1 , length(colChain)) , 1:length(c) , c ) ) %>%
    round %>%
    colDecToHex()
  if(show==TRUE) { showPalettes(c,n) }
  c
}

#' @export
applyPalette <- function(x,colChain=palettePresets$wheel,type="guess",show=F,alpha=NULL){ # colours in the palette defined by the colChain, matched to [`type=`] "discrete" or "continuous" data.
  discrete <- if(type=="guess"){
    is.logical(x) | is.character(x)
  } else if (type=="discrete") {
    TRUE
  } else if (type=="continuous") {
    FALSE
  } else {
    stop("Value for `type=` argument must be \"discrete\" or \"continuous\" or \"guess\" ...")
  }

  isBehx <- isBehaved(x)
  xc <- x[isBehx]#x, cleaned

  if(discrete==TRUE){
    xc <- frank(xc,ties.method = "dense")
  }

  colChain <- postpadChar(colChain,9,"ff")
  c <- colChain %>%
    colHexToDec %>%
    apply( . , 2 , function(c) interpolate( xc %scale_between% c( 1 , length(c)) , 1:length(c) , c ) ) %>%
    round %>%
    colDecToHex()

  if(argGiven(alpha)){ c <- alpha(c,alpha) }

  if(show==TRUE){
    if(discrete==TRUE){
      showPalettes(c[order(xc)] %>% unique)
    } else {
      showPalettes(colChain,gradientN=if(discrete==TRUE){NULL}else{100L})
    }
  }

  x[isBehx] <- c
  x
}

#' @export
showPalettes <- function(colPalettes=palettePresets,gradientN=NULL){

  if(is.character(colPalettes)){
    colPalettes <- list(colPalettes)
  }
  span <- sapply(colPalettes,length)%>%max
  bdr <- if(argGiven(gradientN)){NA}else{"#000000FF"}
  null_plot(0:1,c(0,length(colPalettes)),yaxt="n",xaxt="n")
  l_ply(seq_along(colPalettes), function(j) {
    #dev j <- 1; i <- 2
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

#' @export
alpha <- function(colChain,setAlpha="ff"){
  if(length(setAlpha)==1){
    out <- paste0(substr(colChain,1,7),setAlpha)
  } else if (length(setAlpha)>1){
    paste0( substr(colChain,1,7) , (interpolate(1:length(colChain) %scale_between% c(1,length(setAlpha)),xs=1:length(setAlpha),ys=paste0("0x",setAlpha) %>% strtoi) %>% round %>% as.hexmode())  )
  }
}
