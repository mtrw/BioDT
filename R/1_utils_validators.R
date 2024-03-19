
#' @export
makeValidator <- function(validatorSpecs,wholeObjectValidatorFun=NULL){

  # dev
  vS <- copy(validatorSpecs)[,idx:=1:.N][]
    # dev vS <- hasStartEnd_valSpecs[,idx:=1:.N][]; wholeObjectValidatorFun<-startBeforeEnd; validateMe <- data.table(start=c(1L:5), end=c(1L:5) )
  #valspecs has three cols: colName(character), colClasses(character), colValidatorFuns(list(function(x),function(x),...))
  function(validateMe,error=T){
    if( !is.data.table(validateMe) )                                                { if(error==TRUE){ pstop("object must be a data.table"); return(FALSE) }                                                                        else {return(FALSE)} }
    if( is.null(vS$requisite) ) { vS[,requisite:=TRUE] }
    vS[, c("exists","passClass","passColValidation"):=.(
      colName %in% colnames(validateMe),
      validateMe[,class(get(colName))]!=colClass,
      { colValidatorFun[[1]](validateMe[,get(colName)]) }
    ), by=.(idx) ]

    # Whole object validaton
    if( !is.null(wholeObjectValidatorFun) ){
      if( !wholeObjectValidatorFun(validateMe) )                                    { if(error==TRUE){ pstop("Column ",vS[i]$colName," failed validation function ",deparse(substitute(wholeObjectValidatorFun))); return(FALSE) }    else {return(FALSE)} }
    }

    sapply(1:nrow(vS),function(i){
      if(! vS[i]$colName %in% colnames(validateMe))                               { if(error==TRUE){ pstop("Column ",vS[i]$colName," not found"); return(FALSE) }                                                                 else {return(FALSE)} }
      if(! class(validateMe[,get(vS[i]$colName)])==vS[i]$colClass)                { if(error==TRUE){ pstop("Column ",vS[i]$colName," must be of class ",vS[i]$colClass); return(FALSE) }                                          else {return(FALSE)} }
      if(! all(vS[i]$colValidatorFun[[1]](validateMe[,get(vS[i]$colName)])))      { if(error==TRUE){ pstop("Column ",vS[i]$colName," failed validation function ",deparse(substitute(vS[i]$colValidatorFun))); return(FALSE) }    else {return(FALSE)} }

      TRUE
    }) %>% all
  }
}






#' @export
trueFun <- function(x){TRUE}


#' @export
isBehaved <- function(x){
  all(!is.na(x) & !is.nan(x) & !is.infinite(x))
}

#' @export
isBehavedPositive <- function(x){
  all(isBehaved(x) & x>0)
}

#' @export
containsWhitespace <- function(x){
  all(grepl("[[:space:]]",x))
}

#' @export
startBeforeEnd <- function(x){
  x[,all(start<=end)]
}
