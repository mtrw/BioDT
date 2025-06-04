
#' @export
vTrueFun <- function(x=NULL){TRUE}

is.behaved <- function(x){
  !is.na(x) & !is.nan(x) & !is.infinite(x)
}

#' @export
vIsBehaved <- function(x=NULL){
  if(argNotGiven(x)){return("Must not have NA, NAN, or Infinite values")}
  all(is.behaved(x))
}

#' @export
vIsBehavedSameLengthStrings <- function(x=NULL){
  all(is.behaved(x)) & same(nchar(x))
}

#' @export
vIsBehavedIntlikeFloat <- function(x=NULL){
  if(argNotGiven(x)){return("Must not have NA, NAN, or Infinite values, and have class \"numeric\", and be whole numbers (i.e. no decimals)")}
  all(is.behaved(x) & (x%%1==0))
}

#' @export
vIsBehavedPositiveIntlikeFloat <- function(x=NULL){
  if(argNotGiven(x)){return("Must not have NA, NAN, or Infinite values, and have class \"numeric\", and be positive whole numbers (i.e. no decimals)")}
  all(is.behaved(x) & (x%%1==0) & x>0)
}

#' @export
vIsBehavedPositive <- function(x=NULL){
  if(argNotGiven(x)){return("Must not have NA, NAN, or Infinite values, and all values must be >0")}
  all(is.behaved(x) & x>0)
}

#' @export
vIsBehavedCommaSepInts <- function(x=NULL){
  if(argNotGiven(x)){return("Must not have NA, NAN, or Infinite values, and all values must be a digit [0-9] or a comma (,)")}
  all(is.behaved(x) & !grepl("[^[:digit:],]",x))
}

#' @export
vIsBehavedSimpleText <- function(x=NULL){
  if(argNotGiven(x)){return("Must not have NA, NAN, or Infinite values, and all values must be a 'normal' letter, digit, underscore, or dot")}
  all(is.behaved(x) & x!="" & !grepl("[^[:alnum:]_:\\.]",x))
}

vIsBehavedStrand <- function(x=NULL){
  if(argNotGiven(x)){return("Must not have NA, NAN, or Infinite values, and all entries must be '+' or '-'")}
  all(is.behaved(x) & x%in%c("+","-") )
}

#' @export
vContainsWhitespace <- function(x=NULL){
  if(argNotGiven(x)){return("Must not have NA, NAN, or Infinite values, and values must not have any whitespace (spaces, tab characters, newlines etc.)")}
  all(is.behaved(x) & grepl("[[:space:]]",x))
}



#' @export
makeValidator <- function( specName, requiredCols ){
  # One day export this data to a factory-factory?
  colSpecList <- list(
    seqId =           list( class="character" , valFun=vIsBehaved ), # Every 'valFun' function in this list must  return a text explanation if its `x` arg is NULL
    sSeqId=           list( class="character" , valFun=vIsBehaved ),
    qSeqId=           list( class="character" , valFun=vIsBehaved ),
    seq =             list( class="character" , valFun=vIsBehaved ),
    sSeq =            list( class="character" , valFun=vIsBehaved ),
    qSeq =            list( class="character" , valFun=vIsBehaved ),
    alnSeq =          list( class="character" , valFun=vIsBehavedSameLengthStrings ),
    sAlnSeq =         list( class="character" , valFun=vIsBehaved ),
    qAlnSeq =         list( class="character" , valFun=vIsBehaved ),
    bedFname =        list( class="character" , valFun=vIsBehaved ),
    fastaFname =      list( class="character" , valFun=vIsBehaved ),
    seqFname =        list( class="character" , valFun=vIsBehaved ),
    start =           list( class="numeric" ,   valFun=vIsBehavedPositiveIntlikeFloat ),
    end =             list( class="numeric" ,   valFun=vIsBehavedPositiveIntlikeFloat ),
    sStart =          list( class="numeric" ,   valFun=vIsBehavedPositiveIntlikeFloat ),
    sEnd =            list( class="numeric" ,   valFun=vIsBehavedPositiveIntlikeFloat ),
    qStart =          list( class="numeric" ,   valFun=vIsBehavedPositiveIntlikeFloat ),
    qEnd =            list( class="numeric" ,   valFun=vIsBehavedPositiveIntlikeFloat ),
    strand =          list( class="character" , valFun=vIsBehavedStrand ),
    name =            list( class="character" , valFun=vIsBehavedSimpleText ),
    fastaFname =      list( class="character" , valFun=vIsBehaved ),
    sFastaFname =     list( class="character" , valFun=vIsBehaved ),
    qFastaFname =     list( class="character" , valFun=vIsBehaved ),
    score =           list( class="numeric" ,   valFun=vIsBehaved ),
    minSize =         list( class="character" , valFun=vIsBehavedPositiveIntlikeFloat ),
    maxSize =         list( class="character" , valFun=vIsBehavedPositiveIntlikeFloat ),
    optimalSize =     list( class="character" , valFun=vIsBehavedPositiveIntlikeFloat ),
    productMinSize =  list( class="character" , valFun=vIsBehavedPositiveIntlikeFloat ),
    productMaxSize =  list( class="character" , valFun=vIsBehavedPositiveIntlikeFloat ),
    internal =        list( class="character" , valFun=vIsBehaved ),
    thickStart =      list( class="numeric" ,   valFun=vIsBehavedPositiveIntlikeFloat ),
    thickEnd =        list( class="numeric" ,   valFun=vIsBehavedPositiveIntlikeFloat ),
    itemRgb =         list( class="character" , valFun=vIsBehaved ),
    blockCount =      list( class="integer" ,   valFun=vIsBehaved ),
    blockSizes =      list( class="character" , valFun=vIsBehavedCommaSepInts ),
    blockStarts =     list( class="character" , valFun=vIsBehavedCommaSepInts )
  )

  objSpecList <- list(
    list( c("start","end") , function(vm){all(vm$start <= vm$end)} , "Values in the 'start' column must be smaller than or equal to values in 'end' column" ),
    list( c("qStart","qEnd") , function(vm){all(vm$qStart <= vm$qEnd)} , "Values in the 'qStart' column must be smaller than or equal to values in 'qEnd' column. Plus-to-minus mappings should be reflected by the subject coordinates having sEnd > sStart." ),
    list( c("sStart","sEnd","strand") , function(vm){ vm[,all(fifelse(sEnd<sStart,strand=="-",strand=="+"))] }  , "Where sEnd is less than sStart, this indicates the sequence is reversed in the coordinate system, and hence it should the value in strand should be \"-\"" )
  )

  if(any(!requiredCols %in% names(colSpecList))){stop("Cannot build validator: No information for requiredCols in colSpecList")}

  #colSpecs <- colSpecList[requiredCols]
  objSpecsli <- sapply(seq_len(length(objSpecList)),function(i){
    all(objSpecList[[i]][[1]] %in% requiredCols)
  })
  objSpecs <- objSpecList[objSpecsli]

  #dev library(data.table); validateMe = d.t(seqId="penis",seq=3.0,start=5,end=2); objName=deparse(substitute(validateMe)) ;requiredCols=c("seqId","start","end");croak=T
  validator <- function(validateMe=NULL,objName=deparse(substitute(validateMe)),croak=FALSE,message=NULL,showSpecs=FALSE){
    if(showSpecs==TRUE){
      o <- d.t(
        Required_Column = requiredCols
      )[
        ,idx:=1L:.N
      ][
        ,Class:=get(Required_Column,colSpecList)$class %>% force
        ,by=.(idx)][
          ,Notes:=get(Required_Column,colSpecList)$valFun() %>% force
          ,by=.(idx)][
            ,idx:=NULL
          ][]
      return(o)
    }
    if(argNotGiven(validateMe)){ stop("Argument for `validateMe` must be given unless `showSpecs==TRUE`") }
    if(any(!requiredCols %in% colnames(validateMe))){
      if(croak){
        drawConsoleLine()
        ce("VALIDATION FAILED")
        drawConsoleLine()
        ce("Object `",objName,"` does not meet the specifications for a ",specName,":\n\tThe following columns are required but not present:")
        ce(paste0(requiredCols[!requiredCols %in% colnames(validateMe)],collapse=", "),"\n")
        if(argGiven(message)){ ce("Additional messages: ", message ) }
        drawConsoleLine()
        stop("See validator output above")
      } else {
        return(FALSE)
      }
    }

    oc <- d.t(
      Column_Name = union(requiredCols,colnames(validateMe))
    )[,idx:=1L:.N][]

    oc[,Required_Name:=Column_Name %in% requiredCols]
    oc[,Reserved_Name:=Column_Name %in% names(colSpecList)]

    oc[Reserved_Name==TRUE,c("Required_Class","Class","Pass_Class","Pass_Column_Rules"):={
      # Column_Name = "name"
      rcl <- get(Column_Name,colSpecList)$class
      cl <- validateMe[,class(get(Column_Name))]
      pc <- (rcl==cl)
      pr <- get(Column_Name,colSpecList)$valFun( get(Column_Name,validateMe) )

      .(rcl,cl,pc,pr)
    },by=.(idx)]
    oc[Reserved_Name==FALSE,Pass_Class:=NA]
    oc[Reserved_Name==FALSE,Pass_Column_Rules:=NA]
    #oc[Reserved_Name==FALSE,Notes:="Column not subject to specifications."]
    if(any(oc[Reserved_Name==TRUE,Pass_Class==FALSE | Pass_Column_Rules==FALSE])){
      oc[Reserved_Name==TRUE & (Pass_Class==FALSE | Pass_Column_Rules==FALSE),Notes:=get(Column_Name,colSpecList)$valFun()[1],by=.(idx)]
    }

    oc[,idx:=NULL]

    if(any(oc[Reserved_Name==TRUE,Pass_Class==FALSE | Pass_Column_Rules==FALSE])){
      if(croak){
        drawConsoleLine()
        ce("VALIDATION FAILED")
        drawConsoleLine()
        ce("Object `",objName,"` does not meet the specifications for a ",specName,":\n\tCheck type is correct and see `Notes` column for things to check:")
        print(oc)
        drawConsoleLine()
        stop("See validator output above")
      } else {
        return(FALSE)
      }
    }
    wc <- ldtply(objSpecs,function(spec){
      #browser()
      d.t(
        Pass_Table_Rule = spec[[2]](validateMe),
        Rule = spec[[3]]
      )
    })

    if(any(wc$Pass_Table_Rule!=TRUE)){
      if(croak){
        drawConsoleLine()
        ce("VALIDATION FAILED")
        drawConsoleLine()
        ce("Object `",objName,"` does not meet the specifications for a ",specName,":\n\tSee `Rule` column for things to fix:")
        print(wc)
        drawConsoleLine()
        stop("See validator output above")
      } else {
        return(FALSE)
      }
    }

    return(TRUE)
  }
  return(validator)
}


#' @export
hasCol <- function(dt,colName){
  exists(colName,where=dt,inherits=FALSE)
}



















