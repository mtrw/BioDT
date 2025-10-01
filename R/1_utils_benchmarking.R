#' @export
benchTim <- function(fun,...,n=10000){
  time <- system.time(
    for(i in 1:n){ fun(...) }
  )[["user.self"]]
  names(time) <- "User Time (s)"
  time
}

#' @export
benchTimCmp <- function(funList,...,n=10000,plot=TRUE){
  fll <- as.list(substitute(funList))
  fNames <- fll[2:length(fll)] %>% sapply(function(e){ deparse(e) })
  userTimes <- sapply( funList , function(f){ benchTim(f,...) })
  if(plot==TRUE){
    null_plot(seq_along(fNames) %plusMinus% 1 ,userTimes,xlab="Function",ylab="User Time",xaxt="n")
    title(main=paste("Benchmarking report for",n,"runs"))
    axis(1, at=seq_along(fNames) , labels=fNames )
    points(seq_along(fNames),userTimes,pch=18,cex=2)
    lines(seq_along(fNames),userTimes,type="b")
  }
  names(userTimes) <- fNames
  return(userTimes)
}

# benchTim(sum,1:1000)
# benchTim(mean,1:1000)
# benchTim(prod,1:1000)
# benchTim(exp,1:1000)
#
# benchTimCmp(list(sum,mean,prod,exp),1:1000)

#' @export
startClock <- function(){
  while(exists(varName<-paste0(".",sample(c(LETTERS,letters)),collapse=""),where=globalenv())){}
  assign(varName,Sys.time(),envir = globalenv())
  varName
}

#' @export
stopClock <- function(clockID){
  timeElapsed <- print(Sys.time() - get(clockID,envir = globalenv()))
  rm(clockID)
  timeElapsed
}
