#' Extract the subsequences described by a coordinate table (e.g. read in from a bed file or similar), from a sequence list (e.g. fasta file or dt).
#'
#' @param paramName What it is. Defaults?
#' @returns What the function returns.
#' @description
#' Describe what the function does in detail.
#' @examples
#' Example 1 code
#'
#' Example 2 code
#' @export
extractSeq <- function (
    bedFname=NULL,
    coordDT=NULL,
    fastaFname=NULL,
    seqDT=NULL,
    outFastaFname=NULL,
    stranded=TRUE,
    nameCol=NULL,
    bedToolsBin=system("which bedtools",intern=T)
){
  coordDT <- copy(coordDT)
  if(!is.null(coordDT) & !is.null(nameCol)){
    coordDT[,name:=get(nameCol)]
  }

  strandedArg <- "-s"

  if(!is.null(bedFname) & !is.null(coordDT)){ stop("Coords provided as both file (`bedFname`) and DT (`coordDT`). Choose one please.") }
  if(!is.null(fastaFname) & !is.null(seqDT)){ stop("Sequence provided as both a file (`fastaFname`) and DT (`seqDT`). Choose one please.") }

  if((is.null(bedFname) & is.null(coordDT))){ stop("Coords must be provided either as a file or a DT.") }
  if((is.null(fastaFname) & is.null(seqDT))){ stop("Sequence (fasta) must be provided either as a file or a DT.") }


  if(!is.null(coordDT) & stranded==FALSE){
    coordDT[,strand:="+"]
  } else {
    has_strand(coordDT)
  }
  if(!is.null(bedFname) & stranded==FALSE){
    strandedArg <- ""
  }

  if(!is.null(coordDT) & !is.null(seqDT)){ # All in DTs
    is_coordDT(coordDT)
    is_seqDT(seqDT)
    #return( extractSeqFastaDTBedDT(coordDT,seqDT,stranded=stranded) )
    warning("NOT YET IMPLEMENTED, EMAIL TIM! (Will go ahead and use a stupid method instead)")
  }

  ## Anything from here we are committed to calling bedtools getfasta

  killTfb <- FALSE
  if(is.null(bedFname)){ bedFname <- tempfile(fileext = ".bed"); writeBed(coordDT,bedFname); killTfb <- TRUE }

  killTff <- FALSE
  if(is.null(fastaFname)){ fastaFname <- tempfile(fileext=".fasta") ; writeFasta(seqDT,fastaFname); killTff <- TRUE }

  killTfo <- FALSE
  if( is.null(outFastaFname) ){ outFastaFname <- tempfile(fileext=".fasta"); internArg=TRUE; killTfo <- TRUE }

  nameBedArg <- if(ncol(fread(bedFname,nrows=1,header=F))>=4){"-name"}else{""}
  # outFastaFname <- tfo

  command <- paste0(bedToolsBin," getfasta -fi ",fastaFname," -bed ",bedFname," ",strandedArg," ",nameBedArg," > ",outFastaFname)
  system(command)

  #debugonce(readFasta)
  out <- readFasta(outFastaFname)

  # cleanup
  if(killTfb){unlink(bedFname)}
  if(killTff){unlink(fastaFname)}
  if(killTfo){unlink(outFastaFname)}

  out
}
