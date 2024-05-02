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

  # Validate any given DTs
  if(argGiven(coordDT)){ is_coordDT(coordDT,croak=T,message = "Note: If the strand column is missing, you may still run an unstranded extraction using `stranded=FALSE`") }
  if(argGiven(seqDT))  { is_seqDT(  seqDT  ,croak=T) }

  # Check args make sense
  if(argGiven(bedFname) & argGiven(coordDT)){ stop("Coords provided as both file (`bedFname`) and DT (`coordDT`). Choose one please.") }
  if(sum(argGiven(c(fastaFname,seqDT,coordDT$fastaFname)))>1){ stop("Sequence provided more than one way (`fastaFname`) and DT (`seqDT`). Choose one please.") }
  if(sum(argGiven(c(fastaFname,seqDT,coordDT$fastaFname)))<1){ stop("Sequence (fasta) must be provided either as a file (with `fastaFname=`), or in a seqDT (with `seqDT=`), or as a column of the coordDT named \"fastaFname\" (with `coordDT`).") }

  if((argNotGiven(bedFname) & argNotGiven(coordDT))){ stop("Coords must be provided either as a BED file (with `bedFname=`) or a coordDT (with `coordDT=`).") }

  if(argGiven(bedFname) & argGiven(coordDT$fastaFname)){
    if(nu(coordDT$fastaFname)>1){
      warning("You are using a single BED file to provide coordinates for sequences in multiple files. That seems uncommon, but not inherently wrong. Let's try it anyway ...")
    }
  }

  # If the bed is in a file, read it in. It could be argued this is wasteful since we are going to write it out again. But this does let us check it and manipulate the name col, etc.
  if(argGiven(bedFname)){ coordDT <- readBed(bedFname) }

  # If a different name column is requested, check all is well
  if(argGiven(nameCol)){
    if(!has_name(coordDT)){
      warning(paste0("You requested the column `",nameCol,"` be used to provide the sequence names. The `name` column of the provided coordinate table will be overwritten."))
    }
    if(!nameCol %in% colnames(coordDT)){
      stop(paste0("The requested column `",nameCol,"` does not exist in the coordinate table, which has these columns:",colnames(coordDT)))
    }
  }
  if(argNotGiven(nameCol) & !has_name(coordDT)) { nameCol <- "seqId" }

  if(!has_name(coordDT)){ coordDT[,name:=get(nameCol)] }

  # If strandlessness requested, check it is ok
  if(stranded==FALSE){
    strandedArg <- ""
    if(has_strand(coordDT)){
      warning("You requested strandless extration. The `strand` column of the provided coordDT will be overwritten (with '+').")
    }
    coordDT[,strand:="+"]
  } else { # stranded requested
    has_strand(coordDT,croak=T)
  }

  killTf <- FALSE
  if(argGiven(seqDT)){ # All in DTs --- can be done without bedtools
    #return( extractSeqFastaDTBedDT(coordDT,seqDT,stranded=stranded) )
    warning("NOT PROPERLY IMPLEMENTED, EMAIL TIM! (Will go ahead and use a stupid method instead ... )")
    # Set up fasta if needed, and load its filename into coordDT
    killTf <- TRUE
    writeFasta(seqDT,tfFname<-tempfile(fileext=".fasta"))
    coordDT[,fastaFname:=tfFname]
  }

  ## Anything from here we are committed to calling bedtools getfasta

  killFo <- FALSE
  if(!argGiven(outFastaFname)){
    outFastaFname <- tempfile(fileext = ".fasta")
    killFo <- TRUE
  } else {
    unlink(outFastaFname)
  }
  # Over all target sequences.
  coordDT[,{
    tbFname <- tempfile(fileext = ".bed")
    writeBed(.SD,tbFname)
    command <- paste0(bedToolsBin," getfasta -fi ",fastaFname," -bed ",tbFname," -s >> ",outFastaFname)
    system(command)
    unlink(tbFname)
  },by=.(fastaFname)]

  if(killTf){unlink(tfFname,paste0(tfFname,".fai"))}
  if(killFo){t <- readFasta(outFastaFname); unlink(outFastaFname); return(cbind(coordDT,t[,.(seq)]))} else {ce("Sequences saved as ",outFastaFname)}

}
