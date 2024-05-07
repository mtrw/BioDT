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
    stranded=FALSE,
    nameCol=NULL,
    bedToolsBin=system("which bedtools",intern=T)
){
  # DEV
  # bedFname=NULL
  # fastaFname="/data/gpfs/projects/punim1869/shared_data/misc_sequence/btr_btrLike_queries_Morex_GP.fasta"
  # coordDT=d.t(seqId=readFasta(fastaFname)$seqId)
  # seqDT=NULL
  # outFastaFname=NULL
  # stranded=TRUE
  # nameCol=NULL
  # bedToolsBin=system("which bedtools",intern=T)
  # \ DEV

  coordDT <- copy(coordDT)

  # Validate any given DTs
  if(argGiven(coordDT)){ is_seqListDT(coordDT,croak=T)}
  if(argGiven(seqDT))  { is_seqDT(  seqDT  ,croak=T) }

  # Check args make sense

  if(sum(argGiven(coordDT),argGiven(bedFname))!=1){ stop("Sequence names must be provided either in a coordinate table (`coordDT=`; 'seqId' column) or in a BED file (`bedFname=`; 'chr' column), and only one.") }
  if(sum(argGiven(coordDT),argGiven(bedFname))>1){                             stop("Coordinates can be provided as either a BED file (`bedFname=`) or a coordinate table (`coordDT=`), and only one.") }
  if(sum(argGiven(fastaFname),argGiven(seqDT),argGiven(coordDT$fastaFname))!=1){ stop("Sequences must be provided by either fasta files listed in the coordinate table (`coordDT=`; 'seq' column) or in a sequence table (`seqDT=`; 'seq' column) or in a fasta file (`fastaFname=`), and only one.") }
  #
  if(argGiven(bedFname) & argGiven(coordDT$fastaFname)){
    if(nu(coordDT$fastaFname)>1){
      warning("You are using a single BED file to provide coordinates for sequences in multiple files. That seems uncommon, but not inherently wrong. Let's try it anyway ...")
    }
  }

  # When coordDT isn't there, get it.
  # If the bed is in a file, read it in. It could be argued this is wasteful since we are going to write it out again. But this does let us check it and manipulate the name col, etc.
  if(argGiven(bedFname)){
    coordDT <- readBed(bedFname)
  }
  # Check for a possible error case
  if (!is.null(coordDT$start) | !is.null(coordDT$end)){
    is_coordDT(coordDT,croak=T,message="Presence of 'start'/'end' columns suggest you want to provide a coordinate table and extract sequence from particular coordinates, but validation has failed. Either correct these columns (see error messages) or submit a seqListDT (with a 'seqId' column but no 'start'/'end' columns) to extract whole sequences.")
  }

  # FUTURE: Here is probably where to check for the special case of all info given in DTs

  # Fill fasta file col, writing fasta if necessary

  killTf <- FALSE
  if(argGiven(fastaFname)){
    if(length(fastaFname)>1){ stop("When providing a coordinate table (`coordDT=`) without a 'fastaFname' column, only one fasta file (`fastaFname=`) can be given (otherwise, it could be ambiguous which 'seqId's to extract from which fasta file). To extract from multiple fastas, include them in a 'fastaFnames' column of the sequence list DT.") }
    coordDT[,fastaFname:=fastaFname]
  } else if (argGiven(seqDT)){
    killTf <- TRUE
    writeFasta(seqDT,tfFname<-tempfile(fileext=".fasta"))
    coordDT[,fastaFname:=tfFname]
  }

  # Fill start / end columns
  if(!has_start_end(coordDT)){
    #debugonce(getFai)
    fai <- coordDT[,{
      getFai(fastaFname)
    },by=.(fastaFname)]

    coordDT <- fai[coordDT,on=.(seqId)][,start:=1.0][,end:=as.numeric(seqLength)][,seqLength:=NULL]
  }

  # If a different name column is requested, check all is well
  if(argGiven(nameCol)){
    if(!has_name(coordDT)){
      warning(paste0("You requested the column `",nameCol,"` be used to provide the sequence names. The `name` column of the provided coordinate table will be overwritten."))
    }
    if(!nameCol %in% colnames(coordDT)){
      stop(paste0("The requested column `",nameCol,"` does not exist in the coordinate table, which has these columns:",colnames(coordDT)))
    }
    coordDT[,name:=get(nameCol)]
  } else if(!has_name(coordDT)){
    ce("Note: Column `name` not present in coordDT. Extracted sequences will be named using `seqId` column.")
    coordDT[,name:=seqId]
  }

  # If strandlessness requested, check it is ok
  if(stranded==FALSE){
    strandedArg <- ""
    if(has_strand(coordDT)){
      warning("You requested strandless extration. The `strand` column of the provided coordDT will be overwritten (with '+').")
    }
    coordDT[,strand:="+"]
  } else { # stranded requested
    has_strand(coordDT,croak=T,message = "Stranded extraction was requested (`stranded=TRUE`; default), but no strand column was given. Try re-running with")
  }

  #
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
  if(killTf){unlink(tfFname,paste0(tfFname,".fai")); coordDT[,fastaFname:=NULL][]}
  if(killFo){t <- readFasta(outFastaFname); unlink(outFastaFname); return(coordDT[,seq:=t$seq][])} else {ce("Sequences saved as ",outFastaFname)}
}


