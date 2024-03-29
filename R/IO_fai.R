#' Read .fai files into a data.table
#'
#' @param faiFname A character vector. The .fai filename(s). [no default]
#' @returns A data.table with .fai flavour.
#' @export
readFai <- function( #read a fai file
  faiFname  #fastaFname="/data/gpfs/projects/punim1869/shared_data/misc_sequence/btr_btrLike_queries_Morex_GP.fasta"
){
  out <- ldply(faiFname,function(fn){ #dev fn <- "/data/gpfs/projects/punim1869/shared_data/misc_sequence/btr_btrLike_queries_Morex_GP.fasta.fai"
    fread(fn,select=1:2,header=F,col.names=c("seqName","length"))
  }) %>% setDT
  return(out)
}


#' Generate .fai files for given .fasta files.
#'
#' @param fastaFname A character vector. The input .fasta filename(s). [no default]
#' @param faiFname A character vector. The output .fai filename(s). [paste0(fastaFname,".fai")]
#' @param samtoolsBin A character vector. The .fai filename. [system("which samtools", intern=TRUE)]
#' @returns A data.table with .fai flavour.
#' @description
#' Wrapper for samtools faidx
#' @export
makeFai <- function( #make fai files
  fastaFname,  #fastaFname=c("/data/gpfs/projects/punim1869/shared_data/misc_sequence/btr_btrLike_queries_Morex_GP.fasta","/data/gpfs/projects/punim1869/shared_data/misc_sequence/btr_btrLike_queries_Morex_GP.fasta")
  faiFname=paste0(fastaFname,".fai"),
  samtoolsBin=system("which samtools", intern=TRUE)
){
  stopifnot("FastaFname and faiFname must be the same length" = length(fastaFname)==length(faiFname))
  l_ply(1:length(fastaFname),function(i){
    command <- paste0(samtoolsBin," faidx -o ",faiFname[i]," ",fastaFname[i])
    system(command)
  })
  return(faiFname)
}


#' Create a .fai-flavoured data.table directly from a .fasta file
#'
#' @param fastaFname A character vector. The .fasta filename(s). [no default]
#' @param samtoolsBin Character value. In case namesOnly==TRUE, where to find samtools. [system("which samtools", intern=TRUE)]
#' @returns A data.table with .fai flavour.
#' @description
#' Calls samtools faidx
#' @export
getFai <- function( #from a FASTA file straight into an R fai
    fastaFname,  #fastaFname="/data/gpfs/projects/punim1869/shared_data/misc_sequence/btr_btrLike_queries_Morex_GP.fasta"
    samtoolsBin=system("which samtools", intern=TRUE)
){
  out <- ldply(fastaFname,function(fn){
    command <- paste0(samtoolsBin," faidx -o /dev/stdout ",fn)
    fread(cmd=command,select=1:2,header=F,col.names=c("seqName","length"))
  }) %>% setDT
  return(out)
}
