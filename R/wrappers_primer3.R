#' @export
primer3 <- function(
    seqDT,
    primer3path=system("which primer3",intern=T),
    awkPath=system("which awk",intern=T)
){
  # Expect columns: seqId, seq, (optional)internal(bool), (optional)optimalSize(int), (optional)minSize(int), (optional)maxSize, (optional)productMinSize, (optional)productMaxSize
  #dev seqDT <- readFasta("/data/gpfs/projects/punim1869/shared_data/misc_sequence/btr_btrLike_queries_Morex_GP.fasta")
  #dev primer3path="/home/mrabanuswall/bin/primer3/src/primer3_core"
  #dev awkPath=system("which awk",intern=T)
  seqDT <- copy(seqDT)
  seqDT[,idxBioDT:=1:.N]
  seqDT[,seq:=gsub("-","",seq)]


  #browser()
  itl        <- hasInternalFlag(seqDT,error=F)
  opSize     <- hasOptimalSize(seqDT,error=F)
  mmSize     <- hasMinMaxSize(seqDT,error=F)
  mmProdSize <- hasMinMaxProductSize(seqDT,error=F)


  seqDT[,{

    tf <- tempfile(fileext=".primer3args")

    write(paste0("SEQUENCE_ID=",seqId),tf,append=TRUE)
    seqIn <- if(grepl("[^AGCTNagctn]+",seq)){
      warning(paste0("Ambiguous bases detected in ",seqId," --- replacing with `N`"))
      gsub("[^AGCTNagctn]+","N",seq)
    } else {
      seq
    }
    write(paste0("SEQUENCE_TEMPLATE=",seqIn),tf,append=TRUE)
    write("PRIMER_TASK=generic",tf,append=TRUE)


    #browser()

    if(itl){
      write("PRIMER_PICK_INTERNAL_OLIGO=1",tf,append=TRUE)
    }
    if(opSize){
      write(paste0("PRIMER_OPT_SIZE=",optimalSize),tf,append=TRUE)
    }
    if(mmSize){
      write(paste0("PRIMER_MIN_SIZE=",minSize),tf,append=TRUE)
      write(paste0("PRIMER_MAX_SIZE=",maxSize),tf,append=TRUE)
    }
    if(mmProdSize){
      write(paste0("PRIMER_PRODUCT_SIZE_RANGE=",minProductSize,"-",maxProductSize),tf,append=TRUE)
    }
    write("=",tf,append=TRUE)

    ce("Input file: ")
    system(paste("cat ",tf))

    command <- paste0(primer3path," ",tf," | ",awkPath, " -F'=' " , aScriptPrimer3generic())

    # system(command)
    o <- fread(cmd=command,colClasses=c("character","character","numeric","character","numeric","integer","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","integer","numeric"))
    # print(o)
    #browser()

    unlink(tf)
    o[,seqId:=seqId][]

  },by=.(idxBioDT)][,idxBioDT:=NULL][]
}
