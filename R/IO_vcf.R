
#' @export
readVcfRaw <- function(vcfFname){
  indIds <- fread(cmd=paste0("cat ",vcfFname," | grep '^#CHROM' | head -1"),header=F) %>% as.character() %>% `[`(10:length(.))
  fread(cmd=paste0("cat ",vcfFname,' | grep -v "^#"'))[] %>% setnames(colnames(.),c("seqId","pos","varId","refAllele","altAlleles","varQuality","filter","varInfo","fmt",indIds))
}

#' @export
vcfVarList <- function(vcfDT){
  v <- copy(vcfDT)
  setkey(v,seqId,pos)
  v[,.(alleleLiteral=c(refAllele,strsplit(altAlleles,",")[[1]])),by=.(seqId,pos)][,alleleNumeric:=1:.N,,by=.(seqId,pos)][]
}

#' @export
parseGenotypes <- function(ref,alts,gts,colNamePrefix="Allele_"){
  # colNamePrefix="Allele_"
  # ref <- c("a","b","c")
  # alts <- c("A,aa","B,bb,BB","C,Cc")
  # gts <- c("1/0","2/4","1/.")
  # \ dev
  require(stringr)
  if(!same(c(length(ref),length(alts),length(gts)))){ stop("'ref', 'alts', and 'gts' must have equal lengths") }
  n <- length(ref)
  # establish separator
  sep <- substr(gts[1],2,2)
  ploidy <- stringr::str_count(gts[1],sep)+1
  out <- matrix(character(),nrow=n,ncol=ploidy,dimnames = list(NULL,paste0(colNamePrefix,"_",1:ploidy)))
  data.table(
    alleleList=lapply(1:n,function(i){
      c(ref[i],strsplit(alts[i],",")[[1]])
    }),
    gt=gts
  )[,idx:=1:.N][][,{
    #browser()
    # print(strsplit(gt,sep)[[1]])
    # print(alleleList[[1]][strsplit(gt,sep)[[1]] %>% as.integer])
    out[idx,] <<- alleleList[[1]][strsplit(gt,sep)[[1]] %>% as.integer %>% `+`(1)]
  },by=.(idx)] %>% suppressWarnings()
  as.data.table(out)
}

#' @export
vcf2varDT <- function(vcfFname,keepInfoFields=FALSE,keepGtFields=FALSE,maxAllelesPerSite=2,parseGTs=TRUE){ #2 reminds us it is biallelics only
  # maxAllelesPerSite=2
  # vcfFname="data/calls/Athal_muc_calls.vcf"
  # /DEV
  require(stringr)
  v <- readVcfRaw(vcfFname)
  v <- v[ stringr::str_count(altAlleles,",") < maxAllelesPerSite ]
  if(!same(v$fmt)){ stop("I do not currently work for vcf files with differently-formatted genotype fields across rows. Currently the fields ") }
  fieldNames <- strsplit(v$fmt[1],":")[[1]]
  nFields <- length(fieldNames)
  indIds <- colnames(v) %>% `[`(10:length(.))
  nInds <- length(indIds)
  gtFieldMat <- matrix(character(),ncol=nFields*nInds,nrow=nrow(v))
  colnames(gtFieldMat) <- expand.grid(fieldNames,"_",indIds) %>% apply(1,function(x) paste0(x,collapse = "") )
  for(i_I in 1:nInds){
    #i_I=2
    s <- strsplit(v[,get(indIds[i_I])],":") %>% unlist
    for(i_gf in 1:nFields){
      #i_gf=2
      j_m <- ((i_I-1)*nFields) + i_gf
      gtFieldMat[,j_m] <- s[( 1:length(s) %% nFields == (i_gf%%nFields) )]
    }
  }
  v <- cbind(v,as.data.table(gtFieldMat))
  if( parseGTs==TRUE ){
    if(!"GT" %in% fieldNames){ stop("'GT' is not among the field names in the VCF, i.e., it isn't clear there is any genotype info in here. If you are using some non-standard VCF, then you can parse columns containing VCFesque genotype calls (e.g. '1/1', '0|1', '1:2:2', ...) and then parse them manually with `parseVcfGenotypes()`") }
    for(iId in paste0("GT_",indIds)){
      v %<>% cbind(v[,parseGenotypes(refAllele,altAlleles,get(iId,v),paste0(iId,"_Allele"))]) #utterly horrid. Sorry. Will improve some day.
    }
  }
  v
}
