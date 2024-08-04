
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
parseGenotypes <- function(ref,alts,gts,colNamePrefix="Allele_",gtEncoding=c("literal","altCount")){
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

  if(gtEncoding[1]=="altCount"){
    if(any(grepl("[2-9]",gts))){ stop("When 'gtEncoding'==\"altCount\", a maximum of two alleles per site may be present. Please filter for biallelic sites. If this is being called from vcf2varDT(), then set 'maxAllelesPerSite' to 2") }
    out <- list(strsplit(gts,sep) %>% sapply(function(x)sum(as.integer(x))) %>% suppressWarnings())
    names(out) <- colNamePrefix
    as.data.table(out)
  } else if(gtEncoding[1]=="literal"){
    out <- matrix(character(),nrow=n,ncol=ploidy,dimnames = list(NULL,paste0(colNamePrefix,1:ploidy)))
    data.table(
      alleleList=lapply(1:n,function(i){
        c(ref[i],strsplit(alts[i],",")[[1]])
      }),
      gt=gts
    )[,idx:=1:.N][][,{
      # browser()
      # print(strsplit(gt,sep)[[1]])
      # print(alleleList[[1]][strsplit(gt,sep)[[1]] %>% as.integer])
      out[idx,] <<- alleleList[[1]][strsplit(gt,sep)[[1]] %>% as.integer %>% `+`(1)]
    },by=.(idx)] %>% suppressWarnings()
    return(as.data.table(out))
  } else {
    stop(paste0("'gtEncoding' argument value ",gtEncoding[1]," not recognised. Must be one of \"literal\" or \"altCount\""))
  }
}

#' @export
vcf2varDT <- function(vcfFname,keepInfoFields=FALSE,keepGtFields=FALSE,maxAllelesPerSite=2,parseGTs=c("literal","altCount","no"),IDprefix="ID:"){ #2 reminds us it is biallelics only
  # DEV
  # maxAllelesPerSite=2
  # IDprefix="ID:"
  # vcfFname="data/calls/Athal_muc_varEff.vcf"
  # # vcfFname = "/data/gpfs/projects/punim1869/users/amadhusudans/workspace/covid_data/testing_index_merging/testMerge.vcf"
  # /DEV
  require(stringr)
  v <- readVcfRaw(vcfFname)
  v <- v[ stringr::str_count(altAlleles,",") < (maxAllelesPerSite-1) ]
  indIds <- colnames(v) %>% `[`(10:length(.)) %>% paste0(IDprefix,.)
  colnames(v)[10:ncol(v)] <- indIds
  if(!same(v$fmt)){ stop(paste0("I do not currently work for vcf files with differently-formatted genotype fields across rows. Currently the fields include these: ",unique(v$fmt))) }
  fieldNames <- strsplit(v$fmt[1],":")[[1]]
  nFields <- length(fieldNames)
  nInds <- length(indIds)
  gtFieldMat <- matrix(character(),ncol=nFields*nInds,nrow=nrow(v))
  colnames(gtFieldMat) <- expand.grid(indIds,"_",fieldNames) %>% apply(1,function(x) paste0(x,collapse = "") )
  #browser()
  for(i_I in 1:nInds){
    #i_I=2
    s <- strsplit(v[,get(indIds[i_I])],":") %>% unlist
    for(i_gf in 1:nFields){
      #i_gf=1
      j_m <- ((i_gf-1)*nInds) + i_I
      # ce("ind    ",i_I)
      # ce("field  ",i_gf)
      # ce("outCol ",j_m)
      # ce("named  ",colnames(gtFieldMat)[j_m])
      # ce("entries: ",s[( 1:length(s) %% nFields == (i_gf%%nFields) )][1:10] %>% pastec)
      gtFieldMat[,j_m] <- s[( 1:length(s) %% nFields == (i_gf%%nFields) )]
    }
  }
  v <- cbind(v,as.data.table(gtFieldMat))
  if( parseGTs[1]!="no" ){
    if( parseGTs[1] %!in% c("literal","altCount") ){ stop(paste0("'parseGTs' argument value ",parseGTs[1]," not recognised. Must be one of \"literal\", \"altCount\", or \"no\"")) }
    if("GT" %!in% fieldNames){ stop("'GT' is not among the field names in the VCF, i.e., it isn't clear there is any genotype info in here. If you are using some non-standard VCF, then you can parse columns containing VCFesque genotype calls (e.g. '1/1', '0|1', '1:2:2', ...) and then parse them manually with `parseVcfGenotypes()`") }
    for(iId in paste0(indIds,"_GT")){
      #iId <- paste0(indIds,"_GT")[3]
      v %<>% cbind(v[,parseGenotypes(ref=refAllele,alts=altAlleles,gts=get(iId,v),colNamePrefix=paste0(iId,"_alleleCount"),gtEncoding=parseGTs)]) #utterly horrid. Sorry. Will improve some day.
    }
  }
  v
}
