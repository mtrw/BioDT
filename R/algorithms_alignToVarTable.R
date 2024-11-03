#' Convert homology search/alignment output to a variant table DT.
#'
#' @param aln What it is. Defaults?
#' @returns What the function returns.
#' @description
#' Describe what the function does in detail.
#' @examples
#' Example 1 code
#'
#' Example 2 code

#' #' @export
# CALL_varTableFromAlns <- function(x){ #call internal direct later
#   INTERNAL_varTableFromAlns(x)
# }

#' @export
varTableFromAlns <- function(aln){

  #aln <- copy(bl)
  aln <- copy(aln)

  aln[,qInc:=fifelse(qStart<qEnd,1L,-1L)]
  aln[,sInc:=fifelse(sStart<sEnd,1L,-1L)]
  aln[,qAlnSeq:=fifelse(sInc==1,qAlnSeq,rc(qAlnSeq))]
  aln[,sAlnSeq:=fifelse(sInc==1,sAlnSeq,rc(sAlnSeq))]
  aln[,c("sStart_bl","sEnd_bl"):=.(sStart,sEnd)] # save
  aln[sInc==-1L,c("sStart","sEnd"):=.(sEnd,sStart)]
  aln[sInc==-1L,sInc:=1L]


  #dev
  {
  #   aln <- copy(bl[c(28,33)])
  #   #aln <- copy(bl)[c(6,18)]
  #   aln <- copy(bl)
  #
  #   aln[,qInc:=fifelse(qStart<qEnd,1L,-1L)]
  #   aln[,sInc:=fifelse(sStart<sEnd,1L,-1L)]
  #   aln[,qAlnSeq:=fifelse(sInc==1,qAlnSeq,rc(qAlnSeq))]
  #   aln[,sAlnSeq:=fifelse(sInc==1,sAlnSeq,rc(sAlnSeq))]
  #   aln[,c("sStart_bl","sEnd_bl"):=.(sStart,sEnd)] # save
  #   aln[sInc==-1L,c("sStart","sEnd"):=.(sEnd,sStart)]
  #   aln[sInc==-1L,sInc:=1L]
  #
  #   # cutl <- 1; cutr <- 20
  #   # aln[,sAlnSeq:=substr(sAlnSeq,cutl,cutr)]
  #   # aln[,qAlnSeq:=substr(qAlnSeq,cutl,cutr)]
  #   aln[,idx:=1:.N]
  #   aln[,ss:=sAlnSeq]
  #   aln[,qs:=qAlnSeq]
  #   aln[,si:=sInc]
  #   aln[,qi:=qInc]
  #   aln[,sst:=sStart]
  #   aln[,qst:=qStart]
  #   aln[,se:=sEnd]
  #   aln[,qe:=qEnd]
  #
  #   c <- aln[,CALL_varTableFromAlns(.SD),by=.(qSeqId,sSeqId,idx,ss,qs,si,qi,sst,se,qst,qe)][!grepl(":",qState)]
  #   c[,.N,by=.(sPos,sSeqId,sState)][,.N,by=.(sPos,sSeqId)][N>1]
  #   c[,.N,by=.(sPos,sSeqId,sState)][,.N,by=.(sPos,sSeqId)][N>1][,N:=NULL] -> badPos
  #   (cb <- c[badPos,on=colnames(badPos)][order(sPos)])
  #
  #   # print(cb[1:2])
  #   # printAln(cb[1,c(qs,ss)])
  #   # printAln(cb[2,c(qs,ss)])
  #
  #   # print(cb[3:4])
  #   # printAln(cb[3,c(qs,ss)])
  #   # printAln(cb[4,c(qs,ss)])
  #
  }


  #c <- aln[,CALL_varTableFromAlns(.SD),by=.(qSeqId,sSeqId)]
  c <- aln[,INTERNAL_varTableFromAlns(.SD),by=.(qSeqId,sSeqId)]

  #separate indels
  c[,indel:=FALSE]
  c[grep("i:|d:",sState),indel:=TRUE]
  id <- c[indel==TRUE,][,indel:=NULL][]
  c <- c[indel==FALSE,][,indel:=NULL][]

  #process indels
  id[,qState:=fifelse(grepl("i:",sState),-as.integer(sub(".*:","",sState)),as.integer(sub(".*:","",sState)))][,sState:=NULL][,qPos:=NULL]



  #uniqualise multiple calls for same subject site (caused by overlapping alignments)
  # Where alignments overlap on subject, multiple calls may occur at the same site. If they are the same alt allele, keep the call. If not, delete the site.
  c[,idx:=1:.N]
  #c: kill ambigs
  c[,.N,by=.(qSeqId,sSeqId,sPos)][N>1,][,N:=NULL][] -> flagAmbigs #multiple alignments, possibility for ambiguous calls
  c[flagAmbigs,on=.(qSeqId,sSeqId,sPos)][,.(nS=nu(qState)),by=.(qSeqId,sSeqId,sPos)][nS>1,.(sSeqId,sPos)] -> killAmbigPos
  c <- c[!killAmbigPos,on=.(sSeqId,sPos)]
  #c: trim remaining. These are the same, so just trim the excess
  c[,.N,by=.(qSeqId,sSeqId,sPos)][N>1,][,N:=NULL][] -> restAmbigs
  c[restAmbigs,on=.(qSeqId,sSeqId,sPos)][,.SD[2:.N,],by=.(qSeqId,sSeqId,sPos)][,.(qSeqId,sSeqId,sPos)] -> killIdx
  c <- c[!killIdx,on=.(qSeqId,sSeqId,sPos)]
  # --------------------- #




  id[,idx:=1:.N]
  #id: kill ambigs
  id[,.N,by=.(qSeqId,sSeqId,sPos)][N>1,][,N:=NULL][] -> flagAmbigs #multiple alignments, possibility for ambiguous calls
  id[flagAmbigs,on=.(qSeqId,sSeqId,sPos)][,.(nS=nu(qState)),by=.(qSeqId,sSeqId,sPos)][nS>1,.(qSeqId,sSeqId,sPos)] -> killAmbigPos
  id <- id[!killAmbigPos,on=.(qSeqId,sSeqId,sPos)]
  #id: trim remaining. These are the same, so just trim the excess
  id[,.N,by=.(qSeqId,sSeqId,sPos)][N>1,][,N:=NULL][] -> restAmbigs
  id[restAmbigs,on=.(qSeqId,sSeqId,sPos)][,.SD[2:.N,],by=.(qSeqId,sSeqId,sPos)][,.(qSeqId,sSeqId,sPos)] -> killIdx
  id <- id[!killIdx,on=.(qSeqId,sSeqId,sPos)]

  #sanity
  if(sum(
    c[, .N,by=.(qSeqId,sSeqId,sPos)][N>1] %>% nrow, #empty or bad
    id[,.N,by=.(qSeqId,sSeqId,sPos)][N>1] %>% nrow, #empty or bad
    c[,.N,by=.(sPos,sSeqId,sState)][,.N,by=.(sPos,sSeqId)][N>1] %>% nrow #empty or bad
  )>0){ stop("Bug detected. Fetch the rubber duck (and email Tim)") }


  #fill c
  tmp <- c[,.N,by=.(sPos,sSeqId,sState)][,N:=NULL][,idx:=1:.N][]
  cfiller <- tmp[,.(qSeqId=unique(c$qSeqId)),by=.(sPos,sSeqId,sState,idx)][,idx:=NULL][]
  c <- c[,sState:=NULL][cfiller,on=.(sPos,sSeqId,qSeqId)]
  c[,idx:=NULL]
  rm(tmp,cfiller)

  #fill id
  tmp <- id[,.N,by=.(sPos,sSeqId)][,N:=NULL][,idx:=1:.N][]
  idfiller <- tmp[,.(qSeqId=unique(id$qSeqId)),by=.(sPos,sSeqId,idx)][,idx:=NULL][]
  id <- id[idfiller,on=.(sPos,sSeqId,qSeqId)]
  id[,idx:=NULL]
  rm(tmp,idfiller)

  # detect which NA calls are due to no alignment coverage vs covered but no difference
  aln_ol <- aln[,.(sSeqId,qSeqId,sStart=pmin(sStart,sEnd),sEnd=pmax(sStart,sEnd))]
  setkey(aln_ol,sSeqId,qSeqId,sStart,sEnd) #alignments, filtered)


  #c
  c[,sStart  := sPos]
  c[,sEnd    := sPos]
  setkey(c, sSeqId,qSeqId,sStart,sEnd)
  ol_c <- foverlaps(c,aln_ol,nomatch=NULL,mult="first")[,.(sSeqId,qSeqId,sPos)]
  c[,coveredAlignment:=FALSE]
  c[ol_c,on=.(sSeqId,qSeqId,sPos),coveredAlignment:=TRUE]
  c[coveredAlignment==TRUE & is.na(qState),qState:=sState]
  rm(ol_c)

  id[,sStart := sPos]
  id[,sEnd   := sPos+1]
  setkey(id,sSeqId,qSeqId,sStart,sEnd)
  ol_id <- foverlaps(id,aln_ol,nomatch=NULL,mult="first")[,.(sSeqId,qSeqId,sPos)]
  id[,coveredAlignment:=FALSE]
  id[ol_id,on=.(sSeqId,qSeqId,sPos),coveredAlignment:=TRUE]
  id[coveredAlignment==TRUE & is.na(qState),qState:=0L]
  id[coveredAlignment==TRUE]
  rm(ol_id)
  rm(aln_ol)

  return(list(snps=c,indels=id))

}

