
#' @export
lastz <- function(
    subjectFname=NULL,
    queryFname=NULL,
    subjectSeqDT=NULL,
    querySeqDT=NULL,
    outFastaFname=NULL,
    addAlignments=FALSE,
    addPctIdGaps=FALSE,
    addNgaps=FALSE,
    addColList=NULL,
    outFmtArg=NULL,
    outputColNames=NULL,
    outputColClasses=NULL,
    otherLastzArgs="",
    showCmd=FALSE,
    lastzBinary=system("which lastz",intern=T),
    samtoolsBinary=system("which samtools",intern=T),
    forceBioDToutput=F
){

  makeBioDTOutput <- if(all(argNotGiven(c(outFmtArg,outputColNames,outputColClasses)))){
    outFmtArg <- "general:name1,name2,ncolumn,start1,end1,start2,end2,strand1,strand2,blastid%,score"
    outputColNames <-   c("sSeqId",    "qSeqId" ,    "matchLength" , "sStart" , "sEnd", "qStart"  ,    "qEnd" , "sStrand" , "qStrand" ,   "pctId_noGaps" , "score")
    outputColClasses <- c("character", "character",  "integer",      "numeric", "numeric", "numeric",  "numeric", "character" , "character" , "numeric" ,  "numeric")
    TRUE
  } else if (any(argNotGiven(c(outFmtArg,outputColNames,outputColClasses)))){
    stop("For non-default formats, all three of `outFmtArg`, `outputColNames`, and `outputColClasses` must be provided.\n`outFmtArg` will take a form like 'general:name1,start1,end1,length1,strand1,name2,start2,end2,length2,strand2,id%,blastid%'.\nSee the manual page at https://www.bx.psu.edu/~rsharris/lastz/.")
  } else if (length(outputColNames)!=length(outputColClasses)) {
    stop("`outputColNames` and `outputColClasses` must be the same length character vectors. The latter must contain legitimate class names e.g. \"character\", \"numeric\", \"integer\" ...")
  } else {
    ce("Careful! You're using a non-default format. Confusing errors (if you get them) will likely be caused by:\n\t- Incorrect syntax in `outFmtArg`. See the manual at https://www.bx.psu.edu/~rsharris/lastz/.\n\t- Character classes in `outputColClasses` don't match the output corresponding to `outFmtArgs`.\n\t\t  Note by default lastz's 'id%' field appends a '%' sign to the output necessitating reading as a character string.\n\t- When `forceBioDToutput=F`, the absence of coercing LASTZ output to BioDT (more BLAST-like) conventions:\n\t\t  In constrast to the BLAST-style output of this function's default mode, you will get raw LASTZ output, which (e.g.) always reports start/end positions from lowest to highest (i.e., strand-agnostically), and;\n\t\t  Sequence output is always as per the subject sequence, i.e., with an exact + to - match, the sequences will be the subject, and the reverse-complemented query. (in BLAST it would reverse-complement the subject).")
    FALSE
  }

  if(addAlignments==TRUE){
    outFmtArg %<>% paste0(c(",text1,text2"))
    outputColNames %<>% c("sAlnSeq","qAlnSeq")
    outputColClasses %<>% c("character","character")
  }
  if(addPctIdGaps==TRUE){
    outFmtArg %<>% paste0(c(",id%"))
    outputColNames %<>% c("pctId")
    outputColClasses %<>% c("character")
  }
  if(addNgaps==TRUE){
    outFmtArg %<>% paste0(c(",cgap"))
    outputColNames %<>% c("nGaps")
    outputColClasses %<>% c("integer")
  }
  for(l in addColList){
    # FUTURE documentation / checks. Each list item is a list of 3: arg for lastz fmt string; name for R; class for R
    outFmtArg %<>% paste0(",",l[[1]])
    outputColNames %<>% l[[2]]
    outputColClasses %<>% l[[3]]
  }

  makeSfile <- F
  makeQfile <- F

  if(argNotGiven(subjectFname) & argGiven(subjectSeqDT)){
    is_seqDT(subjectSeqDT,objName=deparse(substitute(subjectSeqDT)),croak=T)
    makeSfile <- T
    writeFasta(subjectSeqDT,subjectFname<-tempfile(fileext = ".fasta"))
  }
  if(argNotGiven(queryFname) & argGiven(querySeqDT)){
    is_seqDT(querySeqDT,objName=deparse(substitute(querySeqDT)),croak=T)
    makeQfile <- T
    writeFasta(querySeqDT,queryFname<-tempfile(fileext = ".fasta"))
  }


  la <- ldtply(subjectFname,function(sFn){ #sseqfile # dev sFn<-subjectFname[1];
    # lastz's unpopularity derives purely from stupid things like its baffling inability to deal with a multi-sequence fasta
    sFaiFname <- paste0(subjectFname,".fai")
    sFai <- if(file.exists(sFaiFname)){
      readFai(sFaiFname)
    } else {
      getFai(subjectFname,samtoolsBinary=samtoolsBinary)
    }

    ldtply(sFai$seqId,function(ssName){ #sseq  # dev ssName<-sFai$seqId[1];
      ldtply(queryFname,function(qFn){ #qseqfile # dev qFn<-queryFname[1];
        qFaiFname <- paste0(queryFname,".fai")
        qFai <- if(file.exists(sFaiFname)){
          readFai(qFaiFname)
        } else {
          getFai(queryFname,samtoolsBinary=samtoolsBinary)
        }

        ldtply(qFai$seqId,function(qsName){ #qseq # dev qsName<-qFai$seqId[1];
          # align sseq in sseqfile to qseq in qseqfile
          #  ssName<-
          sSubsetArg <- paste0("\\[subset=<(echo \"",ssName,"\")\\]")
          qSubsetArg <- paste0("\\[subset=<(echo \"",qsName,"\")\\]")
          cmd <- paste0(lastzBinary," ",sFn,sSubsetArg," ",qFn,qSubsetArg," --format=",outFmtArg," ",otherLastzArgs)#,"\["sSubset,sRange,qFile,qSubset,qRange,outFmtArgs,otherArgs)
          ce(cmd)
          #fread(cmd=cmd,col.names=outputColNames)
          tmp <- fread(cmd=cmd,col.names=outputColNames,colClasses=outputColClasses)
          if(!makeSfile){ tmp[,sFastaFname:=sFn][] }
          if(!makeQfile){ tmp[,qFastaFname:=qFn][] }
          tmp
        })
      })
    })
  })

  if(makeBioDTOutput | forceBioDToutput){
    if(hasCol(la,"sAlnSeq")){ la[sStrand!=qStrand,sAlnSeq:=rc(sAlnSeq)] }
    if(hasCol(la,"qAlnSeq")){ la[sStrand!=qStrand,qAlnSeq:=rc(qAlnSeq)] }
    if(hasCol(la,"sStart") & hasCol(la,"sEnd")){ la[sStrand!=qStrand,c("sStart","sEnd"):=.(sEnd,sStart)]  }
    if(hasCol(la,"sStrand")) { la[,sStrand:=NULL] }
    if(hasCol(la,"qStrand")) { la[,qStrand:=NULL] }
    if(hasCol(la,"pctId")) { la[,pctId:=as.numeric(sub("%","",pctId))] }
  }

  if(argGiven(outFastaFname)){
    write.table(la,outFastaFname,row.names=F,sep="\t",quote=F)
  }

  if(makeSfile){ unlink(Sys.glob(paste0(subjectFname,".*"))) }
  if(makeQfile){ unlink(queryFname) }

  return(la[])
}
