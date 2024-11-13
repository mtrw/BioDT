aScriptPrimer3generic <- function(){ #why is this a function? pfft
  "'

  BEGIN{
    go=0
    OFS=\"\\t\"
    print \"primerPairId\",\"end\",\"penalty\",\"seq\",\"pos\",\"length\",\"tm\",\"gc\",\"selfAnyTh\",\"selfEndTh\",\"hairpinTh\",\"endStability\",\"pairPenalty\",\"pairComplAnyTh\",\"pairComplEndTh\",\"pairProductSize\",\"pairProductTm\"
  }

  $1~/SEQUENCE_ID/ {seqId = $2}

  1 {split($1,a,\"_\")}

  a[3]!~\"^[0-9]+$\" {next}

  1 {go++}

  0 {print $1\"\\n----this:\",thisPrimerPairId,\"--last:\",lastPrimerPairId,\"----------------\"}

  1 {
    lastPrimerPairId=thisPrimerPairId
    thisPrimerPairId=a[3]
  }

  $1~/_LEFT_/ {left = 1}
  $1~/_RIGHT_/ {right = 1}
  $1~/_INTERNAL_/ {internal = 1}

  $1~/PRIMER_LEFT_._PENALTY/ {leftPenalty=$2}
  $1~/PRIMER_RIGHT_._PENALTY/ {rightPenalty=$2}
  $1~/PRIMER_INTERNAL_._PENALTY/ {internalPenalty=$2}

  $1~/PRIMER_LEFT_._SEQUENCE/ {leftSeq=$2}
  $1~/PRIMER_RIGHT_._SEQUENCE/ {rightSeq=$2}
  $1~/PRIMER_INTERNAL_._SEQUENCE/ {internalSeq=$2}

  $1~/PRIMER_LEFT_.$/ {split($2,b,\",\"); leftPos=b[1]; leftLength=b[2];}
  $1~/PRIMER_RIGHT_.$/ {split($2,b,\",\"); rightPos=b[1]; rightLength=b[2];}
  $1~/PRIMER_INTERNAL_.$/ {split($2,b,\",\"); internalPos=b[1]; internalLength=b[2];}

  $1~/PRIMER_LEFT_._TM/ {leftTm=$2}
  $1~/PRIMER_RIGHT_._TM/ {rightTm=$2}
  $1~/PRIMER_INTERNAL_._TM/ {internalTm=$2}

  $1~/PRIMER_LEFT_._GC_PERCENT/ {leftGC=$2}
  $1~/PRIMER_RIGHT_._GC_PERCENT/ {rightGC=$2}
  $1~/PRIMER_INTERNAL_._GC_PERCENT/ {internalGC=$2}

  $1~/PRIMER_LEFT_._SELF_ANY_TH/ {leftSelfAnyTh=$2}
  $1~/PRIMER_RIGHT_._SELF_ANY_TH/ {rightSelfAnyTh=$2}
  $1~/PRIMER_INTERNAL_._SELF_ANY_TH/ {internalSelfAnyTh=$2}

  $1~/PRIMER_LEFT_._SELF_END_TH/ {leftSelfEndTh=$2}
  $1~/PRIMER_RIGHT_._SELF_END_TH/ {rightSelfEndTh=$2}
  $1~/PRIMER_INTERNAL_._SELF_END_TH/ {internalSelfEndTh=$2}

  $1~/PRIMER_LEFT_._HAIRPIN_TH/ {leftHairpinTh=$2}
  $1~/PRIMER_RIGHT_._HAIRPIN_TH/ {rightHairpinTh=$2}
  $1~/PRIMER_INTERNAL_._HAIRPIN_TH/ {internalHairpinTh=$2}

  $1~/PRIMER_LEFT_._END_STABILITY/ {leftEndStability=$2}
  $1~/PRIMER_RIGHT_._END_STABILITY/ {rightEndStability=$2}
  $1~/PRIMER_INTERNAL_._END_STABILITY/ {internalEndStability=$2}

  $1~/PRIMER_PAIR_._PENALTY/ {primerPairPenalty=$2}
  $1~/PRIMER_PAIR_._COMPL_ANY_TH/ {primerPairComplAnyTh=$2}
  $1~/PRIMER_PAIR_._COMPL_END_TH/ {primerPairComplEndTh=$2}
  $1~/PRIMER_PAIR_._PRODUCT_SIZE/ {primerPairProductSize=$2}
  $1~/PRIMER_PAIR_._PRODUCT_TM/ {primerPairProductTm=$2}

  1 {primerPairId = seqId\"_\"thisPrimerPairId}

  (go>1) && (thisPrimerPairId!=lastPrimerPairId) {

    if(left){ print primerPairId,\"left\",leftPenalty,leftSeq,leftPos,leftLength,leftTm,leftGC,leftSelfAnyTh,leftSelfEndTh,leftHairpinTh,leftEndStability,primerPairPenalty,primerPairComplAnyTh,primerPairComplEndTh,primerPairProductSize,primerPairProductTm }

    if(right){ print primerPairId,\"right\",rightPenalty,rightSeq,rightPos,rightLength,rightTm,rightGC,rightSelfAnyTh,rightSelfEndTh,rightHairpinTh,rightEndStability,primerPairPenalty,primerPairComplAnyTh,primerPairComplEndTh,primerPairProductSize,primerPairProductTm }

    if(internal){ print primerPairId,\"internal\",internalPenalty,internalSeq,internalPos,internalLength,internalTm,internalGC,internalSelfAnyTh,internalSelfEndTh,internalHairpinTh,internalEndStability,primerPairPenalty,primerPairComplAnyTh,primerPairComplEndTh,primerPairProductSize,primerPairProductTm }

    left = 0; right = 0; internal = 0;
  }

  END {

    if(left){ print primerPairId,\"left\",leftPenalty,leftSeq,leftPos,leftLength,leftTm,leftGC,leftSelfAnyTh,leftSelfEndTh,leftHairpinTh,leftEndStability,primerPairPenalty,primerPairComplAnyTh,primerPairComplEndTh,primerPairProductSize,primerPairProductTm }

    if(right){ print primerPairId,\"right\",rightPenalty,rightSeq,rightPos,rightLength,rightTm,rightGC,rightSelfAnyTh,rightSelfEndTh,rightHairpinTh,rightEndStability,primerPairPenalty,primerPairComplAnyTh,primerPairComplEndTh,primerPairProductSize,primerPairProductTm }

    if(internal){ print primerPairId,\"internal\",internalPenalty,internalSeq,internalPos,internalLength,internalTm,internalGC,internalSelfAnyTh,internalSelfEndTh,internalHairpinTh,internalEndStability,primerPairPenalty,primerPairComplAnyTh,primerPairComplEndTh,primerPairProductSize,primerPairProductTm }

  }

  '"
}

aScriptFasta2tbl <- '{sub(/\\r/,"")} /^>/ { if(NR>1){printf("\\n")}; printf substr($1,2,length($0))"\t"; next; } 1 {printf $1}'

aScriptCdhitClstr2tbl <- '/^>Cluster/ {c=$2; next}   {sub(">","",$3); sub("\\\\.\\\\.\\\\.$","",$3); print $3,c}'

aScriptTrf2tbl <- '/^@/ {n=substr($1,2); next} {print n,$0}'

aScriptUclust2tbl <- '$1=="H" || $1=="C" {print $9"\\t"$2}'

