
vcf=$1

fSep=$(cat ${vcf} | grep -v "^#" | head -1 | cut -f9 | awk '{print substr($1,3,1)}')

nFs=$(cat ${vcf} | grep -v "^#" | head -1 | cut -f9 | tr -cd $fSep | wc -c | awk '{print $1+1}')

gtSep=$(cat ${vcf} | grep -v "^#" | head -1 | cut -f10 | sed "s/${fSep}.*//" | sed 's/^[[:digit:]]//' | awk '{print substr($1,1,1)}')

ploidy=$(cat ${vcf} | grep -v "^#" | head -1 | cut -f10 | sed "s/${fSep}.*//" | cut -f10 | tr -cd ${gtSep} | wc -c | awk '{print $1+1}')

nGts=$(cat ${vcf} | grep "#CHROM" | awk '{print NF-9; exit}')

cutGtCmd=$(echo 1 | awk -v nFs=${nFs} -v nGts=${nGts} 'BEGIN {   for(i=0;i<=((nGts-1)*nFs);i+=nFs){ printf (i==0)?"":",";  printf i+1; }  }')
cutNotGtCmd=$(echo 1 | awk -v nFs=${nFs} -v nGts=${nGts} 'BEGIN {   for(i=1;i<=(nGts*nFs);i++){ if((i)%nFs){ printf (i==1)?"":","; printf i+1 } }  }')

#echo $fSep $nFs $gtSep $ploidy $nGts $cutGtCmd -- $cutNotGtCmd
echo $nFs $ploidy $nGts 

p1_9=$(mktemp -u); mkfifo ${p1_9}
p10__split=$(mktemp -u); mkfifo ${p10__split}
pGt=$(mktemp -u); mkfifo ${pGt}
pnGt=$(mktemp -u); mkfifo ${pnGt}

cat ${vcf} | grep -v "^#" | cut -f 1-9 > ${p1_9} &
cat ${vcf} | grep -v "^#" | cut -f 10- | sed 's/'"${fSep}"'/\t/g' > ${p10__split} &
cat ${p10__split} | tee ${pnGt} > ${pGt} &

paste ${p1_9} <(cat ${pGt} | cut -f "${cutGtCmd}" | sed 's='${gtSep}'=\t=g') <(cat ${pnGt} | cut -f "${cutNotGtCmd}")

exit 0;
