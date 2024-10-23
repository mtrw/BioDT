# BioDT

Do all your bioinformatics from R, almost exclusively using `data.table`s and tools from that paradigm.

**WORK IN PROGRESS**

Represents common bioinformatics formats like fasta, bed, vcf, as `data.table`s with standardised columns, but which can be decorated or manipulated to the user's desire. Provides wrappers to seamlessly do common bioinformatics tasks that all read to and from `data.table`s. The user enjoys all the power and flexbility of data.table, and is spared the fiddliness and pain of jumping endlessly between primary processing in Linux and interactive exploration in R.

*Includes:*
- Wrappers for common linux bioinformatics tools for things like blast, lastz, bedtools, primer3, clustalo, cd-hit, tandemRepeatsFinder ... the list grows.
- Read/write functions (file ==> data.table; data.table ==> file) and converters (data.table ==> data.table).
- A very flexible and loose 'typing' system for data.tables, e.g. a 'seqDT' has columns "seqId" [character], "seq" ["character"], and as long as they have valid entries, all functions requiring a seqDT will work. A validation system checks data.tables for compliance with BioDT conventions and gives helpful feedback on how to fix any errors.
- **SOON** Interfaces e.g. to pluck out chunks of [bv]cf files as data.tables, using syntax like `vcf.interface["chromosome",start,end]`
- Various plotting and visualisation functions
- Various convenience functions for bioinformatics tasks that keep cropping up again and again (`alleleFrequency()`, `He()`[Effective Heterozygosity], `setRefAlleleMostCommon()`, ...), and just other very useful things (`mostCommonThing()`, )
- A very lightweight and powerful colour palette generator designed for layering custom plots in base::plot.
- **SOON** Documentation and vignettes.

BioDT is being developed by [Tim Rabanus-Wallace](https://safes.unimelb.edu.au/research/cropgem-lab#people) from the [CropGEM lab](https://safes.unimelb.edu.au/research/cropgem-lab) at the University of Melbourne.

*Example:*

\# ${PATH} must include recent versions of blast, bedtools, awk, ...

require(BioDT) # also loads data.table and magrittr

\# Read in some 'seed' sequences to study from a file
seedSeqs <- readFasta("data/seed_seqs.fasta")
\# find homologs in a reference genome
bl <- blast(subjectFname="data/genome.fasta",querySeqDT=seedSeqs)
\# Remind yourself ... what constitutes a valid alignmentDT?
is_alignmentDT(show=TRUE)
\# Plot results
bl[,{
  x <- pmean2(sStart,sEnd) # using BioDT::pmean2()
  y <- frank(sSeqId)
  col <- applyPalette( x=bitscore , colChain = palettePresets$mercedes ) #using BioDT's palette functions
  null_plot(x,y)
  points(x,y,col=col,pch=20,cex=2)
}]
\# Filter for 5 best hits per query
setorder(bl,qSeqId,bitscore)
bestHits <- bl[,last(.SD,5),by=.(qSeqId)] # using data.table::last()
\# Get coords of hits plus surrounding sequence
surround <- 500L
bestHitSurroundCoords <- alignmentDT2coordDT(bestHits,coordsOf = "subject")[,.(seqId,start=start-surround,end=end+surround)]
\# Get the seqs of the hits plus surrounding context
bestHitSurroundSeq <- extractSeq(coordDT=bestHitSurroundCoords,fastaFname = "data/genome.fasta")
\# Multiply align
msa <- MSA(bestHitSurroundSeq)
\# Plot the alignment with custom colours
plotNucleotideAlignment( msa[,.(seqId,seq=alnSeq)] )
