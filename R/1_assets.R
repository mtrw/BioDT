#' @export
palettePresets <- list(
  F1 = list(
    mclaren       = c("#232526","#fa870c","#00daef"),
    mclaren2       = c("#030506","#fa870c"),
    mclaren1995       = c("#f5f8ff","#f02929"),
    ferrari   = c("#151515","#f02929"),
    mercedes      = c("#191919","#76dfc6","#d0d0d0"),
    redbull        = c("#ffa800","#141823","#fe0b13"),
    alphatauri        = c("#0b2945","#f5f8ff","#dd1010"),
    RB             = c("#0a06b9","#e7e6eb","#d30303"),
    astonmartin        = c("#02716c","#e7f1f6","#040300"),
    alfaromeo        = c("#3e373e","#eff0f2","#cd0a24"),
    sauber        = c("#111f28","#00d312"),
    alpine        = c("#fbfeff","#ed98dd","#1766e0ff","#000010"),
    williams        = c("#113950","#032cc6","#04367e","#18fbfe")
  ),
  wheels = list(
    wheel= c("#d41313","#f7760c","#e3c607","#51a321","#0fbfae","#0f3fa6","#871b53"),
    fluoro=c("#FF0000","#fc8d05","#fcf805","#47fc05","#0ffcfc","#0f4efc","#e00ffc"),
    pastel=c("#d16868","#ce9a6f","#ddcc58","#82b565","#7cbab4","#6e88bf","#b77597")
  )
)


#' @export
codonTable <- d.t(
  codon = c("aaa", "aac", "aag", "aat", "aca", "acc", "acg",  "act", "aga", "agc", "agg", "agt", "ata", "atc", "atg", "att",  "caa", "cac", "cag", "cat", "cca", "ccc", "ccg", "cct", "cga",  "cgc", "cgg", "cgt", "cta", "ctc", "ctg", "ctt", "gaa", "gac",  "gag", "gat", "gca", "gcc", "gcg", "gct", "gga", "ggc", "ggg",  "ggt", "gta", "gtc", "gtg", "gtt", "taa", "tac", "tag", "tat",  "tca", "tcc", "tcg", "tct", "tga", "tgc", "tgg", "tgt", "tta",  "ttc", "ttg", "ttt"),
  CODON = c("AAA", "AAC", "AAG", "AAT", "ACA",  "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC",  "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG",  "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",  "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA",  "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAA", "TAC",  "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG",  "TGT", "TTA", "TTC", "TTG", "TTT"),
  aaName = c("Lys", "Asn",  "Lys", "Asn", "Thr", "Thr", "Thr", "Thr", "Arg", "Ser", "Arg",  "Ser", "Ile", "Ile", "Met", "Ile", "Gln", "His", "Gln", "His",  "Pro", "Pro", "Pro", "Pro", "Arg", "Arg", "Arg", "Arg", "Leu",  "Leu", "Leu", "Leu", "Glu", "Asp", "Glu", "Asp", "Ala", "Ala",  "Ala", "Ala", "Gly", "Gly", "Gly", "Gly", "Val", "Val", "Val",  "Val", "Stp", "Tyr", "Stp", "Tyr", "Ser", "Ser", "Ser", "Ser",  "Stp", "Cys", "Trp", "Cys", "Leu", "Phe", "Leu", "Phe"),
  aaSymbol = c("K",  "N", "K", "N", "T", "T", "T", "T", "R", "S", "R", "S", "I", "I",  "M", "I", "Q", "H", "Q", "H", "P", "P", "P", "P", "R", "R", "R",  "R", "L", "L", "L", "L", "E", "D", "E", "D", "A", "A", "A", "A",  "G", "G", "G", "G", "V", "V", "V", "V", "*", "Y", "*", "Y", "S",  "S", "S", "S", "*", "C", "W", "C", "L", "F", "L", "F")
)

#' @export
aaTable <- d.t(
  aaName   = c("Lys", "Asn", "Thr", "Arg", "Ser", "Ile", "Met", "Gln", "His",  "Pro", "Leu", "Glu", "Asp", "Ala", "Gly", "Val", "Stp", "Tyr",  "Cys", "Trp", "Phe"),
  aaSymbol = c("K", "N", "T", "R", "S", "I", "M", "Q", "H", "P", "L", "E",  "D", "A", "G", "V", "*", "Y", "C", "W", "F"),
  aaProp   = c("pos", "pol", "pol", "pos", "pol", "ali", "ali", "pol", "pos",  "pol", "ali", "neg", "neg", "ali", "ali", "ali", "STP", "aro",  "pol", "aro", "aro")
)

#' @export
ntCodesLegal <- c(
  c("-" ,"A", "C", "G", "T", "M", "R", "W", "S", "Y", "K", "V", "H",  "D", "B", "N"),
  c("A", "C", "G", "T", "M", "R", "W", "S", "Y", "K", "V", "H",  "D", "B", "N") %>% tolower()
)

#' @export
aaCodesLegal <- c(
  c("-" ,"K", "N", "T", "R", "S", "I", "M", "Q", "H", "P", "L", "E",  "D", "A", "G", "V", "Y", "C", "W", "F", "N", "*"),
  c("K", "N", "T", "R", "S", "I", "M", "Q", "H", "P", "L", "E",  "D", "A", "G", "V", "Y", "C", "W", "F", "N")  %>% tolower()
)

#' @export
ntCodesLegalNogap <- c(
  c("A", "C", "G", "T", "M", "R", "W", "S", "Y", "K", "V", "H",  "D", "B", "N"),
  c("A", "C", "G", "T", "M", "R", "W", "S", "Y", "K", "V", "H",  "D", "B", "N") %>% tolower()
)

#' @export
aaCodesLegalNogap <- c(
  c("K", "N", "T", "R", "S", "I", "M", "Q", "H", "P", "L", "E",  "D", "A", "G", "V", "*", "Y", "C", "W", "F", "N"),
  c("K", "N", "T", "R", "S", "I", "M", "Q", "H", "P", "L", "E",  "D", "A", "G", "V", "Y", "C", "W", "F", "N")  %>% tolower()
)

#' @export
ntCodesLegalUcOnly <- c(
  c("-" ,"A", "C", "G", "T", "M", "R", "W", "S", "Y", "K", "V", "H",  "D", "B", "N")
)

#' @export
ntCodesLegalUcOnly <- c(
  c("-","A", "C", "G", "T" , "M", "R", "W", "S", "Y", "K", "V", "H",  "D", "B", "N")
)

#' @export
ambig2nt <- d.t(
  nt=c(c("-","A", "C", "G", "T", "M", "R", "W", "S", "Y", "K", "V", "H",  "D", "B", "N"),c("A", "C", "G", "T", "M", "R", "W", "S", "Y", "K", "V", "H",  "D", "B", "N") %>% tolower),
  ambig=c(c("-","A", "C", "G", "T", "AC", "AG", "AT", "GC", "CT", "GT", "ACG", "ACT",  "AGT", "CGT", "ACGT"),c("A", "C", "G", "T", "AC", "AG", "AT", "GC", "CT", "GT", "ACG", "ACT",  "AGT", "CGT", "ACGT") %>% tolower)
)

#' @export
ambig2ntAmbigOnly <- d.t(
  nt=c(c("M", "R", "W", "S", "Y", "K", "V", "H",  "D", "B", "N"),c("M", "R", "W", "S", "Y", "K", "V", "H",  "D", "B", "N") %>% tolower),
  ambig=c(c("AC", "AG", "AT", "GC", "CT", "GT", "ACG", "ACT",  "AGT", "CGT", "ACGT"),c("AC", "AG", "AT", "GC", "CT", "GT", "ACG", "ACT",  "AGT", "CGT", "ACGT") %>% tolower)
)

#' @export
aaCodesLegalUcOnly <- c(
  c("-" ,"K", "N", "T", "R", "S", "I", "M", "Q", "H", "P", "L", "E",  "D", "A", "G", "V", "Y", "C", "W", "F", "N", "*")
)

#' @export
ntCodesLegalUcOnlyNogap <- c(
  c("A", "C", "G", "T", "M", "R", "W", "S", "Y", "K", "V", "H",  "D", "B", "N")
)

#' @export
aaCodesLegalUcOnlyNogap <- c(
  c("K", "N", "T", "R", "S", "I", "M", "Q", "H", "P", "L", "E",  "D", "A", "G", "V", "Y", "C", "W", "F", "N", "*")
)




