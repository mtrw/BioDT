#' @export
timPalettePresets <- function(name,alpha="ff"){
  list(
    wheel     = timPalette(paste0(c("#d41313","#f7760c","#e3c607","#51a321","#0fbfae","#0f3fa6","#871b53"),alpha)),
    yg     = timPalette(paste0(c("#EEDD00FF","#009933FF"),alpha)),
    mclaren       = timPalette(paste0(c("#232526","#f7760c","#00daef"),alpha)),
    ferrari   = timPalette(paste0(c("#151515","#f02929"),alpha)),
    mercedes      = timPalette(paste0(c("#d0d0d0","#191919","#76dfc6"),alpha)),
    redbull        = timPalette(paste0(c("#fec52a","#172a52","#ce0815"),alpha)),
    alphatauri        = timPalette(paste0(c("#0b2945","#e5e8ef","#f40000"),alpha)),
    astonmartin        = timPalette(paste0(c("#02716c","#e7f1f6","#040300"),alpha)),
    alfaromeo        = timPalette(paste0(c("#3e373e","#eff0f2","#cd0a24"),alpha)),
    sauber        = timPalette(paste0(c("#111f28","#00d312"),alpha)),
    alpine        = timPalette(paste0(c("#fbfeff","#ed98dd","#1766e0ff","#000010"),alpha))
  )[[name]]
}
