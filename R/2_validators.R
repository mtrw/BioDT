#' @export
hasSeqFname <- makeValidator(seqFname_valSpecs)
#' @export
hasSeqId <- makeValidator(hasSeqId_valSpecs)
#' @export
hasStrand <- makeValidator(hasStrand_valSpecs)
#' @export
hasName <- makeValidator(hasName_valSpecs)
#' @export
hasScore <- makeValidator(hasScore_valSpecs)
#' @export
hasStartEnd <- makeValidator(hasStartEnd_valSpecs)
#' @export
isSeqDT <- makeValidator(seq_valSpecs)
#' @export
isSeqDT_withSeqFname <- function(validateMe,...){ isSeqDT(validateMe,...) & hasSeqFname(validateMe,...) }
#' @export
isCoordDT <- makeValidator(coord_valSpecs,function(coordDT){coordDT[,all(start<=end)]})
#' @export
isCoordDT_outputsFileStrand <- function(validateMe,...){ isCoordDT(validateMe,...) & hasScore(validateMe,...) & hasName(validateMe,...) & hasStrand(validateMe,...) }
