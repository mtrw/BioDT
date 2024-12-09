% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers_blast.R
\name{readBlastTable}
\alias{readBlastTable}
\title{A convenience function to read in tabular BLAST output in the format requested (under the bonnet) by this wrapper}
\usage{
readBlastTable(
  fileName,
  outputColNames = c("sSeqId", "qSeqId", "matchLength", "sStart", "qStart", "sEnd",
    "qEnd", "pctId_noGaps", "eValue", "bitscore"),
  outputColClasses = c("character", "character", "integer", "numeric", "numeric",
    "numeric", "numeric", "numeric", "numeric", "numeric")
)
}
\description{
A convenience function to read in tabular BLAST output in the format requested (under the bonnet) by this wrapper
}
