\name{compressSpliceMaps}
\alias{compressSpliceMaps}
\title{ Compress Splice Maps }
\description{
  compress Splice Maps for later use with spliced alignments.
}
\usage{
compressSpliceMaps( outPath = NULL, spliceMaps, spliceMapPrefix="SpliceMap")
}
\arguments{
  \item{outPath}{ path for writing the compressed spliceMap files. }
  \item{spliceMaps}{ a vector of filenames of existing splice maps in text format, as
  created by \code{\link{buildSpliceJunctionFiles}}. }
  \item{spliceMapPrefix}{ a character string prefix to combine with the SeqID for naming
  the new compressed splice junction files. }
}
\details{
  This turns one potentially very large splice map into fast compressed maps for each
  SeqID.
}
\value{
  Invoked for its side effect of writing files to 'outPath', one for each SeqID in that
  organism.
}
\author{ Bob Morrison }
\seealso{ 
  \code{\link{spliceToAlign}} { for using these compressed slice maps to resolve splice
  alignments. }
}
