\name{BAM_Files}
\alias{BAM_Files}
\alias{BAM.sort}
\alias{BAM.index}
\alias{BAM.verifySorted}
\alias{BAM.merge}
\alias{BAM.indexFASTA}
\alias{BAM.fieldTable}
\alias{cleanupBAMfiles}
\alias{SAM2BAM}
\title{
Utilities for Manipulating BAM Files.
}
\description{
A variety of functions for operating on BAM and SAM files, as may be needed by various 
DuffyNGS pipeline tools or SAMTOOLS utilities.
}
\usage{
BAM.sort(file, what = c("position", "readID"), index = TRUE, memory = 2.147e+09, threads = 2, verbose = TRUE)
BAM.index(file, verbose = TRUE)
BAM.verifySorted(file, index = TRUE, threads = 2, verbose = TRUE)
BAM.merge(files, newfile, index = TRUE, verbose = TRUE)
BAM.indexFASTA(file, verbose = TRUE)
BAM.fieldTable(file, field, chunkSize = 100000, maxReads = NULL, verbose = TRUE)
cleanupBAMfiles( path = "results")
SAM2BAM(file, verbose = TRUE, delete.SAM = TRUE)
}
\arguments{
  \item{file, files}{
Character string full pathname of one existing BAM file.  Or a character vector of 2+ BAM files.
}
  \item{what}{
Character vector of named fields in the BAM file to do the sorting on.
}
  \item{index}{
Logical, should the sort BAM file have its sort index (.bai) file created as well.
}
  \item{memory}{
Size of local memory in bytes for sort scratch space.
}
  \item{threads}{
Allows use of multithreading to improve sort time performance.
}
  \item{newfile}{
Filename for the newly created BAM file.
}
  \item{field}{
Character string of a single field in the BAM file, to be the value that gets inspected.
}
  \item{chunkSize}{
Numeric buffer size, number of BAM file alignments to read into memory at a time.
}
  \item{maxReads}{
Numeric, maximum number of BAM file alignments to inspect.
}
  \item{path}{
Character string of a top level results folder, that contains various alignment pipeline result subfolders.
}
}
\details{
These utility functions are designed to do typical error checking about file existence, etc., so higher
level pipeline tools don't have to.  Most are just simple wrappers around system calls to SAMTOOLS functions.
}
\value{
For \code{BAM.sort}, when successful, the filename of the created sorted BAM file.

For \code{BAM.index}, when successful, the filename of the created BAM index file.

For \code{BAM.verifySorted}, when successful, the filename of the BAM file.

For \code{BAM.merge}, when successful, the filename of the new merged BAM file.

For \code{BAM.indexFASTA}, when successful, the filename of the created FASTA index file.

For \code{BAM.fieldTable}, a data frame that tallies the counts and percentage for each value
observed for that field in the BAM file.

For \code{cleanupBAMfiles}, a brief summary of file counts and bytes of disk space deleted.

For \code{SAM2BAM}, when successful, the filename of the created BAM file.
}
\author{
Bob Morrison
}
\note{
There are many other low level functions for operating on the contents of BAM files, either
one alignment at a time or in buffered chunks of alignments.  For a brief introduction, 
see \code{\link{bamReader}}.
}
