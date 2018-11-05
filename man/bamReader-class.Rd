\name{bamReader-class}
\Rdversion{1.1}
\docType{class}
\alias{bamReader-class}
\alias{bamClose,bamReader-method}
\alias{bamCopy}
\alias{bamCopy-methods}
\alias{bamCopy,bamReader-method}
\alias{bamSave,bamReader-method}
\alias{bamSort}
\alias{create.index}
\alias{create.index-methods}
\alias{create.index,bamReader-method}
\alias{getHeader,bamReader-method}
\alias{getHeaderText,bamReader-method}
\alias{getNextAlign,bamReader-method}
\alias{getRefCoords}
\alias{getRefCoords-methods}
\alias{getRefCoords,bamReader-method}
\alias{getRefCount}
\alias{getRefCount-methods}
\alias{getRefCount,bamReader-method}
\alias{getRefData}
\alias{getRefData-methods}
\alias{getRefData,bamReader-method}
\alias{filename,bamReader-method}
\alias{bamSort,bamReader-method}
\alias{load.index}
\alias{load.index,bamReader-method}
\alias{initialize,bamReader-method}
\alias{index.initialized}
\alias{index.initialized-methods}
\alias{index.initialized,bamReader-method}
\alias{isOpen,bamReader-method}
\alias{rewind}
\alias{rewind-methods}
\alias{rewind,bamReader-method}

\title{Class \code{"bamReader"}}
\description{A bamReader object encapulates functionality for reading of bamAligns from a BAM-file. It optionally contains a pointer to a bam-index structure which allows fast access to aligns that overlap a specified region (random access). The index is loaded via the load.index function.}
\section{Objects from the Class}{Objects can be created by calls of the form \code{reader<-bamReader(filename)}.}
\section{Slots}{
  \describe{
    \item{\code{filename}:}{Character denoting name of BAM-file.}
    \item{\code{reader}:}{External pointer to opened BAM-file.}
    \item{\code{startpos}:}{Numeric value returned by bam_tell immediately after opening BAM-file. It is used as target position for rewinding by bam_seek.}
    \item{\code{index}:}{External pointer to BAM index (used for retrieving bamRange objects from BAM-file.}
  }
}
\section{Methods}{
  \describe{
    \item{bamClose}{\code{signature(object = "bamReader")}: Closes connection to BAM-file. }
    \item{bamCopy}{\code{signature(object = "bamReader",writer="bamWriter",refids,verbose)}: Copies all aligns for given refids from reader to writer. When no refids are given, aligns for all refids are copied. refids refer to Reference-ID's returned by 'getRefData' (ID column). Missing matches (from refids into Reference-ID's) will terminate the function with an error.}
    \item{bamSave}{\code{signature(object = "bamReader",writer="bamWriter")}: Copies all aligns from reader to writer.}
    \item{create.index}{\code{signature(object = "bamReader", filename="character")}: Creates an index file for opened BAM-file. Therfore the BAM-file must be position-sorted. filename determines the name of the new BAM-Index-file (usually bam_filename.bam.idx) }
    \item{filename}{\code{signature(object = "bamReader")}: Returns filename of opened BAM-file. }
    \item{getHeader}{\code{signature(object = "bamReader")}: Returns object of class bamHeader which contains binary representation of bam-header data. }
    \item{getHeaderText}{\code{signature(object = "bamReader")}: Returns object of class bamHeaderText which contains textual representation of bam-header data. }
    \item{getNextAlign}{\code{signature(object = "bamReader")}: Returns object of class bamAlign which contains data of next Align from file. When EOF is reached the function returns NULL. }
    \item{getNextChunk}{\code{signature(object = "bamReader", n="integer")}: Returns object of class bamRange which contains data of next 'n' Aligns from file. When EOF is reached the function returns a range of size zero. }
   \item{getRefCoords}{\code{signature(object="bamReader",sn="character")}: Helper function takes a sequence name and returns coordinates of entire reference sequence for usage with bamRange, gapList or siteList function. The function returns a vector of length 3. The vector elements are named "refid","start","stop". }
    \item{getRefCount}{\code{signature(object = "bamReader")}: Returns number of reference sequences. }
    \item{getRefData}{\code{signature(object = "bamReader")}: Returns data frame which contains three columns: For each reference sequence, the corresponding row contains the Reference-ID (1st column, refID), the Reference name (2nd column, refName) and the length of the Reference sequence (3rd column, refLength).
}
    \item{isOpen}{\code{signature(object = "bamReader")}: Returns TRUE when file connection is open. }
    \item{rewind}{\code{signature(object = "bamReader")}: Resets current file position. The subsequent call to getNextAlign will return the first align in the BAM-file.}
    \item{bamSort}{\code{signature(object = "bamReader",prefix="character",byName=FALSE,maxmem=1e+9)}: BAM-files must be sorted before an index can be created. During this routine, some content of the source file is copied into the working memory, sorted and written into temporary files. Finally, the sorted files are merged into a single target. The name of the sorted BAM-file will be [prefix].bam.
Number and size of temporary files depend on source file size and used working memory (maxmem). Small
maxmem values result in large numbers of temporary files. The minimum maxmem value is 100000000. Smaller maxmem
values cause an error. }
    \item{initialize}{\code{signature(.Object = "bamReader")}: Initializes object and opens BAM-file for reading. }
	 }
}
\references{
\verb{The SAM Format Specification (v1.4-r985). The SAM Format Specification Working Group. September 7, 2011.}
\url{http://www.samtools.sourceforge.net/SAM1.pdf}
}
\author{Wolfgang Kaisers}
\examples{
library(rbamtools)
bam<- system.file("extdata", "accepted_hits.bam", package="rbamtools")
idx<- system.file("extdata", "accepted_hits.bam.bai", package="rbamtools")
# Open BAM-file for reading
reader<-bamReader(bam)
getNextAlign(reader)
rewind(reader)
getNextAlign(reader)
isOpen(reader)
# Create and load new index
idx<-"index.bam.bai"
create.index(reader,idx)
load.index(reader,idx)
index.initialized(reader)
# Load package provided index
load.index(reader,idx)
index.initialized(reader)
# Read align
align<-getNextAlign(reader)
coords<-as.integer(c(0,0,249250621))
range<-bamRange(reader,coords)
align<-getNextAlign(range)
# Copy all aligns
writer<-bamWriter(getHeader(reader),"newFile1.bam")
bamSave(reader,writer)
bamClose(writer)
# reader2<-bamReader("newFile1.bam")
# bamSort(reader,"newFile1s")
# Copy align for Reference-ID '0'
writer<-bamWriter(getHeader(reader),"newFile2.bam")
bamCopy(reader,writer,0)
bamClose(writer)
# Closing reader
bamClose(reader)
}
\keyword{classes}
\keyword{bamReader}
