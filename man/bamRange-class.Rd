\name{bamRange-class}
\Rdversion{1.1}
\docType{class}
\alias{bamRange-class}
\alias{as.data.frame.bamRange}
\alias{getNextAlign,bamRange-method}
\alias{coerce,bamRange,data.frame-method}
\alias{getPrevAlign}
\alias{getPrevAlign-methods}
\alias{getPrevAlign,bamRange-method}
\alias{initialize,bamRange-method}
\alias{insertPastCurrent}
\alias{insertPastCurrent-methods}
\alias{insertPastCurrent,bamRange-method}
\alias{insertPreCurrent}
\alias{insertPreCurrent-methods}
\alias{insertPreCurrent,bamRange-method}
\alias{moveCurrentAlign}
\alias{moveCurrentAlign-methods}
\alias{moveCurrentAlign,bamRange-method}
\alias{pop_back}
\alias{pop_back-methods}
\alias{pop_back,bamRange-method}
\alias{pop_front}
\alias{pop_front-methods}
\alias{pop_front,bamRange-method}
\alias{push_back}
\alias{push_back-methods}
\alias{push_back,bamRange-method}
\alias{push_front}
\alias{push_front-methods}
\alias{push_front,bamRange-method}
\alias{size}
\alias{size-methods}
\alias{size,bamRange-method}
\alias{stepNextAlign}
\alias{stepNextAlign-methods}
\alias{stepNextAlign,bamRange-method}
\alias{stepPrevAlign}
\alias{stepPrevAlign-methods}
\alias{stepPrevAlign,bamRange-method}
\alias{rewind,bamRange-method}
\alias{writeCurrentAlign}
\alias{writeCurrentAlign-methods}
\alias{writeCurrentAlign,bamRange-method}
\alias{alignQual}
\alias{alignQual-methods}
\alias{alignQual,bamRange-method}
\alias{alignSeq}
\alias{alignSeq-methods}
\alias{alignSeq,bamRange-method}
\alias{readQual}
\alias{readQual-methods}
\alias{readQual,bamRange-method}
\alias{readSeq}
\alias{readSeq-methods}
\alias{readSeq,bamRange-method}
\alias{readID}
\alias{readID-methods}
\alias{readID,bamRange-method}
\alias{refID}
\alias{refID-methods}
\alias{refID,bamRange-method}
\alias{position}
\alias{position-methods}
\alias{position,bamRange-method}
\alias{reverseStrand}
\alias{reverseStrand-methods}
\alias{reverseStrand,bamRange-method}


\title{Class \code{"bamRange"}}
\description{bamRange represents a double linked list of bamAlign objects which overlap with a defined region in a BAM-file. The bamRange-function retrieves all alignments in the depicted Range from BAM-File into a bamRange object. A bamRange object maintains a double-linked list of aligns and a pointer to a current align.

rewind shifts the current align before the first align (=past the last align). When there is no align, nothing happens. The functions writeCurrentAlign, insertPreCurrent, insertPastCurrent can be used to write bamAlign data to bamRange object. The function writeCurrentAlign overwrites the data of the current bamRange align. When no current align is set (i.e. ==NULL), nothing is written. The functions insertPreCurrent and insertPastCurrent insert an bamAlign object into the bamRange list. When no current align is set, insertPreCurrent inserts the object as the last object and insertPastCurrent inserts the object as the first align object in the bamRange list.

push_back (push_front) appends a bamAlign before the first (or behind the last) element of bamRange. pop_front (pop_back) removes the first (last) bamAlign from bamRange.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{range<-bamRange(reader, coords)}.}
\section{Slots}{
  \describe{
    \item{\code{range}:}{External pointer to double linked list of bamAligns.}
  }
}
\section{Methods}{
  \describe{
    \item{as.data.frame}{\code{signature(x = "bamRange")}: Returns data.frame representation of aligns. }
    \item{coerce}{\code{signature(from = "bamRange", to = "data.frame")}: Coercion of bamRange to data.frame.}
    \item{bamSave}{\code{signature(object = "bamRange")}: Saves aligns stored in this list to BAM-file via a bamWriter object. }
    \item{getNextAlign}{\code{signature(object = "bamRange")}: Returns next align from current position and shifts current position to next one. }
    \item{getPrevAlign}{\code{signature(object = "bamRange")}: Returns previous align from current position and shifts current position to previous one. }
    \item{initialize}{\code{signature(.Object = "bamRange")}: Initializes bamRange object. }
    \item{insertPastCurrent}{\code{signature(object = "bamRange")}: Inserts align past current position into list. }
    \item{insertPreCurrent}{\code{signature(object = "bamRange")}: Insert align before current position into list. }
    \item{pop_back}{\code{signature(object = "bamRange")}: Removes last align from list. }
    \item{pop_front}{\code{signature(object = "bamRange")}: Removes first align from list. }
    \item{push_back}{\code{signature(object = "bamRange")}: Adds align at the end of the list. }
    \item{push_front}{\code{signature(object = "bamRange")}: Adds align at the front of the list. }
    \item{rewind}{\code{signature(object = "bamRange")}: Shifts current align to position before first align. }
    \item{size}{\code{signature(object = "bamRange")}: Returns number of aligns in list. }
    \item{writeCurrentAlign}{\code{signature(object = "bamRange")}: Overwrites current align with given align. }

    \item{alignQual}{\code{signature(object = "bamRange")}: Returns quality scores of all aligned reads (one for each query-base).  Will be reverse of raw if aligned to reverse strand. }
    \item{alignSeq}{\code{signature(object = "bamRange")}: Returns base sequence of all aligned reads.  Will be reverse complement of raw if aligned to reverse strand.}
    \item{getTag}{\code{signature(object = "bamRange", tag="character")}: Retrieves the value of the given tag as a character vector, with "" if 'tag' not found.}
    \item{getAllTags}{\code{signature(object = "bamRange", sep="character")}: Retrieves the complete descriptors of all tags as a character vector, with with 'sep' inserted between tags.}
    \item{readQual}{\code{signature(object = "bamRange")}: Returns quality scores of all raw reads prior to alignment (one for each query-base). }
    \item{readSeq}{\code{signature(object = "bamRange")}: Returns base sequence of all raw reads prior to alignment. }
    \item{readID}{\code{signature(object = "bamRange")}: Returns all read Identifiers. }
    \item{refID}{\code{signature(object = "bamRange")}: Returns all reference sequence Identifiers. }
   \item{position}{\code{signature(object = "bamRange")}: Returns all align positions (on Ref-Sequence). }
    \item{reverseStrand}{\code{signature(object = "bamRange")}: Gets flag item denoting alignment to the reverse strand. }
	 }
}
\author{Wolfgang Kaisers}
\examples{
library(rbamtools)
bam <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
idx<-paste(bam,"bai",sep=".")
# Open BAM file
reader<-bamReader(bam)
# create.index(reader,idx)
# Load BAM index file
load.index(reader,idx)
index.initialized(reader)   # Should return 'TRUE'

# +++++++++++++++++++++++++++++++++++++++++++++++
#  Find appropriate refid (=ID)
#  Returns a data.frame with three columns:
#  ID=refid, SN=Sequence Name, LN=Sequence length

rdf<-getRefData(reader)
head(rdf)

# +++++++++++++++++++++++++++++++++++++++++++++++
# The sequence length (LN) also determines valid
# range for start and stop coordinates
# Invalid refid-, start- or stop-coordinates will
# release an error.
# coords: refid=0, start=0, stop=249250621
# +++++++++++++++++++++++++++++++++++++++++++++++

coords<-as.integer(c(0,0,249250621))
range<-bamRange(reader,coords)
size(range)
range<-bamRange(reader,coords,complex=TRUE)
dfr<-as.data.frame(range)
size(range)
align<-getNextAlign(range)
cigarData(align)
rewind(range)
# while(!is.null(align))
# {
#   print(position(align))
#   align<-getNextAlign(range)
# }
bamClose(reader)
}
\keyword{classes}
\keyword{bamReader}
\keyword{bamAlign}
