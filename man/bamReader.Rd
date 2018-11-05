\name{bamReader}
\alias{bamReader}
\title{bamReader}
\description{The bamReader function constructs and returns an S4 object of class bamReader. The bamReader object represents a connection to an opened BAM-file. When an indexname is given or idx=TRUE, the function tries to load an existing BAM-index.}
\usage{bamReader(filename,indexname,idx=FALSE,verbose=0)}
\arguments{
  \item{filename}{Filename of BAM-file to be opened for reading.}
  \item{indexname}{Optional: Name of BAM-index file}
  \item{idx}{Logical}
  \item{verbose}{Numeric: Quantifies the extent of textual feedback (levels: 0,1,2). }
}
\author{Wolfgang Kaisers}
\examples{
bam<-system.file("extdata","accepted_hits.bam",package="rbamtools")
reader<-bamReader(bam)
align<-getNextAlign(reader)
name(align)
position(align)
bamClose(reader)
# Open reader and initialize BAM-index
reader<-bamReader(bam,idx=TRUE)
index.initialized(reader)
}
\keyword{bamReader}
