\name{bamRange}
\alias{bamRange}
\title{bamRange(bamReader, coordinates,complex=FALSE)}
\description{The bamRange function takes a bamReader object, a set of reference coordinates and the 'complex' argument and returns an instance of class 'bamRange'.
}
\usage{bamRange(reader=NULL, coords=NULL,complex=FALSE)}
\arguments{
	\item{reader}{An instance of bamReader. Must be opened and contain initialized index}
	\item{coords}{Integer vector of length 3: coords=c(refid,start,stop)}
	\item{complex}{A logical value (of length 1). Default value: FALSE}

}
\details{The method returns a list of bamAlign's from which overlap with the specified region. When complex is TRUE, the function only retrieves Aligns where nCigar>1 ('complex' aligns, e.g. 45M329N56M). When reader is NULL, an empty range-list ist constructed (can be filled with push_back). When complex is FALSE, the function retrieves all Aligns which fall into the given range.}
\value{An instance of class bamRange which can be accessed sequentially, modified or written to a BAM-file.}
\examples{
library(rbamtools)
bam <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
idx<-paste(bam,"bai",sep=".")
# Open BAM file
reader<-bamReader(bam)

# +++++++++++++++++++++++++++++++++++++++++++++++
#  Create empty range and fill with push_back
range<-bamRange()
for(i in 1:10)
{
  align<-getNextAlign(reader)
  push_back(range,align)
}
size(range)

# +++++++++++++++++++++++++++++++++++++++++++++++
#  
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
range<-bamRange(reader,coords)
size(range)
align<-getNextAlign(range)
cigarData(align)
range<-bamRange(reader,coords,complex=TRUE)
size(range)
align<-getNextAlign(range)
cigarData(align)
# while(!is.null(align))
# {
#   print(position(align))
#   align<-getNextAlign(range)
# }
bamClose(reader)
}
\keyword{bamRange}
\keyword{index}
