\name{alignToWig}
\alias{alignToWig}
\title{
Add BAM Alignments to Wiggle Track data objects.
}
\description{
Lower level function that adds alignments to wiggle tracks of read pileups for one sample.
}
\usage{
alignToWig( filein, reload = TRUE, annotationFile = "Annotation.txt", optionsFile = "Options.txt", 
	results.path = NULL, dataType = "RNA-seq", readSense = "sense",  maxReads = NULL,
	readBufferSize = 1000000)
}
\arguments{
  \item{filein}{
The name of one sorted BAM file.
}
  \item{reload}{
Ligical.  Should a new empty wiggle track data structure be created, or if FALSE add to an existing
wiggle track data structure.
}
  \item{annotationFile}{
File of sample annotation details, which specifies all needed 
sample-specific information about the samples under study.  
See \code{\link{DuffyNGS_Annotation}}.
}
  \item{optionsFile}{
File of processing options, which specifies all processing
parameters that are not sample specific.  See \code{\link{DuffyNGS_Options}}.
}
  \item{results.path}{
The top level folder path for writing result files to.  By default, read from the Options 
file entry 'results.path'.
}
  \item{dataType}{
The type of read data contained in the aligned BAM file.  
}
  \item{readSense}{
One of \code{"sense"} or \code{"antisense"}.  Controls the behavior of how alignments are assigned to
the strands.
}
  \item{maxReads}{
Optional integer.  Sets a maximum count for how many alignments will be loaded from this BAM file.  
}
  \item{readBufferSize}{
Integer.  Sets the buffer sizer for reading in alignments in efficient chunks.
}
}
\details{
This low level function is called repeatedly to combine genomic and splice BAM files.
}
\value{
A family of binary .WIG files is written to the \code{wig} subfolder.  Also the function returns an integer
giving the number of buffers of alignments that it processed, where zero denotes an error.
}
\author{
Bob Morrison
}
\seealso{
  \code{\link{pipe.AlignToWig}} {for the high level pipeline step that invokes this function.}
}
