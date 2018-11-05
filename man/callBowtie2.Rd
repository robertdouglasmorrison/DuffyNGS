\name{callBowtie2}
\alias{callBowtie2}
\title{ 
Call the Bowtie2 Alignment Tool 
}
\description{
  Call Bowtie2 to do an alignment of a file of FASTQ read data.
}
\usage{
callBowtie2(bowtieCommandLine, wait = TRUE, verbose = FALSE)
}
\arguments{
  \item{bowtieCommandLine}{
  the Unix command line that will run Bowtie2.  
  	See \code{\link{buildBowtie2CommandLine}}
  }
  \item{wait}{ 
  should the calling function wait for Bowtie2 to finish 
  }
  \item{verbose}{ 
  prints diagnostic and timing information
  }
}
\details{
  Dispatches the call to Bowtie2 via the R 'system()' function.
}
\value{
  A list of diagnostic information:
  \item{RawReads}{ the number of raw reads in the FASTQ file submitted to Bowtie. }
  \item{Alignments}{ the number of alignments in the main output file written by Bowtie. }
  \item{ReadsAligned}{ the number of reads aligned by Bowtie. This will be less than
  'Alignments' whenever a read reports more that one alignment. }
  \item{NoHits}{ the number of raw reads that failed to align.}
  \item{TooManyHits}{ the number of raw reads not aligned because they matched too
  many locations to satisfy the -m and/or -k Bowtie parameters.}
  \item{Time }{ diagnostic timing information.}
}
\references{ http://bowtie-bio.sourceforge.net/ }
\author{ bob Morrison }
\seealso{ 
  	\code{\link{buildBowtie2CommandLine}}{ for constructing the command. }
  	\code{\link{bowtie2Par}}{ for default settings. }
  	\code{\link{fastqToBAM}}{ for a higher level function for
	doing alignments. }
}
