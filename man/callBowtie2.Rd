\name{callBowtie2}
\alias{callBowtie2}
\alias{callBowtie}
\title{ 
Call the Bowtie Alignment Tool 
}
\description{
  Call Bowtie2 (or older Bowtie) to do an alignment of a file of FASTQ read data.
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
  \item{timeIn }{ diagnostic timing information from \code{proc.time} before the call to Bowtie.}
  \item{timeOut }{ diagnostic timing information from \code{proc.time} after completing the call to Bowtie.}
}
\references{ http://bowtie-bio.sourceforge.net/ }
\author{ bob Morrison }
\seealso{ 
  	\code{\link{buildBowtie2CommandLine}}{ for constructing the command. }
  	\code{\link{bowtie2Par}}{ for default settings. }
  	\code{\link{fastqToBAM}}{ for a higher level function for doing alignments. }
}
