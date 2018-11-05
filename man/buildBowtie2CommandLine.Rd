\name{buildBowtie2CommandLine}
\alias{buildBowtie2CommandLine}
\title{ 
Build a Bowtie2 Command Line 
}
\description{
  Construct one bowtie2 command line from a set of file names, options, and Bowtie
  index information.
}
\usage{
buildBowtie2CommandLine(inputFastqFile, 
	outputAlignFile = sub("fastq$", "align", inputFastqFile), 
	m = 1, k = NULL, noHitsFile = NULL, multiHitsFile = NULL, 
	optionsFile = "Options.txt",
	alignIndex=getOptionValue( optionsFile, "GenomicIndex"),
	alignPolicy=getOptionValue( optionsFile, "GenomicAlignmentPolicy"),
	maxReads=NULL, verbose=TRUE, debug=FALSE) 
}
\arguments{
  \item{inputFastqFile}{ required file name of .fastq reads to be aligned. }
  \item{outputAlignFile}{ name for the to-be-created file of alignments. }
  \item{m}{ maximum number of alingment locations allowed to call a read 'aligned'.
  This is the value passed to the '-m' option in Bowtie. }
  \item{k}{ maximum number of alingment locations to report for an aligned read.
  This is the value passed to the '-k' option in Bowtie. }
  \item{noHitsFile}{ optional file name to collect reads that do not align.}
  \item{multiHitsFile}{ optional file name to collect reads that align to more
  than 'm' locations.}
  \item{alignIndex}{ the basename of the Bowtie index to align against. Note the the path to all
  index files is set via \code{\link{bowtiePar('IndexPath')}}. }
  \item{alignPolicy}{ a character string of additional Bowtie arguments that include any
  alignment policy and threshold setting for calling a read 'aligned'.}
  \item{maxReads}{the maximum number of reads to pass to Bowtie}
}
\details{
  This takes the current 'bowtie2Par()' settings and the given arguments, to
  construct one Unix style command line that can be used to run Bowtie2 to
  generate alignments.

  If \code{inputFastqFile} is compressed, decompression is done via a Unix pipe.
  When \code{maxReads} is not \code{NULL}, the subset of reads is used via 
  Unix 'head'.
}
\value{
  A character string of length 1 that can be passed to \code{callBowtie2} or to the 'system()' 
  command to invoke Bowtie2.
}
\references{ http://bowtie-bio.sourceforge.net/ }
\author{ Bob Morrison }
\note{ 
  Various Bowtie2 command options from the \code{optionsFile} are prepended before the alignment
  policy parameters, so the policy can override as appropriate.
}
\seealso{  \code{\link{callBowtie2}},  \code{\link{bowtie2Par}},  \code{\link{fastqToBAM}}  }
