\name{getBowtie2AlignMetrics}
\alias{getBowtie2AlignMetrics}
\title{
Get Bowtie2 Alignment Metrics
}
\description{
Gets the reads counts for all the various alignment success categories from the
(stderr) output of Bowtie2.
}
\usage{
getBowtie2AlignMetrics(metricsFile, bowtieTiming = NULL, verbose = TRUE)
}
\arguments{
  \item{metricsFile}{
the filename to receive all Bowtie2 output that got send to \code{stderr}.
}
  \item{bowtieTiming}{
a list of 2 \code{\link{proc.time}} objects that are the start and stop time
of the Bowtie alignment step.
}
}
\details{
Since the call to Bowtie2 is a pipe sending its output to samtools, the final
summary counts that Bowtie2 reports get saved to a text file that captures \code{stderr}.  This function
extracts those counts and returns them to the R session that called Bowtie2.
}
\value{
A list of numeric read counts (not the number of alignments).
  \item{RawReads }{ the count of reads in the input FASTQ file. }
  \item{UniqueReads }{ the count of reads that aligned exactly once. }
  \item{MultiReads }{ the count of reads that aligned >1 times. }
  \item{NoHitReads }{ the count of reads that failed to aligned. }
  \item{Time }{ the elapsed time of the call to Bowtie, as a 
  \code{\link{elapsedProcTime}} object }
}
\note{
Dependant on the Bowtie2 version for successful parsing of Bowtie's count output. 
No rigorous testing for any runtime errors is done.
}

