\name{pipe.GenomicAlign}
\alias{pipe.GenomicAlign}
\title{
Genomic Alignment Pileline Step for NGS Data
}
\description{
Runs the genomic alignment step to place reads from FASTQ files onto genomic chromosomal locations.
}
\usage{
pipe.GenomicAlign(inputFastqFile, sampleID, 
	optionsFile = "Options.txt", asMatePairs = FALSE, verbose = TRUE, 
	rawReadCount = NULL)
}
\arguments{
  \item{inputFastqFile}{
character vector of one or more raw FASTQ files
}
  \item{sampleID}{
the SampleID for this sample.  This SampleID keys for a row of annotation
details in the annotation file, for getting sample-specific details.
The SampleID is also used as a sample-specific prefix for all
files created during the processing of this sample.
}
  \item{optionsFile}{
the file of processing options, which specifies all processing
parameters that are not sample specific.  See \code{\link{DuffyNGS_Options}}.
}
  \item{asMatePairs}{
logical flag, should the vector of FASTQ be treated as two mate pair files?
}
  \item{rawReadCount}{
optional argument of the number of reads in the files.  By default this gets calculated on the fly.
}
}
\details{
Aligns the given FASTQ files against the Bowtie2 genomic index specified in the options file.  Reads 
that do align are written to a BAM file in the \code{align} results subfolder.  Reads that fail to align are written
to temporary file(s) in the current working folder for the optional subsequent splice alignment step.
}
\value{
a family of BAM files, FASTQ files, and summary files, written to subfolders under the \code{results.path} 
folder.

Also, a list of alignment counts:
  \item{RawReads}{ the number of reads in the raw FASTQ files}
  \item{UniqueReads}{ the number of reads that aligned to exactly one location}
  \item{MultiReads}{ the number of reads that aligned to two or more locations}
  \item{NoHitReads}{ the number of reads that failed to align}
  \item{Time}{ the time usage of this alignment step, as from \code{\link{proc.time}}}
}
\author{
Bob Morrison
}

\seealso{
\code{\link{pipe.RiboClear}} Ribo clearing alignment against an index of unwanted genomic featurs.
\code{\link{pipe.SpliceAlign}} Splice junction alignment against an index of standard and alternative splice junctions.
}
