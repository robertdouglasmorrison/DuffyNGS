\name{fastqToBAM}
\alias{fastqToBAM}
\title{ 
Convert FASTQ file(s) into a BAM file, by calling Bowtie2.
}
\description{
  Call Bowtie2 to do an alignment of a file of FASTQ read data.  Builds the 
  full Unix command line needed to spawn a call to Bowtie, using all pertinent options file settings.
}
\usage{
fastqToBAM(inputFastqFile, outputFile = sub("(fastq|fq|gz)$", "bam", inputFastqFile[1]), k = NULL, 
	sampleID = "", optionsFile = "Options.txt", noHitsFile = NULL, 
	alignIndex = getOptionValue(optionsFile, "GenomicIndex"), index.path = NULL, 
	alignPolicy = getOptionValue(optionsFile, "GenomicAlignmentPolicy", verbose = F), 
	maxReads = NULL, skipReads = NULL, asMatePairs = FALSE, keepUnaligned = TRUE, 
	wait = TRUE, verbose = FALSE, label = sampleID)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputFastqFile}{
character vector of explicit full path qualified FASTQ filenames.
}
  \item{outputFile}{
character string of the full path qualified name to give to the result BAM file that Bowtie will create.
}
  \item{k}{
the value for Bowtie argument '-k' that finds multiple alignments.
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
  \item{noHitsFile}{
character string of the full path qualified name for the created FASTQ file of reads that fail to align.
Default value is to not save the reads that didn't align.
}
  \item{alignIndex}{
character string of the base filename for the Bowtie index to align against.
}
  \item{index.path}{
character string of the full pathname to a directory of Bowtie indexes, where the above index is located.
}
  \item{alignPolicy}{
character string passed to Bowtie that gives the options to be used for alignment.
}
  \item{maxReads}{
maximum number of reads to have Bowtie align.  See Bowtie option '--qupto'.
}
  \item{skipReads}{
number of reads to ignore at the beginning of the FASTQ file.  See Bowtie option '--skip'.
}
  \item{asMatePairs}{
logical, controls how the vector of FASTQ files are interpreted, and selects paired end or unpaired Bowtie 
alignment mode.
}
  \item{keepUnaligned}{
logical, should unaligned reads be written to the results BAM file.  See Bowtie option '--un'.
}
}
\details{
This is the mid-level wrapper function that pulls the various DuffyNGS alignment settings from 
the options file, and constructs the correct Unix command line to invoke an explicit Bowtie 
alignment step.  This wrapper is used inside all 3 top level alignment steps of ribo clearing, 
genomic, and splice junctions.
}
\value{
  Primarily, one new BAM file of alignments.
  
  Also a list of diagnostic information about the alignment success:
  \item{RawReads}{ the number of raw reads in the FASTQ file submitted to Bowtie. }
  \item{UniqueReads}{ the number of reads that aligned to exactly one location}
  \item{MultiReads}{ the number of reads that aligned to two or more locations}
  \item{NoHitReads}{ the number of reads that failed to align}
  \item{Time}{ the time usage of this Bowtie alignment step, as from \code{\link{proc.time}}}
}
\author{
Bob Morrison
}

\seealso{ 
  \code{\link{pipe.RiboClear}} {Ribo clearing against an index of unwanted transcripts.}
  \code{\link{pipe.GenomicAlign}} {Genomic alignment against an index of target genome(s).}
  \code{\link{pipe.SpliceAlign}} {Splice junction alignment against an index of standard and alternative 
  splice junctions.}
  \code{\link{buildBowtie2CommandLine}} {for constructing the command.}
  \code{\link{bowtie2Par}} {about default Bowtie settings.}
  \code{\link{callBowtie2}} {for the low-level call to Bowtie to do the alignments.}
}
