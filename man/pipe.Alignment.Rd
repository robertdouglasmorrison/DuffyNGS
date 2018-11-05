\name{pipe.RNAalignment}
\alias{pipe.RNAalignment}
\title{ Alignment Stage of DuffyNGS Pipeline for RNA-seq Data }
\description{
Turn a FASTQ file(s) of raw reads into RiboCleared, Genomic, and SpliceJunction alignments.
Encapsulates a 3 stage Bowtie alignment protocol for RNA-seq data, to turn raw reads
into alignments.  First, a RiboClearing stage to remove ribosomal RNA and other
high expression unwanted transcripts. Next a genomic alignment stage to
capture the reads of interest matching the target genome(s). Lastly, a splice alignment
stage to capture reads spanning known splice junctions.  RiboClearing and splice
junction stages are optional steps.
}
\usage{
pipe.RNAalignment(sampleID = NULL, annotationFile = "Annotation.txt", 
	optionsFile = "Options.txt", verbose = TRUE)
}
\arguments{
  \item{sampleID}{
the SampleID for this sample.  This SampleID keys for a row of annotation
details in the annotation file, for getting sample-specific details.
The SampleID is also used as a sample-specific prefix for all
files created during the processing of this sample.
}
  \item{annotationFile}{
the file of sample annotation details, which specifies all needed 
sample-specific information about the samples under study.  
See \code{\link{DuffyNGS_Annotation}}.
}
  \item{optionsFile}{
the file of processing options, which specifies all processing
parameters that are not sample specific.  See \code{\link{DuffyNGS_Options}}.
}
}
\details{
Turn FASTQ file(s) of raw reads into 3 groups of alignments and 2 groups of
unaligned reads.  Details affecting alignment indexes and Bowtie2 parameters
are set in the \code{\link{DuffyNGS_Options}} file.
Main components:

\code{\link{pipe.RiboClear}} Ribo clearing against an index of unwanted
transcripts.

\code{\link{pipe.GenomicAlign}} Genomic alignment against an index of
target genome(s).

\code{\link{pipe.SpliceAlign}} Splice junction alignment against an index of
standard and alternative splice junctions.

This function is called as the alignment stage of the top level \code{\link{pipeline}} tool.
}
\value{
a family of BAM files, FASTQ files, and summary files, written to subfolders
under the top level \code{results.path} folder.

Also, a list of alignment counts:
  \item{nReadsIn}{ the number of reads in the raw FASTQ file}
  \item{nNoHit}{ the number of reads that failed to align to any stage of the alignment pipeline}
  \item{nRibo}{ the number of reads that aligned to the RiboClearing stage}
  \item{nGenomic}{ the number of reads that aligned to the Genomic stage of the alignment pipeline}
  \item{nSplice}{ the number of reads that aligned to the SpliceJunction stage of the alignment pipeline}

Note that in the case of paired end reads, the alignment counts returned are just for the second pairs.
}
\author{
Bob Morrison
}
