\name{pipe.RiboClear}
\alias{pipe.RiboClear}
\title{
Ribo Clearing Pileline Step for NGS Data
}
\description{
Runs the ribo clearing alignment step to remove unwanted genomic feature reads from FASTQ files,
prior to the main genomic alignment and splice alignment steps.  Typically invoked on RNA-seq data
to remove highly abundant unwanted transcripts such as ribosomal RNA, mitochondrial RNA, albumin,
and very high abundance small RNAs; features whose expression abundance might dwarf the expression
levels of typical wanted genes.
}
\usage{
pipe.RiboClear(inputFastqFile, sampleID, annotationFile = "Annotation.txt", 
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
  \item{annotationFile}{
the file of sample annotation details, which specifies all needed 
sample-specific information about the samples under study.  
See \code{\link{DuffyNGS_Annotation}}.
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
Aligns the given FASTQ files against the Bowtie2 ribo clearing index specified in the options file.  Reads 
that do align are written to a BAM file in the \code{riboClear} results subfolder.  Reads that fail to align are written
to temporary file(s) in the current working folder for the subsequent genomic alignment step.

The genomic features to be removed are hard-coded into the \code{RrnaMap} map and the Bowtie ribo clearing 
indexes.  It is not currently something that can be modified by the end user.  See \code{\link{MapSets}}.
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
\note{
the interplay of paired end strand specific reads and ribo clearing is messy.  Typically, ribo clearing finds
both multi-hit reads and un-paired alignments where only one read maps.  Both of these behavoirs break the 
expected convention that paired end read files are still perfect mate pairs after the alignment step.  The
default mode is to turn off paired end strand specific behavior when ribo clearing, but this can be forced 
using option \code{forcePairedEnd}.
}

\seealso{
\code{\link{pipe.GenomicAlign}} {Genomic alignment against an index of target genome(s).}
\code{\link{pipe.SpliceAlign}} {Splice junction alignment against an index of standard and alternative splice junctions.}
}
