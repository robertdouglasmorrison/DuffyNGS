\name{pipe.AlignStats}
\alias{pipe.AlignStats}
\alias{pipe.AlignmentPie}
\title{
Generate Alignment Success Stats Images.
}
\description{
Auxiliary pipeline step that creates a family of alignment statistic images to summarize
all aspects of the alignment pipeline and its metrix.
}
\usage{
pipe.AlignStats( sampleID, annotationFile = "Annotation.txt", optionsFile = "Options.txt", 
	results.path = NULL, banner = "", chunkSize = 500000, maxReads = NULL,
	mode = c( "normal", "QuickQC"), what = NULL, plot = TRUE, fastqFile = NULL)

pipe.AlignmentPie( sampleID, annotationFile = "Annotation.txt", optionsFile = "Options.txt", 
	results.path = NULL, banner = "", mode = c( "normal", "QuickQC"), 
	fastqFile = NULL, useUSR = TRUE)
}
\arguments{
  \item{sampleID}{
The SampleID for this sample. 
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
  \item{banner}{
Optional character string to add to each plot's main heading.
}
  \item{chunkSize}{
Integer.  The buffer size to use for reading in and evaluating alignments. Most statistics are
tallied and images printed after each buffer, to show incremental progress.
}
  \item{maxReads}{
Optional integer to limit the number of alignments evaluated.
}
  \item{mode}{
Controls the behavior of how alignments are interpreted.  Mode \code{"QuickQC"} invokes the 
behavior for preliminary QC analysis.  See \code{\link{pipe.QuickQC}}.
}
  \item{what}{
An optional character string that specifies which types of statistics to monitor.  Default is to
monitor every type of feature, or "SGBIDMA" where
    \item{S}{ Sequences: features about chromosome, like read counts and percentages.}
    \item{G}{ Genes: features about genes, like read counts and percentages for highly detected genes.}
    \item{B}{ Bases: features about base calls, locations of mismatches, and nucleotide usage.}
    \item{I,D}{ Insertions & Deletions:  features about indel locations in the aligned reads.}
    \item{M}{ MARs (Multiply Aligned Reads): features about reads hitting 2+ locations.}
    \item{A}{ Align scores: features about the distribution of Bowtie alignment scores.}
}
  \item{fastqFile}{
Optional character string for the original FASTQ file that was input to the alignment pipeline.
Default is to look it up from annotation file.
}
  \item{useUSR}{
Logical.  Include a survey of USRs (Unique Short Reads) in the pie, to assess presence of empty 
adapters, Poly-N, etc.
}
}
\details{
This pipeline step tries to evaluate every aspect of how well the raw reads aligned to the 
target organism(s).  It generates a large family of plot images, each of which shows some measure of 
alignment success or failure.
}
\value{
A family of files and plot images is created on disk under the subfolder \code{AlignStats}.

Also a list of read counts and percentages as returned from the alignment pie function that
summarizes the alignment status of the entire sample.
}
\author{
Bob Morrison
}
