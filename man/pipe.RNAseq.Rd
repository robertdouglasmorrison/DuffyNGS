\name{pipe.RNAseq}
\alias{pipe.RNAseq}
\alias{pipe.ChIPseq}
\alias{pipe.DNAseq}
\alias{pipeline}
\title{
DuffyNGS Pipelines
}
\description{
The top level all-in-one function for turning raw RNA-seq FASTQ files into
wiggle tracks, transcriptomes, etc.
}
\usage{
pipe.RNAseq(  sampleID = NULL, annotationFile = "Annotation.txt", optionsFile = "Options.txt")
pipe.ChIPseq( sampleID = NULL, annotationFile = "Annotation.txt", optionsFile = "Options.txt")
pipe.DNAseq(  sampleID = NULL, annotationFile = "Annotation.txt", optionsFile = "Options.txt")

pipeline( sampleID = NULL, annotationFile = "Annotation.txt", optionsFile = "Options.txt")
}
\arguments{
  \item{sampleID}{
The SampleID for this sample.  This SampleID keys for one row of annotation
details in the annotation file, for getting sample-specific details.
The SampleID is also used as a sample-specific prefix for all
files created during the processing of this sample.
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
}
\details{
Turn raw sequencing reads into final results for that type of data, building various files of results along
the way.  \code{pipeline} is a wrapper function that looks at the annotation field \code{DataType} to dispatch
the correct pipeline tool.  Main components:

\code{\link{pipe.PreAlignTasks}} Setup.  Create the results directories (and 
cleanup of any previous files) prior to running the pipeline on this sample.

\code{\link{pipe.Alignment}} Alignment.  Turn FASTQ data into the various types of
alignments via the Bowtie2 alignment program.

\code{\link{pipe.Transcript}} Transcription (if RNA-seq data) or ChIP Peaks (if ChIP-seq data).  
Turn alignments into wiggle tracks, and then generate transcriptome(s) or ChIP peak calls for all target species.

\code{\link{pipe.PostAlignTasks}} Optional post-alignment tasks, such as de novo assembly of 
No-Hit reads.

For a typical dataset, this pipeline can take many hours and is usually configured as a batch job
submitted to a compute cluster.
}
\value{
a family of data files written to disk, and a text log of all details sent to standard out.  
}
\author{
Bob Morrison
}
