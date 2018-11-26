\name{DuffyNGS_Annotation}
\alias{DuffyNGS_Annotation}
\title{
Annotation File of Sample-Specific Settings 
}
\description{
The Annotation file defines the SampleID and several sample-specific settings for
each dataset to be processed.  See \code{\link{DuffyNGS_Options}} for processing 
settings that are not specific to each sample.  Each column is a named field, and
all entries are tab-delimited.  See \code{\link{AnnotationTable}} for more details.
}
\section{Annotation Fields}{
  \describe{
  \item{SampleID}{
The SampleID for this sample.  This SampleID keys for one entire row of annotation
details in the annotation file, for getting sample-specific details.
The SampleID is also used as a sample-specific prefix for all
files created during the processing of this sample.
}
  \item{Filename}{
File(s) of raw read data to be aligned; most commonly FASTQ format.  In the case of paired end
data, this field will be 2 filenames, separated by a comma with no intervening spaces.
}
  \item{PairedEnd}{
Logical.  Is this sample a paired end dataset, with 2 separate files of raw reads.
}
  \item{DataType}{
A data type for this sample.  One of:  RNA-seq,  DNA-seq,  ChIP-seq, RIP-seq
}
  \item{Group}{
A GroupID for this sample, for processing tools that combine samples by various traits.
Most such tools allow you to specify the grouping field by a \code{groupColumn} argument.
}
  \item{Color}{
A color for this sample, for processing tools that visualize multiple samples.
Most such tools allow you to specify the color field by a \code{colorColumn} argument.
}
  \item{ReadSense}{
The read orientation for this sample.  One of:  sense,  antisense
This can depend on the details of sample prep, library construction, and sequencer.
It is used if you expect to have strand specific reads, and controls which strand the
reads get assigned to.
}
  \item{StrandSpecific}{
Logical.  Did the details of sample prep, library construction, and sequencing generate
strand specific reads, and should reads only count toward expression totals, etc., if 
they land on the correct (coding) strand.
}
  \item{KeepIntergenics}{
Logical.  Should extra 'non-genes' be explicitly added to capture expression, etc., in
the intergenic spaces.  This field is queried by each tool, so it may be necessary to 
re-run some steps.  (e.g.  differential expression is based on transcription, so adding
non-genes to the DE results requires that the transcription tool was run with 
\code{KeepIntergenics = TRUE}
}
  \item{ExonsOnly}{
  Logical.  Should the entire extent of the gene's location on the chromosome be used when
  RPKM values are calculated, or should the tool restrict to only use regions inside the 
  defined exons.
}
  \item{TargetID}{
  an optional character string to define the TargetID for this sample.  If present, overrides
  the default in the options table.  This is useful when there are multiple mixed organisms
  within one experiment.
}
}}
\note{
The annotation file can contain any other columns as well.  Typically it may have more than
one column of grouping and coloring information, as samples may be compared in multiple ways.
Most sample comparison tools have explicit \code{groupColumn} and \code{colorColumn} arguments
to easily identify how samples will be grouped and colored.
}
