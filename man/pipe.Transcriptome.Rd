\name{pipe.Transcriptome}
\alias{pipe.Transcriptome}
\title{
Turn Alignments into Wiggle Tracks and Transcriptome.
}
\description{
Pipeline step that generates transcriptomes of gene expression levels for one sample.
The alignment BAM files are converted to strand specific and unique/multi-hit wiggle 
tracks, and then read pileup depth is measured for all genes in all target species.
}
\usage{
pipe.Transcriptome( sampleID, annotationFile = "Annotation.txt", optionsFile = "Options.txt", 
	speciesID = NULL, results.path = NULL, dataType = NULL, altGeneMap = NULL, i
	altGeneMapLabel = NULL, loadWIG = FALSE, verbose = TRUE, mode = "normal", exonsOnly = NULL)
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
  \item{speciesID}{
The SpeciesID of the target species to calculate a transcriptome for.  By default, transcriptome 
for all species in the current target are generated.
}
  \item{results.path}{
The top level folder path for writing result files to.  By default, read from the Options 
file entry 'results.path'.
}
  \item{dataType}{
The type of raw data contained in the FASTQ files.  By default, read from the 'DataType' field of 
the Annotation file.
}
  \item{altGeneMap}{
An alternate data frame of gene annotations, in a format identical to \code{\link{getCurrentGeneMap}},
that has the gene names and locations to be measured for read pileup depth.  By default, use the 
standard built-in gene map for each species.
}
  \item{altGeneMapLabel}{
A character string identifier for the alternate gene map, that becomes part of all created file names
to indicate the gene map that produced those transcriptomes.
}
  \item{loadWIG}{
Logical, should all wiggle track data structures be rebuilt from the alignment BAM files.  By default,
the wiggle files are only created if they do not yet exist for this sample.
}
  \item{verbose}{
Logical, send progress information to standard out.
}
  \item{mode}{
Controls the behavior of how alignments are loaded into the the wiggle track data objects.
}
  \item{exonsOnly}{
Logical, controls which regions of the gene footprint get integrated into the gene's expression
calculation.  When FALSE, the entire region of the gene's extent, without regard to exons/introns/UTRs
is used.  This is better for poorly annotated genomes like plasmodium.  When TRUE, only the regions 
inside exons are used for measuring the gene's abundance.  This mode is better for carefully annotated 
genomes with overlapping genes; but it runs slower.  When NULL, looks up the TRUE/FALSE value for this
sample from the Annotation Table argument 'ExonsOnly'.
}
}
\details{
This pipeline step turns the BAM files of alignment results into tabular files of gene expression data.
}
\value{
A family of files is created:
  \item{Wiggles}{A file and subfolder of wiggle track data objects for this sample may be
  written to the 'wig' subfolder of results.  These contain the strand and unique/multihit 
  details of all chromosomes for all species in the sample.}
  \item{Transcriptomes }{A file of gene expression for all genes in the gene annotation, for
  each species.  These are written to the 'transcript' subfolder of results.}
}
\author{
Bob Morrison
}
\seealso{
  \code{\link{pipe.Alignment}} for the previous pipeline step that turns raw FASTQ data to BAM alignments.
  \code{\link{pipe.TranscriptToHTML}} for turning the transcriptome into a web page with gene pileup images.
  \code{\link{pipe.DiffExpression}} for turning a set of transcriptomes into files of differential gene 
  expression.
}
