\name{pipe.AlignToWig}
\alias{pipe.AlignToWig}
\title{
Turn BAM Alignments into Wiggle Tracks.
}
\description{
Pipeline step that generates the wiggle tracks of read pileups for one sample.
The alignment BAM files are converted to R data objects that contain strand specific 
and unique/multi-hit wiggle track details to nucleotide resolution for all species and 
chromosomes in the current target.
}
\usage{
pipe.AlignToWig( sampleID, annotationFile = "Annotation.txt", optionsFile = "Options.txt", 
	results.path = NULL, dataType = "RNA-seq", 
	mode = c( "normal", "QuickQC", "spliceOnly"))
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
  \item{dataType}{
The type of read data contained in the aligned BAM file.  
}
  \item{mode}{
Controls the behavior of how alignments are loaded into the the wiggle track data objects.  Optional 
mode \code{"QuickQC"} invokes the behavior for preliminary QC analysis.  See \code{\link{pipe.QuickQC}}.
}
}
\details{
This pipeline step turns the BAM files of alignment results into wiggle track data objects.  It combines
the genomic alignments BAM file with the optional splice junction alignments BAM file, to create the 
nucleotide resolution read depth information needed for later transcriptomes and read pileup image 
creation.
}
\value{
A family of files is created on disk:
  \item{Wiggles}{A file and subfolder of wiggle track data objects for this sample is
  written to the 'wig' subfolder of results.  These contain the strand and unique/multihit 
  details of all chromosomes for all species in the sample.}
}
\author{
Bob Morrison
}
\seealso{
  \code{\link{alignToWig}} {for the low level function that does the real work.}
}
