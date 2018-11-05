\name{pipe.PostAlignTasks}
\alias{pipe.PostAlignTasks}
\title{
Post-Alignment Tasks
}
\description{
Do secondary processing steps after running the alignment pipeline for a sample.
}
\usage{
pipe.PostAlignTasks(sampleID, annotationFile = "Annotation.txt", optionsFile = "Options.txt", 
	dataType = c( "RNA-seq", "ChIP-seq", "DNA-seq"), verbose = TRUE)
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
  \item{dataType}{
the type of NGS data to be aligned, determines the folders/files that will
be generated.
}
}
\details{
Runs any optional secondary processing steps after the main alignment pipeline is
finished.   Currently, the only optional step is 'NoHits-Analysis'.
}
