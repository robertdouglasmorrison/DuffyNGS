\name{pipe.PreAlignTasks}
\alias{pipe.PreAlignTasks}
\title{
Pre-Alignment Tasks
}
\description{
Do setup prior to running the alignment pipeline for a sample.
}
\usage{
pipe.PreAlignTasks(sampleID, annotationFile = "Annotation.txt", optionsFile = "Options.txt", 
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
Creates sub-folders for all results, and pre-deletes any old previous results files 
for that SampleID prior to alignment.
}
