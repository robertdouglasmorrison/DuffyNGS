\name{pipe.TranscriptCluster}
\alias{pipe.TranscriptCluster}
\title{
Hierarchical Clustering of Transcriptomes.
}
\description{
Combine a set of sample transcriptomes into a matrix and pass it to the hierarchical clustering tool.
}
\usage{
pipe.TranscriptCluster(sampleIDset = NULL, annotationFile = "Annotation.txt", 
		optionsFile = "Options.txt", results.path = NULL, speciesID = getCurrentSpecies(), 
		intensityColumn = "RPKM_M", useLog = FALSE, label = "")
}
\arguments{
  \item{sampleIDset}{
The set of SampleIDs to compare by clustering.
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
  \item{intensityColumn}{
Character string of one column name in the standard transcriptome files.
}
  \item{useLog}{
Logical, should the intensity values be log transformed prior to clustering.
}
  \item{label}{
Character string passed to the plot routine, added to the main header.
}
}
\details{
A simple wrapper function that first builds a matrix of gene expression, by calling
\code{\link{expressionFileSetToMatrix}}, and then passing that on the the clustering
tool \code{\link{expressionCluster}}.
}
\value{
Invisibly, a list with 2 elements:
  \item{expressionMatrix}{a data frame of gene expression, with columns for GeneID, Product, and then a
  column of intensity values for each sample.}
  \item{cluster}{the results object passed back from the clustering tool.}
}
\author{
Bob Morrison
}
