\name{pipe.PlotGene}
\alias{pipe.PlotGene}
\title{
Plot Gene Expression as Wiggle Track Pileups
}
\description{
Plot gene expression for one or more genes from one or more samples,
and optionally save them as PNG images.
}
\usage{
pipe.PlotGene(sampleIDs, genes, annotationFile = "Annotation.txt", 
	optionsFile = "Options.txt", results.path = NULL, PLOT.FUN = NULL, 
	asPNG = FALSE, path = ".", keepShortGeneName = NULL, ...)
}
\arguments{
  \item{sampleIDs}{
A vector or one or more SampleIDs.
}
  \item{genes}{
A vector of one or more genes, from any species in the current target.
}
  \item{annotationFile}{
The file of sample annotations.
}
  \item{optionsFile}{
The file of program options.
}
  \item{results.path}{
The folder of pipeline results, to find existing wiggle track objects for each SampleID.
}
  \item{PLOT.FUN}{
Optional function to call for generating the gene expression plot.  By default, call either
\code{\link{plotWIGgene}} or \code{\link{plotMultiWIGgene}} depending on the number of SampleIDs.
}
  \item{asPNG}{
Logical, should PNG image files be made for each gene.
}
  \item{path}{
Folder to write .PNG image files to.
}
  \item{keepShortGeneName}{
When not NULL, an integer value passed to \code{\link{shortGeneName}} to trim the PNG image 
file name down to a shortened version of the gene name.
}
  \item{\dots}{
Other arguments to be passed to the plot function, including 'type', 'useLog', 'tailWidth', etc.
}
}
\details{
This function gives a simpler interface wrapper to the lower level gene expression plot functions, 
handling all details about SpeciesIDs, WIG objects, etc.
}
\value{
If \code{asPNG}, then .PNG files are written to disk.
}
