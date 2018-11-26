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
	optionsFile = "Options.txt", results.path = NULL, targetID = NULL, colorColumn = "Color",
	plotFormat = c( "", "png", "pdf"), plot.path = "Gene_Plots", keepShortGeneName = NULL, 
	pileup.col=c(4,2,1), png.width=1000, png.height=700, pdf.width=10, pdf.height=7, ...)
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
  \item{targetID}{
optionally specify the current target organism(s).  Affects the display of genome uniqueness and
detectability on default read pileup plots.
}
  \item{PLOT.FUN}{
Optional function to call for generating the gene expression plot.  By default, call either
\code{\link{plotWIGgene}} or \code{\link{plotMultiWIGgene}} depending on the number of SampleIDs.
}
  \item{plotFormat}{
Controls the plotting behaviour.  Default is to write to the current graphics device. Optionally 
it can create PNG or PDF files on disk, using the name of the gene as the filename.
}
  \item{plot.path}{
Folder to write PNG or PDF image files to.
}
  \item{keepShortGeneName}{
When not NULL, an integer value passed to \code{\link{shortGeneName}} to trim the image 
filename down to a shortened version of the gene name.
}
  \item{pileup.col}{
a vector of 3 colors, for drawing of forward strand, reverse strand, and "combo strand" colors, respectively. 
(Combo strand only applies to ChIP-seq plots)
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
If \code{plotFormat} is not empty, then PNG or PDF files are written to disk. Otherwise creates a gene 
expression plot image in the current graphics device.
}
\note{
The default behavoir of showing full read pileup images for many samples at one time consumes lots of 
memory and time.  When the number of samples gets large, consider providing a simpler plot function 
based on some other display metric such as RPKM expression values.  Secondly, the pileup images are 
not ideal for quantitative comparisons. 
}
