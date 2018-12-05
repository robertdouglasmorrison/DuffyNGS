\name{pipe.PlotSNP}
\alias{pipe.PlotSNP}
\title{
Plot Single Nucleotide Polymorphism (SNP) Sites at Nucleotide Resolution.
%%  ~~function to do ... ~~
}
\description{
Plot aligned read base call distributions for one or more SNP sites from one or more samples,
and optionally save them as PNG images.
}
\usage{
pipe.PlotSNP(sampleIDs, seqID, position, geneID = NULL, groups = sampleIDs, 
		annotationFile = "Annotation.txt", optionsFile = "Options.txt", results.path = NULL, 
		tailWidth = NULL, plotFormat = c("", "png", "pdf"), plotFileName = NULL, 
		plot.path = "SNP_Plots", label = "", SNPtablePath = "~/SNPs/", mf = NULL, 
		start = NULL, stop = NULL, verbose = TRUE, ...)
}
\arguments{
  \item{sampleIDs}{
A vector or one or more SampleIDs.
}
  \item{seqID}{
Character string of one chromosome name.
}
  \item{position}{
Integer vector of one or more chromosomal locations, in nucleotides.
}
  \item{geneID}{
Optional, give the GeneID of the one SNP site (just to save time).
}
  \item{groups}{
Character vector of group identifiers, of the same length as \code{sampleIDs}, used in the
plot headers to show sample grouping information when making multi-sample plots.
}
  \item{annotationFile}{
The file of sample annotations.
}
  \item{optionsFile}{
The file of program options.
}
  \item{results.path}{
The folder of pipeline results, to find existing BAM files for extracting read pileup details.
}
  \item{tailWidth}{
Integer number of extra nucleotides of context, to include in the SNP plot image on both sides
of the one center position of interest.  There is a tradeoff between showing so much context that no
SNPs are discernable, and showing so little context that the overall alignment quality in the
area cannot be assured.  Further, the plot function lays down so many extra supporting features
that it is very hard to keep all content clear and distinct.  Adjusting the tail width can
help find the optimal tradeoff between clear details and adequate context.
}
  \item{plotFormat}{
Controls the plotting behaviour.  Default is to write to the current graphics device. Optionally 
it can create PNG or PDF files on disk, using the name of the gene and nucleotide position as the 
default filename.
}
  \item{plotFileName}{
Optional character string of the filename to use for the created plot image file.
}
  \item{plot.path}{
Folder to write PNG or PDF image files to.
}
  \item{label}{
Optional vector of text, same length as \code{sampleIDs}, for adding to plot headers.
}
  \item{SNPtablePath}{
Pathname to optional tables of known SNP sites for this chromosome, for annotating the lower margin
of the SNP plot with previously identified known SNP sites. (this functionality is only available 
for some organisms, and how it's implemented varies with the size of the genome...).
}
  \item{mf}{
An optional vector of length 2, as used by graphics parameters \code{mfcol, mfrow}, to control the
layout of composite SNP plots that show 2 or more plots.  A default layout is chosen based on the 
number of plots to be made, but finer control can give a more informative plot layout by careful ordering
of the sample IDs and explicitly setting 'mf'.
}
  \item{start}{}
  \item{stop}{
Optional explicit boundaries, in nucleotide units, of the region around the center SNP site to plot.  
Overrides \code{tailWidth}.
}
  \item{\dots}{
Other arguments to be passed to the low level SNP plot function.
}
}
\details{
This function gives a simpler interface wrapper to the lower level SNP plot function, 
handling all details about SpeciesIDs, BAM files, etc.; as well as accomodating getting multiple 
SNP plot images into one single composite plot, as for comparing between samples.
}
\value{
If \code{plotFormat} is not empty, then PNG or PDF files are written to disk. Otherwise creates a SNP 
plot image in the current graphics device.
}
\note{
Showing SNP read pileup images for many samples at one time is quite slow, as it consumes lots of 
memory and time reading from multible BAM files.  When the number of samples gets large, the available
resolution per sample gets so small that making usable images becomes impossible.
}
\seealso{
  \code{\link{pipe.PlotGene}} for aligned read plot images that show the depth of read coverage over
  entire genes, in one or more samples.
}
\author{
Bob Morrison
}
