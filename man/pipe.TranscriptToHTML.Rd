\name{pipe.TranscriptToHTML}
\alias{pipe.TranscriptToHTML}
\title{
Turn transcriptome into HTML table with gene plot hyperlinks.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pipe.TranscriptToHTML(sampleID, speciesID = NULL, annotationFile = "Annotation.txt", 
	optionsFile = "Options.txt", results.path = NULL, N = 100, minYmax = 10, 
	tailWidth = 1000, pause = NULL, altGeneMap = NULL, altGeneMapLabel = NULL, 
	targetID = NULL, plotType = "box", useLog = FALSE, label = "", 
	PLOT.FUN = NULL, geneSet = NULL, ...)
}
\arguments{
  \item{sampleID}{
the SampleID for this sample.  This SampleID keys for a row of annotation
details in the annotation file, for getting sample-specific details.
The SampleID is also used as a sample-specific prefix for all
files created during the processing of this sample.
}
  \item{speciesID}{
The SpeciesID of the target species to operate on.  By default, all species in the current 
target are turned into HTML tables.
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
  \item{results.path}{
The top level folder path for writing result files to.  By default, read from the Options 
file entry 'results.path'.
}
  \item{N}{
The number of genes to include in the HTML table, and make gene plot images for.
}
  \item{minYmax}{
Plotting option passed to low level default gene plot function, sets the lowest maximum Y axis limit.
Useful for making low and high expressing genes share similar plot bounds.
}
  \item{tailWidth}{
Plotting option passed to low level default gene plot function, sets the extent of intergenic space 
outside the gene boundaries to include in the plot, to enlarge the X axis limits.
Useful for seeing the genomic neighborhood around genes.  Units are in nucleotides.
}
  \item{plotType}{
Plotting option passed to low level default gene plot function, controls the type of read pileup
images generated.  May be one of "boxes" for color-filled boxes, "lines" for colored lines, and
"segments" for colored line segments.
}
  \item{useLog}{
logical, should the Y axis be shown as log scaled.
}
  \item{altGeneMap}{
An alternate data frame of gene annotations, in a format identical to \code{\link{getCurrentGeneMap}},
that has the gene names and locations to be measured for read pileup depth.  By default, use the 
standard built-in gene map for each species.
}
  \item{altGeneMapLabel}{
A character string identifier for the alternate gene map, that becomes part of all created path and file names
to indicate the gene map that produced those transcriptomes.
}
  \item{label}{
optional text to append as part of the "main" plot legend.
}
  \item{PLOT.FUN}{
the name of an optional function to call instead of the default read pileups function.  Only requirement
is the that the first argument be the name of the gene to be drawn.
}
  \item{geneSet}{
an optional character vector of GeneIDs, that limits the HTML table to only include genes in this set.
}
  \item{\dots}{
other arguments passed down to the low level gene plotting function.
}
}
\details{
This is an optional pipeline tool to convert static transcriptome results into more interactive HTML.  
Writes files and gene plot images to the \code{html} subfolder under the main results path.  Including 
too many genes makes for unresponsive HTML and overly large plot folders, so use discretion with "N".
}
\value{
Writes HTML files and hyperlinked gene plot images to the \code{htlm} subfolder.
}
\author{
Bob Morrison
}
\seealso{
  \code{\link{pipe.Transcriptome}} { for creating transcriptome files from BAM files.}
  \code{\link{pipe.PlotGene}} { for making default gene plot images.}
  \code{\link{table2html}} { for the low level tool that turns data frames to HTML.}
}
