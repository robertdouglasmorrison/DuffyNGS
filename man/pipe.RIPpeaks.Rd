\name{pipe.RIPpeaks}
\alias{pipe.RIPpeaks}
\title{
Turn RIP-seq Data Alignments into Wiggle Tracks and RIP Peaks.
}
\description{
Pipeline step that generates files of RIP peaks for one sample.
The alignment BAM files are converted to strand specific and unique/multi-hit wiggle 
tracks, and then the read pileup depth is interogated for the presence of RIP peaks.
}
\usage{
pipe.RIPpeaks(sampleID, annotationFile = "Annotation.txt", optionsFile = "Options.txt", 
	speciesID = NULL, results.path = NULL, loadWIG = FALSE, storage.mode = "normal", 
	doPeakSearch = TRUE, peak.type = "auto", cutoff.medians = 3, 
	min.width = 50, max.width=1000,
	use.multicore = TRUE, p.value = 0.075, verbose = FALSE, 
	controlPeaksFile = NULL, scaledReadCount = NULL, visualize = interactive())
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
  \item{loadWIG}{
Logical, should all wiggle track data structures be rebuilt from the alignment BAM files.  By default,
the wiggle files are only created if they do not yet exist for this sample.
}
  \item{storage.mode}{
Controls the behavior of how alignments are loaded into the the wiggle track data objects.
}
  \item{doPeakSearch}{
Logical, controls whether the full peak search algorithm is run, or to instead just 
rerun the post-search summarization and plotting steps of this pipeline.
}
  \item{peak.type}{
Controls the type of peak shapes the picker evaluates against the raw pileup data.  Choices 
are 'gaussian', 'gumbel', 'pulse' (a singleton square wave), and 'auto'.   The default 'auto' lets the peak picker 
use all known peak shapes and select the one that best fits the raw data at each location.
}
  \item{cutoff.medians}{
Defines the cutoff threshold between noise and potential true peaks, as a multiple of the 
strand's median raw read pileup depth. Only locations having higher read depth than the 
cutoff threshold will be passed to the lower level peak picking algorithm.  Setting the 
cutoff too low will slow computation and introduce false peaks, while setting it too 
high will cause some true but smaller peaks to be missed.
}
  \item{min.width}{}
  \item{max.width}{
Used to define the expected 'half-width at half-height' (HWHH) of an expected RIP peak in the
sample.  As the RIP peaks themselves are a result of read pileups at various locations along the
genome, the best minimum estimate for the canonical width will be exactly the length in
nucleotides of the raw reads in the FASTQ data.  Since RIP peaks are a function of mRNA expression,
there is no robust way to know the upper limit on a RIP peak's width. Peak scores will get down-graded
when the observed width is outside these bounds.
}
  \item{use.multicore}{
Logical, should the two strands of raw read pileup data (forward, reverse) be 
run in parallel as two separate child peakpick processes. As of R version 3, plotting during multicore
is no longer supported, so care should be taken with this argument.
}
  \item{p.value}{
The limiting P-value for determining how many peaks to include in the final tables of RIP 
peak results.  Peaks with poor scoring final P-values are omitted.
}
  \item{verbose}{
Logical, send progress information to standard out.
}
  \item{controlPeaksFile}{
Character filename of a set of negative control peaks from samples that should not have real
RIP peaks.  This set of 'non-peaks' are used for calculating the P-values for peaks in the 
given sample.  Not yet implemented for RIP-seq data...
}
  \item{scaledReadCount}{
Total raw read count that all samples will be scaled to prior to peak picking, to assure
that all samples will be scored in an equivalent manner.  Many of the peak shape scoring
metrics are influenced by read depth, so this helps give uniform peak scoring between samples.
}
  \item{visualize}{
Logical, should the pipeline generate plot images of peak evaluation and final peak calls
during the peak search run.  Plotting can only occur in an interactive session.  Assumes
the program is running on a machine that supports X11 graphics.
}
}
\details{
This is the high level pipeline step for finding RIP peaks in a sample.  If needed, it first
turns the BAM alignments into wiggle track data; then calls the peak peaker, and finally 
generates several files of RIP peak results.

The algorithm was written and tuned for the transcription factor induction data from 
the ISB's Baliga lab for M.tuberculosis.  It may need some retuning for
other genomes and/or sequencing data variabilities.

This function operates on a single sample.  In a typical experiment, there will be multiple 
samples that need to be evaluated as a set.  See \code{\link{pipe.RIPpeaksToAltGeneMap}} for
combining RIP peaks from multiple samples into one set of features to be interogated for differential
expression abundance.
}
\value{
A family of files are written to this sample's subfolder in the 'RIPpeaks' folder of results.
These include:
  \item{Strand specific peak files}{ Two text files of single strand peak calls, for the 
  forward and reverse wiggle tracks.}
  \item{Strand specific ROC curve plots}{ If interactive graphics were generated, ROC curves 
  of the peak pick score distributions are made.}
  \item{Final peaks text file}{ One text file of final peak calls, after merging the 2 strand-specific
  results.  There is no joining based on forward/reverse peak pairs, as that behavior is not 
  required or expected in RIP-seq data.}
  \item{Final peaks .CSV file}{ One Excel readable file of final peak calls, after 
  appending a column of the genomic DNA sequence under the peak.  This column is
  configured as FASTA sequences that can be submitted to a motif detection algorithm to
  find the binding motif of that sample's transcription factor.}
  \item{Final peaks .html file}{ One web browser readable file of final peak calls, with
  hyperlinks into a subfolder of plot images of all high scoring final RIP peaks.}
  \item{Peaks scoring details}{ One text file of all scoring metric details for all identified
  peaks.  Useful for tuning the algorithm or just understanding why peaks got scored as
  they did.}
}
\references{
  RIPSeeker: a statistical package for identifying protein-associated transcripts from RIP-seq experiments.
  Yue Li, Dorothy Yanling Zhao, Jack Greenblatt, Zhaolie Zhang
  Nucleic Acids Research, 2013, Vol. 41, No 8.
}
\author{
Bob Morrison
}
\seealso{
  \code{\link{calcWigRIPpeaks}} for the intermediate level implementation, and
  \code{\link{findRIPpeaks}} for the low level single strand read pileup depth peak pick tool.
  \code{\link{pipe.RIPpeaksToAltGeneMap}} for the next step of joining peak sites into gene features for
  differential detection.
}
