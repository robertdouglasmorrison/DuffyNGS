\name{pipe.BAMpileup}
\alias{pipe.BAMpileup}
\alias{BAM.mpileup}
\title{
Extract Pileup Details from BAM file.
}
\description{
High level and low level wrappper functions to invoke SAMTOOLS MPILEUP to get alignment details
for a region in a BAM file.
}
\usage{
pipe.BAMpileup(sampleID, geneID = NULL, seqID = NULL, start = NULL, stop = NULL, 
	annotationFile = "Annotation.txt", optionsFile = "Options.txt", 
	results.path = NULL, summarize.calls = TRUE, verbose = FALSE)

BAM.mpileup(files, seqID, fastaFile, start = NULL, stop = NULL, 
	min.depth = 3, max.depth = 10000, min.gap.fraction = 0.25, 
	mpileupArgs = "", summarize.calls = FALSE, verbose = TRUE)
}
\arguments{
  \item{sampleID}{
the SampleID for this sample.  This SampleID keys for a row of annotation
details in the annotation file, for getting sample-specific details.
The SampleID is also used as a sample-specific prefix for all
files created during the processing of this sample.
}
  \item{geneID}{
Optional character string of one GeneID.
}
  \item{seqID}{
Optional character string of one SeqID.  Required explicitily for the low level function.
}
  \item{start}{
Optional integer of a starting nucleotide.
}
  \item{stop}{
Optional integer of a stopping nucleotide.
Note that any valid combination of the above that is sufficient to define a region of a chromosome is allowed.
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
The top level folder where all results have been written to. Default is taken from the options file.
}
  \item{summarize.calls}{
Logical.  Controls whether the very low level details of the base calls at each nucleotide are summarized
into tabular form.  FALSE leaves the raw large cryptic strings as made by SAMTOOLS.  When TRUE, they are 
converted into a short text string count summary of A/C/G/T/Indels.
}
  \item{fastaFile}{
Full pathname to the genomic DNA FASTA file that made the Bowtie index that the BAM file was aligned against.
}
  \item{min.depth}{}
  \item{max.depth}{}
  \item{min.gap.fraction}{
Bounds passed to SAMTOOLS MPILEUP, to control how it performs and reports its assessment of the BAM file.
}
  \item{mpileupArgs}{
Other optional arguments to pass to SAMTOOLS MPILEUP.
}
}
\details{
These functions give two different levels of control for extracting aligned read pileup information from
BAM files, via the SAMTOOLS MPILEUP utility.  The high level function gives an easy way to extract details for
single genes or regions, while the low level function gives finer control.
}
\value{
A data frame of details about read coverage and base depth in the specified region of a chromosome.
Only base locations with 1+ read of coverage are returned, so there can be gaps with no information.
Columns include:
  \item{SEQ_ID}{ chromosome name as a character string}
  \item{POSITION}{ chromosome location as an integer}
  \item{REF_BASE}{ the expected base call of the reference genome}
  \item{DEPTH}{ the depth of coverage at this position, as an integer}
  \item{CALL_BASE}{ when summarize is TRUE, the one base that was most frequently observed.
  If FALSE, the the full cryptic string of base calls as made by MPILEUP}
  \item{CALL_SCORE}{ only when summarize is FALSE, the full cryptic string of Phred quality scores}
  \item{BASE_TABLE}{ only when summarize is TRUE, the short sorted table of observed base counts}
}
\references{
SAMTOOLS  www.htslib.org/doc/samtools.html
}
\author{
Bob Morrison
}
\note{
Note that the high level function only operates on a single sample.  While the low level function can
be given a vector of BAM files, and then returns depth and call details for multiple files at one time.
}
