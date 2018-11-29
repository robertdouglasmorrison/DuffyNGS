\name{pipe.RIPpeaksToAltGeneMap}
\alias{pipe.RIPpeaksToAltGeneMap}
\title{
Combine RIP Peaks into an Alternate Gene Map
}
\description{
Merge RIP peaks from multiple samples into one set of gene-like features
that can be quantified for abundance and differential expression.
}
\usage{
pipe.RIPpeaksToAltGeneMap(sampleIDset, annotationFile = "Annotation.txt", 
		optionsFile = "Options.txt", speciesID = NULL, results.path = NULL, 
		max.pvalue = 0.05)
}
\arguments{
  \item{sampleIDset}{
Vector of SampleIDs for merging RIP peak data.  
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
  \item{max.pvalue}{
The limiting P-value for determining which peaks from each sample to include in the final 
set of peaks.
}
}
\details{
This is second pipeline step for doing RIP peak analyses.  Given a set of samples each with
their own files of RIP peak results already generated, this function combines them all into 
one merged set of RIP peaks. And then returns those final peaks in a "Gene Map" format that is
directly usable by all the transcriptome and differential expression pipes.

Peaks from different samples that overlap on the same strand are merged such the the final 
peak encloses all original smaller peaks.  A peak need only occur in any one sample for it to
be represented in the final set.  In this way, any RIP peak feature that was detected in any
sample will be a location that gets quantified and compared across sample groups in the 
downstream analyses.
}
\value{
A data frame of RIP peak locations, in the format of a Gene Map, directly usable by downstream
tools like \code{\link{pipe.Transcriptome}} and \code{\link{pipe.DiffExpression}} and 
\code{\link{pipe.MetaResults}}.
Columns include:
  \item{GENE_ID}{the unique identifier for this RIP peak.  It combines the most proximal GeneID 
  with the center point of this RIP peak}
  \item{POSITION}{the starting edge of the RIP peak region to be evaluated}
  \item{END}{the stopping edge of the RIP peak region to be evaluated}
  \item{SEQ_ID}{the chromosome name where this peak is located}
  \item{STRAND}{which chromosome strand the peak was found on.  Note that RIP peaks are always
  strand specific, so there may be 2 separate peaks covering the same genomic location}
  \item{NAME}{the name of the most proximal GeneID}
  \item{PRODUCT}{the gene product term from the most proximal GeneID}
  \item{N_EXON_BASES}{the size of the peak in nucleotides}
  \item{also}{other columns that identify which one sample this RIP peak was most detected in,
  with all the peak scoring details from that sample's RIP peak results.  Intended just for 
  information purposes to tie final peaks back to their source sample, not used by any downstream
  tools}
}
\author{
Bob Morrison
}
\seealso{
  \code{\link{pipe.RIPpeaks}} to call RIP peaks for a single sample.
  \code{\link{MapSets}} for background on Gene Maps.
}
