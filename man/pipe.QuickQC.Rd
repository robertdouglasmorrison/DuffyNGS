\name{pipe.QuickQC}
\alias{pipe.QuickQC}
\title{
Quick Quality Control 'A to Z' Alignment & Transcriptome.
}
\description{
Complete alignment pipeline for a quick QC look at a sample.  Runs every pipeline step on a
subset of the raw data, to generate initial assessment metrics about every stage of the
alignment and transcriptome process. Useful for verifying all runtime options, sample 
annotation settings, needed environment variables, etc.
}
\usage{
pipe.QuickQC(sampleID, annotationFile = "Annotation.txt", optionsFile = "Options.txt", 
		results.path = "QuickQC", banner = "QuickQC", chunkSize = 1e+06, 
		maxReads = 5e+06, pause = 0, verbose = T, nUSRkeep = 1e+06, 
		mode = c("all", "alignOnly", "pieOnly"))
}
%- maybe also 'usage' for other objects documented here.
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
  \item{results.path}{
The top level folder path for writing result files to.  By default, creates a new folder called 
'QuickQC'.
}
  \item{mode}{
Controls the behavior of which parts of the QC pipeline get run.
}
  \item{banner}{
Character string that gets appended to the headers of all created plot images.
}
  \item{chunkSize}{
Numeric size of the buffer used for processing aligned reads.  Smaller values give faster more 
frequent updates and QC images, but slow down overall job runtime.
}
  \item{maxReads}{
Numeric, maximum number of FASTQ reads to process. If the data is paired end, does half that many 
reads from each mate pair file separately.
}
  \item{pause}{
Optional delay in seconds, between plotted images, to give more time to visually inspect QC plots.
}
  \item{nUSRkeep}{
Number of Unique Short Reads (USR) to process when interrogating the 'NoHits' reads for discovery 
of foreign contaminants, etc.
}
}
\details{
Since the full pipeline can take 10 to 36 hours to complete for one sample, that is far too long to wait to 
verify that all option and annotation settings are optimal.  This QC pipeline is intended to give
fast feedback about all aspects of the sample dataset, runtime settings, alignment success, and the 
user's computation environment.

Typically, with each new batch of experimental data coming off the sequencer, it is wise to run 
one or two samples through the QC pipe to verify that the sequencing step worked as expected.  Typical 
things that can vary between experiments include the need for base trimming at the 5' and 3' ends,
the amount of ribosomal RNA, raw read sense, etc.

Once a few QC runs have completed without error, and the transcriptomes and statistics look as expected,
then all samples can confidently be submitted to the cluster as batch jobs.
}
\value{
A large variety of subfolders and files are written to disk. Each subfolder holds a family of plot images 
that convey aspects of the quality of the raw FASTQ reads, the alignments, gene expression levels, and an
assessment of possible foreign contaminants.

  \item{FastqReadStat }{files of plot images that assess the raw base calls and quality scores in 
  the FASTQ files.  Metrics about the A/C/G/T/N distribution and variance, as well as metrics about the
  tile to tile variance on the slide are reported as PNG files.}

  \item{AlignStats }{files of plot images that assess many aspects of read alignment details from the BAM files.
  Metrics about:  base mismatches, including location within the read and A/C/G/T/N identities;  gene and 
  chromosome coverage and distribution; insertions and deletions within the reads; Mulitply Aligned Reads (MARs);
  species proportions if mixed organism target; and the distribution of alignment quality scores.}

  \item{transcript }{transcriptome files of gene expression.}

  \item{html }{HTML files of the top expressing genes in each species, with hyperlinks to read pileup images
  that display read alignment strand calls and chromosomal coverage.}

  \item{CR }{'Consensus Reads' de novo assembly of No Hit reads.  Raw FASTQ reads that failed to align get
  analyzed in a two step process:   1. grow the reads into larger 'consensus read' contigs that share 
  overlapping DNA sequences, to create larger pieces of DNA that represent unknown material that was 
  clearly present in the FASTQ file yet failed to align to the target genome(s);  2.  Blast these DNA contigs
  against the giant NCBI 'NT' database, to get an unbiased call for what those DNA contigs most likely truly
  are.  The result is text and HTML files with hyperlinks to CR pileup images, and a PNG image of the overall
  proportions of foreign material identified in the sample.}
}
\author{
Bob Morrison
}
\note{
The 'CR' step assumes there is a locally installed version of BLAST and an NT Blast database, as configured 
in the Options file and environment variables.
}
