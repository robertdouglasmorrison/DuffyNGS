\name{bamReadStats}
\alias{bamReadStats}
\alias{fastqReadStats}
\alias{readStats}
\title{
Quality Control Statistics of Raw Sequencing Reads
}
\description{
Assess the various quality metrics of a file of raw sequencing reads, looking at
read length, Phred scores, nucleotide distributions, flow cell tile variation, etc.
}
\usage{
bamReadStats(filein, sampleID, statsPath = "BamReadStats", 
	calcStats = TRUE, plotAllTiles = FALSE, baseOrder = "ACGTN", 
	chunkSize = 1e+05, maxReads = NULL, pause = 0)

fastqReadStats(filein, sampleID, statsPath = "FastqReadStats", 
	calcStats = TRUE, plotAllTiles = FALSE, baseOrder = "ACGTN", 
	chunkSize = 1e+05, maxReads = NULL, pause = 0)

readStats(filein, sampleID, ...)

}
\arguments{
  \item{filein}{
Full pathname to an existing file of raw sequencing reads, may be 'gzip' 
compressed for FASTQ.
}
  \item{sampleID}{
SampleID for this file, used as a prefix on the names of created files and plots.
}
  \item{statsPath}{
Destination folder to receive the created plots and summary data.
}
  \item{calcStats}{
Logical, either calculate the statistics or just replot using a previous file
of statistics.
}
  \item{plotAllTiles}{
Create separate plots of each tile in the file, in addition to the one overall
plot per flow cell lane.
}
  \item{baseOrder}{
Order to display the nucleotides in the base call bar plots.
}
  \item{chunkSize}{
Buffer size, in reads, for processing the file.  Small buffer size yields
rapid plots and updates of progress, slower overall performance, and uses
less memory.
}
  \item{maxReads}{
The maximum number of reads to process, \code{NULL} means use all.  Note that
sorted BAM files place un-aligned reads at the end.
}
  \item{pause}{
Delay in seconds for viewing each plot.
}
}
\details{
This function analyzes several metrics about the raw sequencing file, to assess the quality
of the sequencing run, and to help select suitable alignment parameters for tuning the
alignment pipeline.  With default arguments, it generates about 4 plots per sample,
focusing on Phred scores, nucleotide distributions, and the variance between tiles.

There are separate functions for BAM and FASTQ files, and a wrapper function \code{readStats} that
uses the file extension to dispatch based on the given file type.
}
\value{
In addition to the plots created, one rather complex data object is written to disk,
containing the details that are used to generate the plots.  It can be loaded to extract
those numeric details.
}
\author{
Bob Morrison
}
\note{
With each new release of sequencing machine software, the details about the coordinates
of each read change inside the ReadID.  Parsing out the "lane:tile:X:Y" terms is a
perpetual work in progress.  This tool may break/fail on new ReadID formats.
}
