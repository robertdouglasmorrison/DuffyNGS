\name{buildSpliceJunctionFiles}
\alias{buildSpliceJunctionFiles}
\title{
Build Splice Junction Index and Maps
}
\description{
Create the Bowtie Index and maps for the Splice Junctions portion
of the alignment pipeline.
}
\usage{
buildSpliceJunctionFiles(speciesID, genomicFastaFile, 
	outPath = ".", spliceMapPrefix = "spliceMap", 
	fragmentSize = 60, maxExonSkip = 2, intronSplicesToo = TRUE, 
	verbose = TRUE)
}
\arguments{
  \item{speciesID}{
the SpeciesID for this organism
}
  \item{genomicFastaFile}{
the file of FASTA sequences, or a folder of FASTA files, for this organism.
See \code{\link{getFastaSeqFromFilePath}}
}
  \item{outPath}{
destination folder for the created SpliceJunction map files
}
  \item{spliceMapPrefix}{
a character string prefix for each created map file
}
  \item{fragmentSize}{
the number of bases in each exon end fragment, to become half of each
splice junction entry in the index.  Typically about 90\% of the length
of a typical raw read.
}
  \item{maxExonSkip}{
maximum number of skipped exons, to add exon skipping splice junctions
to the index
}
  \item{intronSplicesToo}{
logical, should intron-exon and intron-intron alternate splices also be
added to the index
}
}
\details{
The Splice Junction index is composed of lots of tiny sequences that
represent the standard and alternative splicing of exons, as defined
by the \code{\link{ExonMap}}.  Each sequence extends \code{fragmentSize}
bases into each of the two exons(introns) being joined.
}
\value{
a Splice Junction Bowtie index and a splice map for each SeqID (chromosome)
are written to disk.
}
\note{
Indexes and maps may need to be moved to the folder of Bowtie indexes.
}
\author{
Bob Morrison
}
