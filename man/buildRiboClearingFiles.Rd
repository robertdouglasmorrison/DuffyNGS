\name{buildRiboClearingFiles}
\alias{buildRiboClearingFiles}
\title{
Build Ribo Clearing Index and Map
}
\description{
Create the Bowtie Index and map for the RiboClearing portion
of the alignment pipeline.
}
\usage{
buildRiboClearingFiles(speciesID, genomicFastaFile, 
	outPath = ".", mapFile = "riboMap.txt", 
	tailSize = 24, verbose = TRUE)
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
destination folder for the created RiboClearing map file
}
  \item{mapFile}{
name for the newly created RiboClearing map file
}
  \item{tailSize}{
number of bases to extend beyond the gene edges for constructing the
Bowtie index of the genes to be cleared.  This allows reads that span the
gene boundary to be matched by the index and thus cleared.
}
}
\details{
Exactly which genes go into the RiboClearing index is determined by the
\code{RrnaMap}.  See \code{\link{MapSets}}.
}
\value{
a Ribo Clearing Bowtie index and RiboClear map are written to disk
}
\note{
Indexes and maps may need to be moved to the folder of Bowtie indexes.
}
\author{
Bob Morrison
}
