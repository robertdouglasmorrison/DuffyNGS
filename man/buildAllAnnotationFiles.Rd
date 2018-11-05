\name{buildAllAnnotationFiles}
\alias{buildAllAnnotationFiles}
\title{
Create All Annotation-dependant Files for a DuffyNGS Target Organism(s)
}
\description{
Build all needed files for running the various Bowtie pipeline pieces for
a target organism(s).   This will construct genomic Bowtie indexes, RiboClearing
maps and indexes, Splice Junction maps and indexs, and cross species detectability
objects.
}
\usage{
buildAllAnnotationFiles(speciesIDSet, genomicFastaFileSet, 
	fastaPatternSet = "\\.fasta$", optionsFileSet = "Options.txt", 
	outPath = ".", steps = 1:4, riboTailSize = 36, 
	exonFragmentSize = 36, maxExonSkipSet = 2, intronSplicesTooSet = FALSE, 
	detectableReadSize = 32, chunkSize = 1e+06, 
	otherSpeciesIDSet = "Hs_grc", otherGenomicIndexSet = "Hs.genomic_idx", 
	verbose = TRUE, debug = FALSE)
}
\arguments{
  \item{speciesIDSet}{
a character vector of the SpeciesIDs that comprise one new target
}
  \item{genomicFastaFileSet}{
a character vector (same length as \code{speciesIDset}) of the FASTA files or folders
that are the genomic definition of each species.
}
  \item{fastaPatternSet}{
a character vector (same length again) of regular expressions for finding FASTA files
in the folder, for species that have each chromosome in separate files.
}
  \item{optionsFileSet}{
a character vector of OptionsFile names.  If more than one species, this vector must
have an options file for each, plus one extra for the combined target.
}
  \item{outPath}{
destination folder for local copies of all created files
}
  \item{steps}{
See details.  A vector in integers in the range 1:4, to specify which componment files
to build.  Since each step can be quite time consuming (or error prone), they can be
done one at a time or just as new annotation files (MapSets) dictate.
}
  \item{riboTailSize}{
the number of bases to extend past the ends of genes being targeted for RiboClearing.  This
allows reads that overlap the gene ends to still align, and thus be cleared. The value
should be about 75\% of a typical read length.
}
  \item{exonFragmentSize}{
the number of bases into each exon to use for building the Splice Junction index.  The value
should be about 90\% of a typical read length, so any aligned reads must be spanning 2 exons.
}
  \item{maxExonSkipSet}{
the maximum number of exons to be skipped, for the alternative splicing exon pairs.
}
  \item{intronSplicesTooSet}{
logical, build alternative splice junctions using both introns and exons.
}
  \item{detectableReadSize}{
number of bases in a canonical read, for measuring in-species uniqueness and cross-species
detectability.
}
  \item{chunkSize}{
buffer size, for the detectability step
}
  \item{otherSpeciesIDSet}{
a character vector of 'second' SpeciesIDs, for building cross-detectability files
}
  \item{otherGenomicIndexSet}{
a character vector (same length as \code{otherSpeciesIDset}) of the Bowtie genomic indexes
for each 'second' species.
}
  \item{verbose}{
logical, print lots of messages during progress
}
  \item{debug}{
logical, perform sanity checks during processing
}
}
\details{
This is an umbrella function that attempts to build all needed files for a new target
organism(s) and/or new read length data sets.  It builds all needed Bowtie indexes 
(Genomic, RiboClearing, Splice Junctions), all needed alignment decoding maps (RiboClearing,
Splice Junctions), and all 2-species detectability information.

The key focus is to easily combine multiple organisms into one joint target index family
in a way that seemlessly creates all needed files for future analysis of these organisms.

For the non-genomic indexes, the typical raw read lengths impacts the settings for building
these indexes, to optimize the align-ability of the reads. For this reason, it is often
wise to make separate indexes if you will work with very different read lengths (i.e.
separate riboClear and splice junction indexes for short (~40bp) and long (~100bp) reads.)

Step=1:  Build Genomic indexes for each organism, and a joint genomic index if 
\code{speciesIDSet} is longer than one.

Step=2:  Build Ribo Clearing indexes and maps for each organism, and a joint 
RiboClearing index if \code{speciesIDSet} is longer than one.

Step=3:  Build Splice Junction indexes and maps for each organism, and a joint 
SpliceJunction index if \code{speciesIDSet} is longer than one.

Step=4:  Build detectability files.  
}
\value{
A family of files written to disk, that may need to be moved to the appropriate folders
for use in subsequent pipeline runs.
}
\note{
Indexes and maps need to be moved to the folder of Bowtie indexes, and
detectability files need to be moved to the /data subfolder of the the
DuffyNGS package installation location.
}
\author{
Bob Morrison
}
