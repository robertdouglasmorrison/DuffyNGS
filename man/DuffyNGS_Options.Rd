\name{DuffyNGS_Options}
\alias{DuffyNGS_Options}
\title{
Options File of Processing Settings 
}
\description{
The Options file defines the various processing settings for
all datasets to be processed.  See \code{\link{DuffyNGS_Annotation}} for 
settings that are specific to each sample.  There are 2 columns 
\code{ OptionName, Value} and all entries are tab-delimited.
See \code{\link{OptionsTable}} for more information.
}
\details{
The motivation behind the Options file is to allow the end user a large degree
of control over the run time behavior of the DuffyNGS pipeline.  A fixed set of
parameters and command arguments can be altered by changing the entries in this
file.  If there are addititional fields that you want finer control over, contact the
authors.
}
\section{Options}{
  \describe{
  \item{targetID}{
The TargetID for processing the samples.  This specifies the set of species that the raw
data will be aligned against, etc.  See \code{\link{Targets}}.   Default is \code{HsPf}
}
  \item{nCores}{
Enables parallel execution using multiple cores.  Sets the maximum number of cores to
use simultaneously.  Most low level functions that operate on chromosomes are written
to do chromosomes in parallel. See \code{\link{multicore.setup}}.
}
  \item{genomicFastaFile}{
Full pathname to the file (or folder) of genomic DNA for this set of species.  In the
case of very large genomes, it can be the folder name that contains FASTA files for
each chromosome.
}
  \item{fastqData.path}{
Full pathname of the folder that contains the raw read datasets to be aligned.
This will be combined with the actual filename from the Annotation file to build the
complete pathname for a sample's raw FASTQ files.
}
  \item{results.path}{
Full pathname of the folder that will receive all generated results files for all
processed samples.  It will be created if it does not yet exist.  Default is \code{"."}
}
  \item{tmp.path}{
Explicit full pathname to temporary storage on the current compute node.  Used for
temporary large scratch disk space.
}
  \item{bowtie2Program}{
Full pathname of the Bowtie2 executable to be used for alignments.
}
  \item{bowtie2Index.path}{
Full pathname to the folder of Bowtie2 alignment indexes.  The final filename for an index will
be constructed from this path and the named index for each alignment step.  See below.
}
  \item{bowtie2InputOptions}{
A character string of Bowtie parameters affecting how the raw reads are interpreted.
This includes hard trimming, fastq format, quality score type, etc.
}
  \item{bowtie2PairOptions}{
A character string of Bowtie parameters affecting paired end behavior.
}
  \item{trim5}{}
  \item{trim3}{
Integer.  Hard trimming of bases.  These will be passed to Bowtie2, and are used by other 
DuffyNGS tools that need to treat the raw reads exactly as Bowtie did.
}
  \item{maxMultiHits}{
Integer.  Controls the maximum number of chromosomal locations reported by Bowtie for an aligned read.
Passed to Bowtie as its '-k' options.
}
  \item{readBufferSize}{
Integer.  Controls memory usage for pipeline steps that operate on reads in a buffered manner.  Larger
size will run faster but consume more memory.
}
  \item{DE.minimumRPKM}{
Numeric.  Used to prevent dividing by very low or zero expression values.  Added as an offset to most 
ratio and fold change calculations.  Sensible values vary by the size of the genome.  If set too high then
true fold change will be under reported.  If set too low, then unexpressed genes may be reported as 
having extremely large fold change.
}
  \item{DE.weightFoldChange}{}
  \item{DE.weightPvalue}{
Numeric weights that set the relative importance of these two primary metrics for ordering differentially 
expressed genes.  See \code{\link{diffExpressRankOrder}}.
}
}}
