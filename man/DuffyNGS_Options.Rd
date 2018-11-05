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
}}
