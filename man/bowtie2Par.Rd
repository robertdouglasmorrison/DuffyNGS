\name{bowtie2Par}
\alias{bowtie2Par}
\alias{bowtie2Par.defaults}
\alias{getBowtie2Version}
\title{
Set or Query Bowtie2 Parameters
}
\description{
Set and retrieve named parameters that control the Bowtie2 alignment
program.
}
\usage{
bowtie2Par(...)
bowtie2Par( <name>=<value>)

bowtie2Par.defaults( optionsFile)

getBowtie2Version()
}
\arguments{
  \item{\dots}{
parameters to query, or name-value pairs for setting parameters
}
  \item{optionsFile}{
name of a \code{\link{DuffyNGS_Options}} that contains Bowtie2 parameters to be
set
}
}
\details{
Bowtie2 uses several command line parameters that affect its runtime behaviour,
and this utility provides a common interface to setting and retrieving those
Bowtie2 parameters.  To see all current parameters, use \code{bowtie2Par()}.  

Required parameters at present:

  \code{Program} is the full pathname to the Bowtie2 executable program

  \code{IndexPath} is the pathname to the folder of all needed Bowtie2
index files

  \code{GenomicIndex} is the basename of the Bowtie2 Index to be used for 
the genomic alignment phase of the pipeline

Calling \code{bowtie2Par.defaults()} with no \code{optionsFile} will use values from
the environment variables:  BOWTIE2_PROGRAM, BOWTIE2_INDEX_PATH, and BOWTIE2_INDEX_FILE, 
respectively.
}
\value{
For \code{ bowtie2Par(parameter)} the currently set value of that parameter, or an error
if that parameter is not set.

For \code{getBowtie2Version()}, the current version of Bowtie2 as a character string
}

\seealso{
  \code{\link{fastqToBAM}} for a higher level wrapper to calling Bowtie
  \code{\link{buildBowtie2CommandLine}} for the utility that constructs valid Bowtie2 command lines.
  \code{\link{callBowtie2}} for envoking the Bowtie2 alignment program.
  \code{\link{buildBowtie2BuildCommandLine}} for the utility that constructs valid Bowtie2-build command lines.
  \code{\link{callBowtie2Build}} for envoking the Bowtie2-build index creation program.
}

