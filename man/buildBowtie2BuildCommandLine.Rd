\name{buildBowtie2BuildCommandLine}
\alias{buildBowtie2BuildCommandLine}
\title{
Build a Bowtie2-build Command Line
}
\description{
Construct a Bowtie2-build command line that can create a Bowtie2 index file.
}
\usage{
buildBowtie2BuildCommandLine(inputFastaFile, outputIndexFile, 
	optionsFile = "Options.txt", verbose = TRUE, debug = FALSE)
}
\arguments{
  \item{inputFastaFile}{
a character string that contains comma separated full pathnames to existing FASTA files
}
  \item{outputIndexFile}{
base name for the set of index files to be created
}
  \item{optionsFile}{
the options file.  This invokes \code{\link{bowtie2Par.defaults}} to set the location
of the Bowtie2-build executable.
}
}
\note{
This blindly assumes that the 'bowtie2-build' executable is in the same folder as
'bowtie2'.
}
\value{
a character string, suitable for passing to \code{\link{callBowtie2Build}} or to
\code{\link{system}}.
}
