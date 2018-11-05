\name{checkBowtie2Version}
\alias{checkBowtie2Version}
\alias{getBowtie2Version}
\title{
Check the Bowtie2 executable
}
\description{
Checks that the Bowtie2 progam is found, executable, and returns its version number.
}
\usage{
checkBowtie2Version()
getBowtie2Version()
}
\details{
The Bowtie executable can be specified as an environment variable \code{BOWTIE2_PROGRAM} 
or with an \code{bowtie2Program} entry in the \code{\link{DuffyNGS_Options}} file.

\code{checkBowtie2Version} does the complete test of findability, executability, and
returns the version number from Bowtie.

\code{getBowtie2Version} just returns the currently stored version number.
}
\value{
The version number as a character string, if successful.
}
