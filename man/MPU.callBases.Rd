\name{MPU.callBases}
\alias{MPU.callBases}
\alias{MPU.strandDepth}
\alias{MPU.callStringToTable}
\alias{MPU.callTableToString}
\alias{MPU.callStringsToMatrix}
\alias{MPU.callMatrixToBaseCounts}
\title{
MPILEUP Base Call Utility Functions
}
\description{
A collection of low level functions that manipulate the cryptic base call information 
extracted from BAM files by the SAMTOOLS MPILEUP utility.  Most are intended to convert
the cryptic notation into more R friendly formats.
}
\usage{
MPU.callBases(baseCalls, referenceBase, minDepth = 3, mode = c("counts", "percents"))
MPU.strandDepth(baseCalls, pos = NULL)
MPU.callStringToTable(callStrings)
MPU.callTableToString(callTables)
MPU.callStringsToMatrix(callStrings)
MPU.callMatrixToBaseCounts(m, referenceBase, indelsToo = FALSE, normalize = FALSE)
}
\arguments{
  \item{baseCalls}{
Vector of character strings from MPILEUP, that are in the cryptic \code{",,,..A,,..A,...,,"} format.
One string for each nucleotide location.
}
  \item{referenceBase}{
Vector of single character bases, same length as \code{baseCalls}, that are the expected
reference genome calls for those nucleotide locations.
}
  \item{minDepth}{
Integer minimum depth of base calls needed to make an Indel consensus.  Only applies to locations
that report at least one insertion or deletion.  Indel nucleotides with too few base calls will
disregard those calls and return the reference base.
}
  \item{mode}{
Returns the base distribution as either integer counts, or numeric percentages that sum to 100\%
}
  \item{pos}{
Option vector of locations to use as the rownames for the resulting data frame.
}
  \item{callStrings}{
Character vector of base call strings in the compact \code{"A=25;C=1;T=1"} format.
}
  \item{callTables}{
A table or list of tables, where each table is a vector of named integer counts, and the names are 
the base calls.
}
  \item{m}{
A matrix of integer base call counts, as returned by \code{MPU.callStringsToMatrix}.
}
}
\details{
Note that most all MPILEUP base call representations use the comma "," as a generic marker for 
"matches the genomic base", and only uses the typical A/C/G/T characters to denote SNPs, alternate 
base calls that do not match the expected genomic base.
}
\value{
Function \code{MPU.callBases} returns a list:
  \item{calls }{ a vector of consensus base calls.  Note that when the consensus is an Indel, that element
  may have zero or more than one character.}
  \item{depth.table }{ a list of tables, where each element is the full distribution of observed base calls 
  at that nucleotide location}

Function \code{MPU.strandDepth} returns a data frame with 2 columns of integers, \code{FWD_DEPTH, REV_DEPTH},
giving the read depth for each strand separately.

Function \code{MPU.callStringToTable} returns a list of base call tables.

Function \code{MPU.callTableToSring} returns a character vector of base call strings.

Function \code{MPU.callStringsToMatrix} returns a matrix of base call depths, with columns for each possible
base call \code{ , A C G T N Indel}, and a row for each nucleotide position.

Function \code{MPU.callMatrixToBaseCounts} returns a matrix of final base call depths, after merging reference
base information with the MPILEUP calls.  The special comma "," notation is removed, the column names are
the typical base calls \code{ A C G T Indel}, and a row for each nucleotide position.  The matrix values may
be either integer counts, or normalized numeric values that sum each row to 100\%
}
\author{
Bob Morrison
}
