\name{DuffyNGS}
\alias{DuffyNGS}
\alias{DuffyNGS-package}
\title{
Duffy Lab Next Gen Sequencing "NGS" package Overview
}
\description{
A suite of R tools from the lab of Patrick Duffy, providing an integrated
infrastucture for the processing and analysis of NextGen Sequencing datasets.
}
\details{
This package provides a suite of utilities and wrapper functions that harness
the power of industry standard NGS tools such as Bowtie2, SAMTOOLS, Velvet, etc.,
into a self-contained R package that leverages the genome annotation features and mixed
organism support provived by the \code{\link{DuffyTools}} utilities package.  

The goal is to provide simple R functions that encapsulate and implement complex 
bioinformatic computational workflows into small managable 'pipes' that streamline
the overall processing workflow into discrete steps with minimal input from the 
scientist or bioinformatician.

One primary focus is on providing a standardized yet flexible structure for naming and
managing file names and folder layouts, to simplify the organization of the myriad
types of result files that eminate from NGS datasets.  This becomes essential in mixed
organism experiments where the sample contains both host and pathogen, and each NGS data
will produce transcriptomes and variant calls for multiple species.

A second main focus is to simplify the combining of the strengths and weaknesses of many published
bioinformatic tools, to give more robust measures of gene and pathway differential
expression.  Included are pipes that implement many public algorithms, such as EdgeR,
DESeq, SAM, RankProduct, QuSage, Enrichment, and then various averaging tools that merge 
and combine those disparate results to produce 'meta' results that better reflect the 
overall consensus of multiple algorithmic approaches.
}
\author{
Bob Morrison  bob.morrison@seattlebiomed.org   robert.morrison@nih.gov
}
