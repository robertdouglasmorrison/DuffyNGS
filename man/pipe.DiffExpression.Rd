\name{pipe.DiffExpression}
\alias{pipe.DiffExpression}
\alias{runDiffExpression}
\title{
Turn Pairs of Transcriptomes into Differential Gene Expression Ratio Files
}
\description{
This pipe takes transcriptome files as input, and creates ratio files that measure the
differential gene expression as both fold change and P-values.  Ratio files are then used as
input to the Round Robin differential expression pipe.
}
\usage{
pipe.DiffExpression(sampleIDset, groupSet = NULL, annotationFile = "Annotation.txt", 
		optionsFile = "Options.txt", results.path = NULL, speciesID = NULL, 
		altGeneMap = NULL, altGeneMapLabel = NULL, minRPKM = NULL, 
		missingOnly = FALSE, verbose = TRUE)
}
\arguments{
  \item{sampleIDset}{
Vector of SampleIDs, giving the set of transcriptome samples to use for making ratio files.
}
  \item{groupSet}{
Vector of grouping labels, the same length as \code{sampleIDset}, giving the group membership
for each sample.  Only pairs of samples from different groups are turned into ratio files.  Pairs
of samples from the same group are silently ignored.
}
  \item{annotationFile}{
File of sample annotation details, which specifies all needed 
sample-specific information about the samples under study.  
See \code{\link{DuffyNGS_Annotation}}.
}
  \item{optionsFile}{
File of processing options, which specifies all processing
parameters that are not sample specific.  See \code{\link{DuffyNGS_Options}}.
}
  \item{results.path}{
The top level folder path for writing result files to.  By default, read from the Options 
file entry 'results.path'.
}
  \item{speciesID}{
The SpeciesID of the target species to calculate ratio files for.  By default, transcriptomes 
from all species in the current target are used to generate ratio files for all species.
}
  \item{altGeneMap}{
An alternate data frame of gene annotations, in a format identical to \code{\link{getCurrentGeneMap}},
that has the gene names and locations to be measured for read pileup depth.  By default, use the 
standard built-in gene map for each species.
}
  \item{altGeneMapLabel}{
A character string identifier for the alternate gene map, that becomes part of all created path and file names
to indicate the gene map that produced those transcriptomes.
}
  \item{minRPKM}{
Numeric offset that is added to all gene expression values before calculating the fold change ratio,
to prevent division by zero and extremely high ratios from extremely low expression values.  By default, 
takes its value from the \code{DE.minimumRPKM} option. Note the option can use a species specific 
suffix to allow automatic use of different minimums for different species in a mixed organism experiment.
}
  \item{missingOnly}{
Logical, controls whether all ratio files are re-calculated, or if only those that are not yet found in
the ratios folder get created.
}
}
\details{
This pipe is perhaps slightly misnamed, in that it only does mnay simple 2-sample comparisons, to create 
the ratio files needed by the Round Robin pipe, the oldest and originally only DE tool in the DuffyNGS suite.
True group-wise differential expression pipes are mentioned below.
}
\value{
A family of files is created in the \code{ratio} subfolder, containing the 2-way fold change differential
expression comparison between transcriptomes from 2 different group categories.
}
\author{
Bob Morrison
}
\seealso{
  \code{\link{pipe.DESeq}} {}
  \code{\link{pipe.EdgeR}} {}
  \code{\link{pipe.RankProduct}} {}
  \code{\link{pipe.RoundRobin}} {}
  \code{\link{pipe.SAM}} { for the 5 separate Differential Expression tools that do group versus group
  comparisons.}

  \code{\link{pipe.MetaResuts}} { for the dispatching and final merging & averaging of all 5 DE tools.}
}
