\name{pipe.ExtractPipelineSummaryDetails}
\alias{pipe.ExtractPipelineSummaryDetails}
\alias{extractPipelineSummaryDetails}
\title{
Summarize Alignment Metrics
}
\description{
Tabulate the read counts and percentages for all alignment categories.
}
\usage{
pipe.ExtractPipelineSummaryDetails(sampleIDset, optionsFile = "Options.txt", results.path = NULL, verbose = FALSE)
extractPipelineSummaryDetails(sampleIDset, optionsFile = "Options.txt", results.path = NULL, verbose = FALSE)
}
\arguments{
  \item{sampleIDset}{
The set of SampleIDs to evaluate.
}
  \item{optionsFile}{
File of processing options, which specifies all processing
parameters that are not sample specific.  See \code{\link{DuffyNGS_Options}}.
}
  \item{results.path}{
The top level folder path for writing result files to.  By default, read from the Options 
file entry 'results.path'.
}
}
\details{
Collects the read counts for all possible alignment outcomes, from raw input reads to NoHits, Ribo Cleared,
Genomic, and Splices.  Summarizes by both count and percentage.

When the samples are from a mixed species target, also gives a breakdown of aligned reads per species.

Reported counts are alway for unpaired reads, even if the data was processed as paired end data.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (sampleIDset, optionsFile = "Options.txt", results.path = NULL, 
    verbose = FALSE) 
{
    ans <- extractPipelineSummaryDetails(sampleID = sampleIDset, 
        optionsFile = optionsFile, results.path = results.path, 
        verbose = verbose)
    write.table(ans, "Pipeline.SummaryDetails.csv", sep = ",", 
        quote = T, row.names = F)
    return(invisible(ans))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
