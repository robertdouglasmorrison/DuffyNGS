\name{pipe.VariantSummary}
\alias{pipe.VariantSummary}
\title{
Summarize all Chromosomes of SNP Calls into one File
}
\description{
Combines all the separate chromosome SNP call results for a sample into a single file, with optional
filtering by exon, score, etc.
}
\usage{
pipe.VariantSummary(sampleID, speciesID = getCurrentSpecies(), annotationFile = "Annotation.txt", 
		optionsFile = "Options.txt", results.path = NULL, seqIDset = NULL, min.depth = 1, 
		min.score = 5, exonOnly = FALSE, snpOnly = FALSE)
}
\arguments{
  \item{sampleID}{
Character string of one SampleID that already has SNP calls done.
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
  \item{speciesID}{
The SpeciesID of the target species to call SNPs for.  By default, use the current species.
}
  \item{results.path}{
The top level folder path for writing result files to.  By default, read from the Options 
file entry 'results.path'.
}
  \item{seqIDset}{
Optional character vector of SeqIDs. Default is to combine and summarize SNPs from all chromosome.
}
  \item{min.depth}{}
  \item{min.score}{
  Filtering arguments, to only keep SNP calls that meet minimum criteria about depth of read coverage
  and the BCFTOOLS CALL score metric.
}
  \item{exonOnly}{
Optional filtering to remove all SNP calls that lie in intergenic or intron regions.  Intergenic SNPs are
not traditionally as interesting.  And SNPs in introns are vary likely false SNPs due to alignment
methods and/or poor genome annotation of true exon boundaries.
}
  \item{snpOnly}{
Options filtering to remove all SNP calls that do not see an alternate allele as the most frequent
observed base.  There may be cases, like mixed infections or due to the limitations of alignment precision, 
where a minor allele gets flagged as a SNP, but has less the 50\% of the observed base calls.  So the site
is both called a SNP and yet still matches the reference base call. This option removes these 'paradoxical'
SNP calls.
}
}
\value{
One file is written under the \code{VariantCalls} subfolder:
  \item{Summary.VCF.txt }{One final file of SNP sites, after merging all chromosomes and cleaning
  up much of the BCFTOOLS details.  Includes a column "ALT_AA" that tries to suggest if the SNP
  changes the amino acid sequence of the protein.}
}
\author{
Bob Morrison
}
\seealso{
  \code{\link{pipe.VariantCalls}} {for doing the initial SNP calls on each chromosome.}
  
  \code{\link{pipe.VariantComparison}} {for finding SNPs that are diffentially detected between groups.}
}
