# pipe.RemoveDiscordantAlignsFromNoHits.R

# when doing paired end alignments, Bowtie sends all non-concordant pairs to the no hit files
# even if one of the mates did align discordantly.  Thus some reads are in both the 'aligned'
# and the 'NoHits' results.   Remove those from the 'NoHits'...

`pipe.RemoveDiscordantAlignsFromNoHits` <- function( sampleID, file=NULL, 
						annotationFile="Annotation.txt", optionsFile="Options.txt",
						mode=c("pairedNoHits", "mergedNoHits"), verbose=TRUE) {

	if (verbose) cat( "\nRemoving 'Discordant Alignments' from NoHits..:  ", sampleID)

	# set up file names for final results and temp Bowtie call
	resultsPath <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	noHitFileIn <- paste( sampleID,"noHits.fastq.gz", sep=".")
	finalNoHitFileOut <- paste( sampleID,"noHits.fastq.gz", sep=".")
	finalNoHitFileOut <- file.path( resultsPath, "fastq", finalNoHitFileOut)
	tempResultAlignFile <- paste( sampleID,"DiscordantAligns.bam", sep=".")
	tempResultNoHitsFile <- paste( sampleID,"noHits.Temp.fastq.gz", sep=".")

	# input 'No Hits' file(s) may be already named explicitly...
	mode <- match.arg( mode)

	if ( mode == "pairedNoHits") {
		# it may be paired end data that made paired end nohits
		# files are in current working directory
		if ( is.null( file)) {
			noHitFileIn <- paste( sub( ".gz$", "", noHitFileIn), 1:2, "gz", sep=".")
		} else {
			noHitFileIn <- file
		}
	} else {
		if ( is.null( file)) {
			# otherwise its already the final format/folder
			noHitFileIn <- file.path( resultsPath, "fastq", noHitFileIn)
		} else {
			noHitFileIn <- file
		}
	}
	if ( ! all( file.exists( noHitFileIn))) {
		cat( "\nInput file(s) of discordant reads not found: ", noHitFileIn)
		return( NULL)
	}

	# call bowtie against the genomic index
	ans <- fastqToBAM( noHitFileIn, tempResultAlignFile, sampleID=sampleID, optionsFile=optionsFile,
			annotationFile=annotationFile, noHitsFile=tempResultNoHitsFile, asMatePairs=FALSE)
	if (verbose) {
		cat( "  Done.\n")
	}

	# clean up...
	file.delete( finalNoHitFileOut)
	file.delete( tempResultAlignFile)
	file.delete( paste( sampleID, "Bowtie2.AlignMetrics.txt",sep="."))
	file.delete( noHitFileIn)
	file.rename( tempResultNoHitsFile, finalNoHitFileOut)

	return( ans)
}

