# pipe.FileCleanup.R

# do the file combining, renaming, deleting temps, etc., after the alignment pipeline has run.

`pipe.FileCleanup` <- function( sampleID,  optionsFile="Options.txt", verbose=FALSE) {

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\n\nStarting 'File Cleanup' for Sample:     ", sampleID, "\n\n")
	}

	results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)

	# delete files we don't need
	if (verbose) cat("\nDeleting temp files...")
	file.delete( paste( sampleID, "not", c("ribo","genomic","splice"), "fastq", sep="."))
	file.delete( paste( sampleID, "not", c("ribo","genomic","splice"), "fq", sep="."))
	file.delete( paste( sampleID, "not", c("ribo","genomic","splice"), "fastq.gz", sep="."))
	file.delete( paste( sampleID, "not", c("ribo","genomic","splice"), "fq.gz", sep="."))
	file.delete( paste( sampleID, "not", c("ribo","genomic","splice"), "fastq", rep( 1:2, each=3), "gz", sep="."))
	file.delete( paste( sampleID, "not", c("ribo","genomic","splice"), "fq", rep( 1:2, each=3), "gz", sep="."))

	# see if any BAM conversions went longer than their genomic parts...
	if (verbose) cat("\nDeleting BAM conversion log files...")
	file.delete( paste( sampleID, "convertRiboBAM.log.txt", sep="."))
	file.delete( paste( sampleID, "convertingRiboBAM.done", sep="."))
	file.delete( paste( sampleID, "convertSpliceBAM.log.txt", sep="."))
	file.delete( paste( sampleID, "convertingSpliceBAM.done", sep="."))

	# see if any paired end stranded result file pairs need to be re-joined
	genomeFilesIn <- file.path( results.path, "align", paste( sampleID, "_", 1:2, ".genomic.sorted.bam", sep=""))
	genomeFileOut <- file.path( results.path, "align", paste( sampleID, ".genomic.sorted.bam", sep=""))
	if ( all( file.exists( genomeFilesIn))) {
		if (verbose) cat("\nMerging separate mate genomic BAM files...")
		BAM.merge( genomeFilesIn, genomeFileOut, index=TRUE, verbose=T)
		lc1 <- quickFileLineCountLookup( genomeFilesIn[1], sampleID=paste(sampleID,"1",sep=""), what="lineCount")
		lc2 <- quickFileLineCountLookup( genomeFilesIn[2], sampleID=paste(sampleID,"2",sep=""), what="lineCount")
		nreads <- nlines <- lc1 + lc2
		quickFileLineCountRecord( genomeFileOut, sampleID=sampleID, lineCount=nlines, readCount=nreads)
		file.delete( genomeFilesIn)
	}

	spliceFilesIn <- file.path( results.path, "splicing", paste( sampleID, "_", 1:2, ".splice.converted.sorted.bam", sep=""))
	spliceFileOut <- file.path( results.path, "splicing", paste( sampleID, ".splice.converted.sorted.bam", sep=""))
	if ( all( file.exists( spliceFilesIn))) {
		if (verbose) cat("\nMerging separate mate splice BAM files...")
		BAM.merge( spliceFilesIn, spliceFileOut, index=TRUE, verbose=T)
		lc1 <- quickFileLineCountLookup( spliceFilesIn[1], sampleID=paste(sampleID,"1",sep=""), what="lineCount")
		lc2 <- quickFileLineCountLookup( spliceFilesIn[2], sampleID=paste(sampleID,"2",sep=""), what="lineCount")
		nreads <- nlines <- lc1 + lc2
		quickFileLineCountRecord( spliceFileOut, sampleID=sampleID, lineCount=nlines, readCount=nreads)
		file.delete( spliceFilesIn)
	}
	if ( verbose) cat( "\n...Done.\n")
}

