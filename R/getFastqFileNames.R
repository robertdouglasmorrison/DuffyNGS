# getFastqFileNames.R -- utility to get the various Fastq file names


# the raw fastq files can be multiples, and perhaps stranded pairs
getRawFastqFileNames <- function( sampleID, annotationFile, optionsFile, verbose=T) {

	annT <- readAnnotationTable( annotationFile)
	rawFastqFile <- getAnnotationValue( annT, sampleID, "Filename")

	# grab and show a few files, to force auto-mount of file system
	fastq.path <- getOptionValue( optionsFile, "fastqData.path", verbose=verbose)
	dirList <- dir( fastq.path)
	if (verbose) cat( "\n\n'Fastq' Folder: ", fastq.path, "\nN_Files Found: ", length(dirList),  
			"\n", head(dirList), "  etc.")

	# allow paired end and strand specific reads
	pairedEnd <- getAnnotationTrue( annT, key=sampleID, "PairedEnd", notfound=FALSE, verbose=verbose)
	strandSpecific <- getAnnotationTrue( annT, key=sampleID, "StrandSpecific", notfound=FALSE, verbose=verbose)
	doPairs <- (pairedEnd && strandSpecific)

	# separate the files if more than one
	rawFastqFileSet <- strsplit( rawFastqFile, split=", *")[[1]]
	if (pairedEnd) {
		nFastq <- 2
		if ( length( rawFastqFileSet) != nFastq) 
				stop( "Fastq filename should be 2 comma separated names if paired end reads")
	} else {
		nFastq <- length( rawFastqFileSet)
	}

	# allow compression
	rawFastqFileSet <- file.path( fastq.path, rawFastqFileSet)
	for ( i in 1:nFastq) rawFastqFileSet[i] <- allowCompressedFileName( rawFastqFileSet[i])
	
	out <- list( "files"=rawFastqFileSet, "asMatePairs"=doPairs)
	return( out)
}


`getNotRiboFastqFileNames` <- function( sampleID, asMatePairs=FALSE) {

	if (asMatePairs) {
		notRiboFiles <- paste( sampleID, "not.ribo.fq", 1:2, "gz", sep=".")
	} else {
		notRiboFiles <- paste( sampleID, "not.ribo.fq.gz", sep=".")
	}
	return( notRiboFiles)
}


`getNotGenomicFastqFileNames` <- function( sampleID, asMatePairs=FALSE) {

	if (asMatePairs) {
		notGenomicFiles <- paste( sampleID, "not.genomic.fq", 1:2, "gz", sep=".")
	} else {
		notGenomicFiles <- paste( sampleID, "not.genomic.fq.gz", sep=".")
	}
	return( notGenomicFiles)
}


`getNotSpliceFastqFileNames` <- function( sampleID, asMatePairs=FALSE) {

	if (asMatePairs) {
		notSpliceFiles <- paste( sampleID, "not.splice.fq", 1:2, "gz", sep=".")
	} else {
		notSpliceFiles <- paste( sampleID, "not.splice.fq.gz", sep=".")
	}
	return( notSpliceFiles)
}

