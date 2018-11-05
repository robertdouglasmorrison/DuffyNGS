# pipe.FastqReadStats.R

`pipe.FastqReadStats` <- function( sampleID, annotationFile="Annotation.txt",
				optionsFile="Options.txt", results.path=NULL,
				chunkSize=200000, maxReads=NULL, fastqFile=NULL, ...) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		resultsPath <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	} else {
		resultsPath <- results.path
	}
	fastqPath <- getOptionValue( optT, "fastqData.path", notfound=".")

	statsPath <- file.path( resultsPath, "FastqReadStats")
	if ( ! file.exists( statsPath)) dir.create( statsPath, recursive=TRUE)

	if ( is.null( fastqFile)) {
		fqFileIn <- getAnnotationValue( annotationFile, key=sampleID, columnArg="Filename")
	} else {
		fqFileIn <- fastqFile
	}

	# allow for paired end file pair
	fqFileSet <- strsplit( fqFileIn, split=", *")[[1]]
	sampleIDs <- sampleID
	NS <- length(fqFileSet)
	if ( NS > 1) {
		sampleIDs <- paste( sampleID, 1:NS, sep="_")
	}

	for ( iSample in 1:NS) {

		fqFile <- fqFileSet[iSample]
		sampleID <- sampleIDs[iSample]
	
		# if not given the fastq file explicitly, turn it to its full pathname
		if ( is.null( fastqFile)) {
			fqFile <- file.path( fastqPath, fqFile)
			fqFile <- allowCompressedFileName( fqFile)
		}

		fastqReadStats( fqFile, sampleID=sampleID, statsPath=statsPath, chunkSize=chunkSize, 
				maxReads=maxReads, ...)

	}  # end of all sampleIDs for this possibly paired end sample

	return()
}
