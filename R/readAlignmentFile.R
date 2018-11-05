# readAlignmentFile.R

# efficient reading of aligments into a data.frame

`readAlignmentFile` <- function( filein, type="alignID", sampleID="", 
				checkLineCount=TRUE, verbose=TRUE) {

	fileToUse <- allowCompressedFileName( filein)

	nLines <- -1
	if ( checkLineCount) {
		nLines <- getFileLineCount( fileToUse, sampleID=sampleID)
		if ( nLines < 1) {
			warning( "readAlignmentFile:  unable to read alignments.")
			return( data.frame())
		}
	}

	# catch compressed files
	if ( filein != fileToUse) {
		con1 <- openCompressedFile( fileToUse, open="r")
		conIn <- file( "")
		writeLines( readLines(con1), con=conIn)
		close( con1)
	} else {
		conIn <- file( fileToUse, open="r")
	}

	# set up based on type of file
	myColClasses <- myColNames <- NA
	switch( type,
		"bowtie" = {
			myColClasses <- c( "character", "character", "character", "integer", "character", 
					"character", "integer", "character")
			myColNames <- BOWTIE_COLUMNS
			},
		"alignID" = {
			myColClasses <- c( "character", "character", "character", "integer", "character", 
					"character", "character")
			myColNames <- ALIGNID_COLUMNS
			},

		"SAM"= {}
	)
	#cat( "\nColumns: ", myColClasses, myColNames)

	# scan for existence of header line
	hasHeader <- FALSE
	if ( type == "alignID" ) {
		tmp <- readLines( con=conIn, n=1)
		hasHeader <- ( base::substr( tmp, 1,7) == "READ_ID")
		pushBack( tmp, conIn)
	}
	if ( type == "SAM" ) {
		repeat {
			tmp <- readLines( con=conIn, n=1)
			if ( base::substr( tmp, 1,1) != "@") {
				pushBack( tmp, conIn)
				break
			}
		}
	}

	# read in the alignment file
	if ( type != "SAM") {
		alignDF <- read.delim( file=conIn, header=hasHeader, colClasses=myColClasses, col.names=myColNames,
				comment.char="", quote="", stringsAsFactors=FALSE, nrows=nLines)
	} else {
		txt <- readLines( con=conIn)
		if (verbose) cat( "\nconverting SAM format..")
		alignDF <- Sam2Bowtie( txt)
	}
	nIn <- nrow(alignDF)

	if( verbose) {
		cat( "\nLoaded file: \t", fileToUse, "\nN_Alignments: \t", prettyNum( nIn, big.mark=","),"\n")
	}

	close( conIn)
	return( alignDF)
}


`test.readAlignmentFile` <- function() {

	# re-create a .align text file on disk
	tmpFile <- build.testAlignFile()
	zz <- readAlignmentFile( tmpFile, type="bowtie", verbose=FALSE)
	checkEquals( dim(zz), c(500,8))
	checkEquals( colnames(zz), BOWTIE_COLUMNS)
	remove.testFile( tmpFile)

	# the 'assigned' alignment after adding GeneIDs
	tmpFile <- build.testAlignIDFile()
	zz <- readAlignmentFile( tmpFile, type="alignID", verbose=FALSE)
	checkEquals( dim(zz), c(999,7))
	checkEquals( colnames(zz), ALIGNID_COLUMNS)
	checkEquals( typeof( zz$START), "integer")
	remove.testFile( tmpFile)
}


`setupBufferedReadAlignmentFile` <- function( filein, type="alignID", sampleID="", 
					bufferSize=1000000, verbose=TRUE) {

	bufferedReaderBufferSize <<- bufferSize
	bufferedReaderType <<- type
	fileToUse <- allowCompressedFileName( filein)

	nLines <- getFileLineCount( fileToUse, sampleID, verbose=verbose)
	if ( nLines == 0) {
		warning( "bufferedReadAlignmentFile:  unable to read alignments.")
		return( 0)
	}
	nBuffers <- ceiling( nLines / bufferedReaderBufferSize)

	if( verbose) {
		cat( "\nReading file: \t", fileToUse, "\nN_Alignments: \t", prettyNum( nLines, big.mark=","),
			"\nN_Buffers:    \t", nBuffers, "\n")
	}

	# catch compressed files
	if ( regexpr( "\\.gz$", fileToUse) > 0) {
		con1 <- gzfile( fileToUse, open="r")
		conIn <- file( "")
		writeLines( readLines(con1), con=conIn)
		close( con1)
	} else {
		conIn <- file( fileToUse, open="r")
	}
	bufferedReaderCon <<- conIn

	# set up based on type of file
	myColClasses <- myColNames <- NA
	switch( type,
		"bowtie" = {
			myColClasses <- c( "character", "character", "character", "integer", "character", 
					"character", "integer", "character")
			myColNames <- BOWTIE_COLUMNS
			},
		"alignID" = {
			myColClasses <- c( "character", "character", "character", "integer", "character", 
					"character", "character")
			myColNames <- ALIGNID_COLUMNS
			}
	)
	bufferedReaderColClasses <<- myColClasses
	bufferedReaderColNames <<- myColNames

	# scan for existence of header line
	tmp <- readLines( con=conIn, n=1)
	hasHeader <- ( base::substr( tmp, 1,7) == "READ_ID")
	pushBack( tmp, conIn)
	bufferedReaderHeader <<- hasHeader
	if ( type == "SAM" ) {
		repeat {
			tmp <- readLines( con=conIn, n=1)
			if ( base::substr( tmp, 1,1) != "@") {
				pushBack( tmp, conIn)
				break
			}
		}
		bufferedReaderHeader <<- FALSE
	}

	# we are all ready to read
	return( nBuffers)
}


`readAlignmentBuffer` <- function( verbose=FALSE) {

	# read in the alignment file
	BUF_SIZE <- bufferedReaderBufferSize
	if ( verbose) cat( "  bufferRead..")
	if ( bufferedReaderType != "SAM") {
		alignDF <- read.delim( file=bufferedReaderCon, header=bufferedReaderHeader, 
				colClasses=bufferedReaderColClasses, 
				col.names=bufferedReaderColNames, comment.char="", quote="", 
				stringsAsFactors=FALSE, nrows=BUF_SIZE)
	} else {
		txt <- readLines( con=bufferedReaderCon, n=BUF_SIZE)
		if ( verbose) cat( "  convertSAM..")
		alignDF <- Sam2Bowtie( txt)
	}
	bufferedReaderHeader <<- FALSE
	return( alignDF)
}


`cleanupBufferedReadAlignmentFile` <- function() {

	close( bufferedReaderCon)
	rm( bufferedReaderCon, bufferedReaderColClasses, bufferedReaderColNames, bufferedReaderHeader, 
		bufferedReaderType, bufferedReaderBufferSize, envir=.GlobalEnv)
}
