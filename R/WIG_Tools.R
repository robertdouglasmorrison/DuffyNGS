# WIG_Tools.R

# assorted Wiggle object functions and constants


# set up default constants
`WIG.defaults` <- function() {
	assign( "RecentSamples", vector(mode="character"), envir=WIG_Env)
	assign( "RecentSeqIDs", vector(mode="character"), envir=WIG_Env)
	assign( "RecentWiggleChunks", vector(mode="list"), envir=WIG_Env)
}


# get the stored read counts in the entire data structure
`WIG_getTotalRawReads` <- function( WIG) {

	out <- list( "Unique"=WIG$Info$UniqueReads, "Multi"=WIG$Info$TotalReads)
	return( out)
}

# special form for doing any RPKM calls, to not excessively inflate very low count samples...
`WIG_getTotalReadsForRPKM` <- function( WIG, minReadsPerSample=100000) {

	minUniques <- max( WIG$Info$UniqueReads, minReadsPerSample)
	minMultis <- max( WIG$Info$TotalReads, minReadsPerSample)
	out <- list( "Unique"=minUniques, "Multi"=minMultis)
	return( out)
}


`WIG_getWigglesOneSeq` <- function( WIG, seqID, verbose=FALSE) {

	# in 'combined' mode, the wiggle chunk is local to the object
	if ( seqID %in% names(WIG)) return( WIG[[ seqID]])

	if ( ! (seqID %in% rownames( WIG$Info$RawReads))) {
		cat( "\nSeqID not in WIG object: ", seqID)
		cat( "\nFound: ", rownames( WIG$Info$RawReads), "\n") 
		return(NULL)
	}

	# if not already found, this must be a separate storage WIG object
	if ( (! "StorageMode" %in% names(WIG)) || (WIG$StorageMode != "separate")) {
		cat( "\nInvalid WIG object:   cannot find storage mode")
		return(NULL)
	}

	# see if recently 'gotten'
	sampleID <- WIG$Info$SampleID
	recentSamples <- WIG_Env[[ "RecentSamples"]]
	recentSeqIDs <- WIG_Env[[ "RecentSeqIDs"]]

	# we are storing 1 wiggle track (chromosome) for each sample
	where <- match( sampleID, recentSamples, nomatch=0)
	if ( where > 0) {
		if ( seqID == recentSeqIDs[where]) {
			return(  WIG_Env[[ "RecentWiggleChunks"]][[where]])
		}
	}

	# load from disk, is called 'wiggleChunk'
	subwigPath <- WIG$SubWigFolder
	thisfile <- file.path( subwigPath, paste( sampleID, seqID, "WIG.rda", sep="."))
	if ( ! file.exists( thisfile)) {
		cat( "\nFailed to find WIG file for: ", seqID, "\nTried: ", thisfile)
		return(NULL)
	}
	if (verbose) cat( "  load", sub( ".WIG.rda", "", basename(thisfile)))
	load( thisfile)

	# save a copy 
	if ( where == 0) {
		where <- length( recentSamples) + 1
		WIG_Env[[ "RecentSamples"]][ where] <- sampleID
	}
	WIG_Env[[ "RecentSeqIDs"]][where] <- seqID
	WIG_Env[[ "RecentWiggleChunks"]][[where]] <- wiggleChunk

	return( wiggleChunk)
}


`WIG_updateWigglesOneSeq` <- function( WIG, seqID, wiggleChunk, errorNotFound=TRUE) {

	if ( ! (seqID %in% rownames( WIG$Info$RawReads))) {
		cat( "\nSeqID not in WIG object: ", seqID)
		cat( "\nFound: ", rownames( WIG$Info$RawReads), "\n") 
		cat( "\nWIG Update failed..")
		return(NULL)
	}

	# in 'combined' mode, the wiggle chunk is local to the object
	if ( seqID %in% names(WIG)) {
		WIG[[ seqID]] <- wiggleChunk
		sampleID <- WIG$Info$SampleID
		path <- sub( "align.+", "", WIG$Info$FileName[1])
		prefix <- getCurrentSpeciesFilePrefix()
		thisfile <- file.path( path, "wig", paste( sampleID, prefix, "WIG.rda", sep="."))
		wiggles <- WIG
		cat( "\rRe-writing WIG file: ", seqID, thisfile)
		save( wiggles, file=thisfile)
		return( WIG)
	}

	# find in the sparse files
	if ( (! "StorageMode" %in% names(WIG)) || (WIG$StorageMode != "separate")) {
		cat( "\nInvalid WIG object:   cannot find storage mode")
		cat( "\nWIG Update failed..")
		return(NULL)
	}

	# we are storing 1 wiggle track (chromosome) for each sample
	sampleID <- WIG$Info$SampleID
	subwigPath <- WIG$SubWigFolder
	thisfile <- file.path( subwigPath, paste( sampleID, seqID, "WIG.rda", sep="."))
	if ( ! file.exists( thisfile)) {
		if ( errorNotFound) {
			cat( "\nFailed to find WIG file for: ", seqID, "\nTried: ", thisfile)
			cat( "\nWIG Update failed..")
			return(NULL)
		}
		if ( ! file.exists( subwigPath)) dir.create( subwigPath, recursive=TRUE)
	}
	save( wiggleChunk, file=thisfile)
	recentSamples <- WIG_Env[[ "RecentSamples"]]
	recentSeqIDs <- WIG_Env[[ "RecentSeqIDs"]]
	where <- match( sampleID, recentSamples, nomatch=0)
	if ( where > 0) {
		WIG_Env[[ "RecentSeqIDs"]][where] <- seqID
		WIG_Env[[ "RecentWiggleChunks"]][[where]] <- wiggleChunk
	}
	return( WIG)
}


# add up number of reads in a region of one chromosomes wiggles object.
`WIG_getReadCountsInRegion` <- function( wiggleChunk, start, stop, readLength=32)  {

	# given a seqID's chunk from a wiggles object, get the subset of base depths for
	# the given range and turn into a read count
	N <- length( wiggleChunk)
	out <- rep( 0, times=N)

	for ( i in 1:N) {

		totalBases <- sum.baseDepthTableSubset( wiggleChunk[[i]], start, stop)
		out[i] <- totalBases / readLength
	}

	names( out) <- names( wiggleChunk)
	return( as.list.default(out))
}


`writeWIGfiles` <- function( WIG, sampleID, seqID=NULL, path=".", mode=c("combined", "separate"), 
				chromo.prefix="chr", chromo.name=NULL) {


		mode <- match.arg( mode)
		wigInfo <- WIG$Info
		species <- wigInfo$Species
		setCurrentSpecies( species)
		prefix <- getCurrentSpeciesFilePrefix()
		myFileName <- paste( sampleID, prefix, sep=".")
		seqSet <- getCurrentSeqMap()$SEQ_ID

		myStrands <- c( "Plus", "Minus")
		mySigns <- c( "(+)", "(-)")
		myColors <- c( "color=0,0,255", "color=255,0,0")

		# create the file for this set of wiggles
		if ( mode == "combined") {
			fileOut <- paste( myFileName, "wig", sep=".")
			fileOut <- file.path( path, fileOut)
			conOut <- file( fileOut, open="w")
		}

		for( i in 1:2) {
			strand <- myStrands[i]
			sign <- mySigns[i]
			colr <- myColors[i]

			if ( mode == "separate") {
				fileOut <- paste( myFileName, strand, "wig", sep=".")
				fileOut <- file.path( path, fileOut)
				conOut <- file( fileOut, open="w")
			}

			headerExtras <- "gridDefault=on yLineOnOff=on visibility=full "
			headerName <- paste( "name=\"", sampleID, "_", strand, "\"", sep="")
			headerDesc <- paste( "description=\"", sampleID, " Read Coverage ", sign, " Strand\"", sep="")
			cat( "\n", species, "\tRead Coverage: ", sign)

			# start it with the required header
			trackDefLine <- paste( "track type=wiggle_0", headerName, headerDesc, colr, headerExtras, sep=" ")
			writeLines( trackDefLine, con=conOut)
	
			nWig <- 0
			if ( ! is.null(seqID)) seqSet <- seqID
			for( s in seqSet) {
				#if ( names( wiggles)[iseq] != thisSeq) stop( "WIG object naming error")
				thisWIG <- WIG_getWigglesOneSeq( WIG, s)
				# get the strand we want
				thisWIG <- thisWIG[[ strand]]
				if ( nrow( thisWIG) < 1) next

				# the seqID we write out is trimmed down a bit
				seqIDout <- sub( base::paste( "^", species, sep=""), "", s)
				seqIDout <- sub( "^_", "", seqIDout)
				seqIDout <- sub( "^0", "", seqIDout)
				if ( ! is.na( chromo.prefix)) {
					if ( base::substr( seqIDout, 1,nchar(chromo.prefix)) != chromo.prefix) {
						seqIDout <- base::paste( chromo.prefix, seqIDout, sep="")
					}
				}
				if ( ! is.null( chromo.name)) {
					seqIDout <- chromo.name
				}
				declareLine <- paste( "variableStep chromo=", seqIDout, sep="")
				writeLines( declareLine, con=conOut)
				cat( "\nDeclareLine: ", declareLine)

				# convert this wiggle set to the form needed
				#1.  if the base depth is less than 1, drop it
				drops <- which( round( thisWIG$DEPTH) < 1)
				if ( length(drops)) thisWIG <- thisWIG[ -drops, ]

				#2.  convert to locations & depths
				ans <- baseDepthTableToVector( thisWIG)
				loc <- as.integer( names(ans))
				deep <- as.integer( round( ans))

				#3.  tiny chance that some locations are not in increasing order...
				bad <- which( diff( loc) < 1)
				if ( length(bad)) {
					drops <- bad + 1
					cat( "\nTrap bad WIG loci: ", strand, length(bad), "|", loc[drops], "|", deep[drops])
					loc <- loc[ -drops]
					deep <- deep[ -drops]
				}
				sml <- base::paste( loc, deep, sep="\t")
				writeLines( sml, con=conOut)

				nWig <- nWig + nrow( thisWIG)
			}
			if (mode == "separate") close( conOut)
		}
		if (mode == "combined") close( conOut)
		cat( "\nWrote WIG file: ", fileOut, "\nN_Wiggles: ", nWig)
		return()
}
