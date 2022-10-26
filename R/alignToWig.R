# alignToWig.R

# turn alignments into seperate .WIG file for each species, chromosome, and strand
# doing it a buffer at a time, adding wiggles to small buffer files for each chromo/strand/etc.
# using new 'baseDepthVector' method instead of looping...


`alignToWig` <- function( filein, sampleID, reload=TRUE, annotationFile="Annotation.txt", optionsFile="Options.txt", 
			readSense="sense", dataType="RNA-seq", results.path=NULL, maxReads=NULL, readBufferSize=1000000, 
			verbose=TRUE) {

	if ( verbose) {
		cat( verboseOutputDivider)
		cat( "\n\nConverting alignments to Wiggles for:  ", filein, 
			"\nRead Sense:      ", readSense, "\n")
	}

	DataType <- dataType
	FileNameSet <- vector()

	# let's better account for wildly different read lengths
	RawReadLength <- 0
	RawReadCount <- 0
	isSPLICEfile <- (regexpr( "splic", filein) > 0)

	nbuf <- 0

	# set up to repeatedly call local functions
	optT <- readOptionsTable( optionsFile)
	if( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	wigPath <- file.path( results.path, "wig")
	if ( ! file.exists( wigPath)) dir.create( wigPath, recursive=T, showWarnings=F)

	# create temp storage
	tmp.path <- getOptionValue( optT, "tmp.path", notfound="/tmp")
	tempWigPath <- tempFolder( sampleID, tmpRoot=tmp.path)
	on.exit( removeTempFolder( tempWigPath))

	# get the set of all species we expect to see
	saveSpecies <- getCurrentSpecies()
	saveTarget <- getCurrentTarget()
	target <- getAndSetTarget( optionsFile, sampleID=sampleID, annotationFile=annotationFile, verbose=F)
	mySpeciesSet <- getCurrentTargetSpecies()
	myReadSense <- tolower(readSense)

	# track total wiggle count to decide what format to use for storage
	totalWiggleCount <- 0
	storageMode <- "combined"
	subwigPath <- NULL

	# when we fill wiggles by buffer, we don't see every chromosome every time, so keep a global tally
	allSeqWiggleCnts <- vector()
	for (spec in mySpeciesSet) {
		setCurrentSpecies( spec)
		seqIDs <- getCurrentSeqMap()$SEQ_ID
		smlCnts <- rep( 0, times=length(seqIDs))
		names(smlCnts) <- seqIDs
		allSeqWiggleCnts <- c( allSeqWiggleCnts, smlCnts)
	}


	PASTE <- base::paste

	# several local functions...

	# count up the wiggles per species
	`tallyWiggleCountsBySpecies` <- function( cnts) {
		cnts <- unlist( cnts)
		countsOut <- vector( mode="list")
		for ( s in mySpeciesSet) {
			setCurrentSpecies(s)
			seqSet <- getCurrentSeqMap()$SEQ_ID
			mycount <- 0
			who <- which( names(cnts) %in% seqSet)
			if ( length( who) > 0) {
				#mycount <- sum( sapply( cnts[who], function(x) x[1]), na.rm=T)
				mycount <- sum( cnts[who], na.rm=T)
			}
			countsOut[[ s]] <- mycount
		}
		names( countsOut) <- mySpeciesSet
		return( countsOut)
	}


	# create the name for a single temp file that holds one strand of one chromo of one species
	`scratchFileName` <- function( sampleID, seq, strand) {
		return( file.path( tempWigPath, PASTE( sampleID, seq, strand, "tmp.rda", sep=".")))
	}


	# erase all scratch files we may need
	`eraseAllScratchFiles` <- function( sampleID) {
		for ( s in mySpeciesSet) {
			setCurrentSpecies(s)
			seqSet <- getCurrentSeqMap()$SEQ_ID
			for (seq in seqSet) {
				thisfile <- scratchFileName( sampleID, seq, "Plus")
				if ( file.exists( thisfile)) file.remove( thisfile)
				thisfile <- scratchFileName( sampleID, seq, "Minus")
				if ( file.exists( thisfile)) file.remove( thisfile)
				thisfile <- scratchFileName( sampleID, seq, "Plus.unique")
				if ( file.exists( thisfile)) file.remove( thisfile)
				thisfile <- scratchFileName( sampleID, seq, "Minus.unique")
				if ( file.exists( thisfile)) file.remove( thisfile)
			}
		}
	}



	# erase all final result .wig files we will make anew
	`eraseAllOutputWigFiles` <- function( sampleID) {
		for ( s in mySpeciesSet) {
			setCurrentSpecies(s)
			prefix <- getCurrentSpeciesFilePrefix()
			thisfile <- file.path( wigPath, PASTE( sampleID, prefix, "WIG.rda", sep="."))
			if ( file.exists( thisfile)) file.remove( thisfile)
			thisfile <- file.path( wigPath, PASTE( sampleID, prefix, "Plus.wig", sep="."))
			if ( file.exists( thisfile)) file.remove( thisfile)
			thisfile <- file.path( wigPath, PASTE( sampleID, prefix, "Minus.wig", sep="."))
			if ( file.exists( thisfile)) file.remove( thisfile)
			thisfile <- file.path( wigPath, PASTE( sampleID, prefix, "Plus.unique.wig", sep="."))
			if ( file.exists( thisfile)) file.remove( thisfile)
			thisfile <- file.path( wigPath, PASTE( sampleID, prefix, "Minus.unique.wig", sep="."))
			if ( file.exists( thisfile)) file.remove( thisfile)
			
			# also remove any seqID level .rda files
			subwigPath <- file.path( wigPath, sampleID)
			if ( ! file.exists( subwigPath)) next
			seqSet <- getCurrentSeqMap()$SEQ_ID
			for (seq in seqSet) {
				thisfile <- file.path( subwigPath, PASTE( sampleID, seq, "WIG.rda", sep="."))
				if ( file.exists( thisfile)) file.remove( thisfile)
			}
		}
	}


	# turn all the sequence tmp files into 1 WIG.rda object/file
	`compressWIGfromScratchFiles` <- function( sampleID) {
	
		# we will store the current best estimate of the average read length
		storedReadLength <- round( RawReadLength, digits=2)

		for ( s in mySpeciesSet) {
			setCurrentSpecies(s)
			prefix <- getCurrentSpeciesFilePrefix()
			seqSet <- getCurrentSeqMap()$SEQ_ID
			wiginfo <- vector( mode="list", length=length(seqSet))
			wigdepth <- vector( mode="list", length=length(seqSet))

			# use the number of wiggles to decide how to pack
			if ( totalWiggleCountPerSpecies[[ s]] > 5000000) {
				storageMode <- "separate"
				cat("\nWIG storage mode for ", s, ": ", storageMode)
				subwigPath <- file.path( wigPath, sampleID)
				if ( ! file.exists( subwigPath)) dir.create( subwigPath, recursive=T)
			} else {
				storageMode <- "combined"
				cat("\nWIG storage mode for ", s, ": ", storageMode)
				subwigPath <- ""
			}

			wiggles <- vector( mode="list")

			for (iseq in 1:length(seqSet)) {
				seq <- seqSet[iseq]
				thisfile <- scratchFileName( sampleID, seq, "Plus")
				if ( ! file.exists( thisfile)) {
					objP <- EMPTY_BASE_DEPTH_TABLE
					infoP <- 0
				} else {
					load( thisfile)
					objP <- if (savedAsTable) saveDF else baseDepthVectorToTable( saveDF)
					infoP <- saveNaligns
					if (length( saveFileSet) > length( FileNameSet)) {
						FileNameSet <- saveFileSet
					}
				}
				thisfile <- scratchFileName( sampleID, seq, "Minus")
				if ( ! file.exists( thisfile)) {
					objM <- EMPTY_BASE_DEPTH_TABLE
					infoM <- 0
				} else {
					load( thisfile)
					objM <- if (savedAsTable) saveDF else baseDepthVectorToTable( saveDF)
					infoM <- saveNaligns
					if (length( saveFileSet) > length( FileNameSet)) {
						FileNameSet <- saveFileSet
					}
				}
				thisfile <- scratchFileName( sampleID, seq, "Plus.unique")
				if ( ! file.exists( thisfile)) {
					objPU <- EMPTY_BASE_DEPTH_TABLE
					infoPU <- 0
				} else {
					load( thisfile)
					objPU <- if (savedAsTable) saveDF else baseDepthVectorToTable( saveDF)
					infoPU <- saveNaligns
					if (length( saveFileSet) > length( FileNameSet)) {
						FileNameSet <- saveFileSet
					}
				}
				thisfile <- scratchFileName( sampleID, seq, "Minus.unique")
				if ( ! file.exists( thisfile)) {
					objMU <- EMPTY_BASE_DEPTH_TABLE
					infoMU <- 0
				} else {
					load( thisfile)
					objMU <- if (savedAsTable) saveDF else baseDepthVectorToTable( saveDF)
					infoMU <- saveNaligns
					if (length( saveFileSet) > length( FileNameSet)) {
						FileNameSet <- saveFileSet
					}
				}

				wiggleChunk <- list( "Plus"=objP, "Minus"=objM, "PlusUnique"=objPU,
							"MinusUnique"=objMU)
				wiginfo[[ iseq]] <- list( "Plus"=infoP, "Minus"=infoM, "PlusUnique"=infoPU,
							"MinusUnique"=infoMU)
				if (DataType == "DNA-seq") {
					wigdepth[[ iseq]] <- list( "Plus"=median(objP$DEPTH), "Minus"=median(objM$DEPTH))
				}

				if ( storageMode == "combined") {
					wiggles[[ iseq]] <- wiggleChunk
				} else {
					thisFile <- file.path( subwigPath, PASTE( sampleID, seq, "WIG.rda", 
								sep="."))
					save( wiggleChunk, file=thisFile)
					# show some life...
					cat( ".")
				}
			}

			names(wiginfo) <- seqSet
			names(wigdepth) <- seqSet
			if ( storageMode == "combined") names(wiggles) <- seqSet
			wigspecies <- s

			# now that everyting is ready, repackage the number of alignments into one 'Info' object
			countsM <- matrix( base::unlist( wiginfo), nrow=4, ncol=length(seqSet))
			rownames(countsM) <- c( "Plus", "Minus", "PlusUnique", "MinusUnique")
			colnames(countsM) <- seqSet
			countsM <- t( countsM)
			if (DataType == "DNA-seq") {
				depthM <- matrix( base::unlist( wigdepth), nrow=2, ncol=length(seqSet))
				rownames(depthM) <- c( "Plus", "Minus")
				colnames(depthM) <- seqSet
				depthM <- t( depthM)
			}

			sums <- apply( countsM, MARGIN=2, sum, na.rm=T)
			totalReads <- sum( sums[1:2])
			totalUniques <- sum( sums[3:4])

			info <- list( "SampleID"=sampleID, "Species"=wigspecies, 
					"FileName"=FileNameSet, "RawReads"=countsM, 
					"ReadLength"=storedReadLength, "TotalReads"=totalReads,
					"UniqueReads"=totalUniques, "DataType"=DataType)
			if (DataType == "DNA-seq") {
				info <- list( "SampleID"=sampleID, "Species"=wigspecies, 
					"FileName"=FileNameSet, "RawReads"=countsM, 
					"ReadLength"=storedReadLength, "TotalReads"=totalReads,
					"UniqueReads"=totalUniques, "DataType"=DataType,
					"MedianDepth"=depthM)
			}
			wigFile <- PASTE( sampleID, prefix, "WIG.rda", sep=".")
			wigFile <- file.path( wigPath, wigFile)
			wiggles$Info <- info
			wiggles$StorageMode <- storageMode
			wiggles$SubWigFolder <- subwigPath
			save( wiggles, file=wigFile)
			cat( "\nWrote compressed WIG file:  ", wigFile)
		}
		return()
	}


	`expandWIGtoScratchFiles` <- function( sampleID) {

		fileNamesSeen <- vector()
		totalReadsSeen <- 0
		readLengthSeen <- 0
	
		for ( s in mySpeciesSet) {
			setCurrentSpecies(s)
			prefix <- getCurrentSpeciesFilePrefix()
			seqSet <- getCurrentSeqMap()$SEQ_ID
			wigFile <- PASTE( sampleID, prefix, "WIG.rda", sep=".")
			wigFile <- file.path( wigPath, wigFile)
			if ( ! file.exists( wigFile)) next
			load( file=wigFile)

			# grab what we need from the stored info header
			info <- wiggles$Info
			thisLength <- info$ReadLength
			thisCount <- info$TotalReads
			readLengthSeen <- max( readLengthSeen, thisLength)
			totalReadsSeen <- totalReadsSeen + thisCount
			saveFileSet <- info$FileName

			storageMode <- "combined"
			if ( "StorageMode" %in% names(wiggles)) {
				storageMode <- wiggles$StorageMode
				subwigPath <- wiggles$SubWigFolder
			}
			countsM <- info$RawReads

			# all are in table mode at this point
			savedAsTable <- TRUE


			expandOneWiggleChuckToScratchFile <- function( iseq) {

				seq <- seqSet[iseq]
				if (storageMode == "separate") {
					thisFile <- file.path( subwigPath, PASTE( sampleID, seq, "WIG.rda", sep="."))
					load( thisFile)
				} else {
					wiggleChunk <- wiggles[[iseq]]
				}

				thisfile <- scratchFileName( sampleID, seq, "Plus")
				saveDF <- wiggleChunk$Plus
				saveNaligns <- countsM[ iseq, "Plus"]
				if ( nrow( saveDF) < 1) {
					file.delete( thisfile)
				} else {
					save( saveDF, saveNaligns, saveFileSet, savedAsTable, file=thisfile)
				}

				thisfile <- scratchFileName( sampleID, seq, "Minus")
				saveDF <- wiggleChunk$Minus
				saveNaligns <- countsM[ iseq, "Minus"]
				if ( nrow( saveDF) < 1) {
					file.delete( thisfile)
				} else {
					save( saveDF, saveNaligns, saveFileSet, savedAsTable, file=thisfile)
				}

				thisfile <- scratchFileName( sampleID, seq, "Plus.unique")
				saveDF <- wiggleChunk$PlusUnique
				saveNaligns <- countsM[ iseq, "PlusUnique"]
				if ( nrow( saveDF) < 1) {
					file.delete( thisfile)
				} else {
					save( saveDF, saveNaligns, saveFileSet, savedAsTable, file=thisfile)
				}

				thisfile <- scratchFileName( sampleID, seq, "Minus.unique")
				saveDF <- wiggleChunk$MinusUnique
				saveNaligns <- countsM[ iseq, "MinusUnique"]
				if ( nrow( saveDF) < 1) {
					file.delete( thisfile)
				} else {
					save( saveDF, saveNaligns, saveFileSet, savedAsTable, file=thisfile)
				}
			}
			multicore.lapply( 1:length(seqSet), FUN=expandOneWiggleChuckToScratchFile)

			# all species have the exact same list of file names, so only collect them once
			if ( s == mySpeciesSet[1]) {
				fileNamesSeen <- c( fileNamesSeen, info$FileName)
			}
		}

		# stuff the info about read length back into the global locations
		RawReadLength <<- readLengthSeen
		RawReadCount <<- totalReadsSeen

		return( fileNamesSeen)
	}


	`OneSeqOneStrand2wiggles` <- function( x) {

		# given a small subset of pointers to the rows that are one strand of one chromosome
		if ( length(x) < 1) return(0)

		thisSeqID <- bufDF$SEQ_ID[ x[1]]

		if ( myReadSense == "sense") {
			thisStrand <- if ( bufDF$STRAND[ x[1]] == "-") "Plus" else "Minus"
		} else {
			thisStrand <- if ( bufDF$STRAND[ x[1]] == "+") "Plus" else "Minus"
		}

		# extract just the data we will need
		startSet <- bufDF$START[ x]
		stopSet <- bufDF$STOP[ x]
		wtSet <- bufDF$DEPTH[ x]

		# make a baseDepth vector from the new data
		thisBaseVec <- baseDepthVector( startSet, stopSet, wtSet)
	
		# get the chunk we already have for this seq/strand
		#saveDF <- EMPTY_BASE_DEPTH_TABLE
		saveDF <- EMPTY_BASE_DEPTH_VECTOR
		saveNaligns <- 0
		saveFile <- scratchFileName( sampleID, thisSeqID, thisStrand)
		savedAsTable <- FALSE
		if ( file.exists( saveFile)) load( saveFile)
		if (savedAsTable) saveDF <- baseDepthTableToVector( saveDF)
		savedAsTable <- FALSE

		# merge them
		if ( length( saveDF) > 1) {
			bigBaseVec <- mergeIntegerTables( thisBaseVec, saveDF)
		} else {
			bigBaseVec <- thisBaseVec
		}

		# leave this result in the 1D table form
		#saveDF <- baseDepthVectorToTable( bigBaseVec)
		saveDF <- bigBaseVec
		saveNaligns <- round( sum( wtSet)) + saveNaligns
		saveFileSet <- FileNameSet
		save( saveDF, saveNaligns, saveFileSet, savedAsTable, file=saveFile)

		#totalWiggles <<- totalWiggles + length(saveDF)
		curWiggleCount <- length(saveDF)

		# some data types do not make a unique vs multi distinction
		if ( DataType == "ChIP-seq") return( curWiggleCount)
		# since RIP-seq is at it's heart RNA-seq type reads, let's treat it more lake RNA than ChIP
		#if ( DataType == "RIP-seq") return( curWiggleCount)

		# extract just the unique read data we will need
		uniqs <- x[ bufDF$DEPTH[x] == 1]
		startSet <- bufDF$START[ uniqs]
		stopSet <- bufDF$STOP[ uniqs]
		wtSet <- bufDF$DEPTH[ uniqs]
		saveNaligns <- round( sum( wtSet))
		savedAsTable <- FALSE
		thisBaseVec <- baseDepthVector( startSet, stopSet, wtSet)
		saveFile <- scratchFileName( sampleID, thisSeqID, PASTE( thisStrand, "unique", sep="."))
		#if ( ! reload) {

			# get the chunk we already have for this unique seq/strand
			saveDF <- EMPTY_BASE_DEPTH_VECTOR
			saveNaligns <- 0
			savedAsTable <- FALSE
			if ( file.exists( saveFile)) load( saveFile)
			if (savedAsTable) saveDF <- baseDepthTableToVector( saveDF)
			savedAsTable <- FALSE

			# merge them
			if ( length( saveDF) > 1) {
				bigBaseVec <- mergeIntegerTables( thisBaseVec, saveDF)
			} else {
				bigBaseVec <- thisBaseVec
			}

			# leave this result in the 1D vector form
			saveDF <- bigBaseVec
			saveNaligns <- round( sum( wtSet)) + saveNaligns
			saveFileSet <- FileNameSet
		#} # end of what may be different as 'uniques'..

		save( saveDF, saveNaligns, saveFileSet, savedAsTable, file=saveFile)

		# done with the unique-only copy
		return( curWiggleCount)
	}


	`writeWIGfiles` <- function( sampleID, species) {


		setCurrentSpecies( species)
		prefix <- getCurrentSpeciesFilePrefix()
		myFileName <- paste( sampleID, prefix, sep=".")
		seqSet <- getCurrentSeqMap()$SEQ_ID

		for( strand in c( "Plus", "Minus")) {

		headerExtras <- "color=100,50,150 gridDefault=on yLineOnOff=on visibility=full maxheightPixels=40:40:12"

		if ( strand == "Plus") {
			headerName <- paste( "name=\"", myFileName, "_plus\"", sep="")
			headerDesc <- paste( "description=\"", myFileName, "Read Coverage (+) Orientation\"", sep="")
			cat( "\n\n", species, "\tRead Coverage (+) Orientation")
		} else {
			headerName <- paste( "name=\"", myFileName, "_minus\"", sep="")
			headerDesc <- paste( "description=\"", myFileName, "Read Coverage (-) Orientation\"", sep="")
			cat( "\n\n", species, "\tRead Coverage (-) Orientation")
		}

		# create the file for this set of wiggles
		fileOut <- paste( myFileName, strand, "wig", sep=".")
		fileOut <- file.path( results.path, "wig", fileOut)
		conOut <- file( fileOut, open="w")

		# start it with the required header
		trackDefLine <- paste( "track type=wiggle_0", headerName, headerDesc, headerExtras, sep=" ")
		writeLines( trackDefLine, con=conOut)
	
		# get that species wiggle object, its always called 'wiggles'
		load( file.path( results.path, "wig", paste( sampleID, prefix, "WIG.rda", sep=".")))

		nWig <- 0
		for( iseq in 1:length( seqSet)) {
			thisSeq <- seqSet[iseq]
			if ( names( wiggles)[iseq] != thisSeq) stop( "WIG object naming error")
			thisWIG <- wiggles[[ iseq]]
			# get the strand we want
			thisWIG <- thisWIG[[ strand]]
			if ( nrow( thisWIG) < 1) next

			# the seqID we write out is trimmed down a bit
			seqIDout <- sub( PASTE( "^", species, "_", sep=""), "", thisSeq)
			seqIDout <- sub( "^0", "", seqIDout)
			if ( base::substr( seqIDout, 1,3) != "chr") seqIDout <- PASTE( "chr", seqIDout, sep="")

			# convert this wiggle set from actual 'base-1' positions to 'base-0, half-open'
			seqOut <- rep( seqIDout, times=nrow( thisWIG))
			startsOut <- thisWIG$START - 1
			endsOut <- thisWIG$STOP
			dataOut <- thisWIG$DEPTH
			sml <- PASTE( seqOut, as.integer(startsOut), as.integer(endsOut), dataOut, sep="\t")
			writeLines( sml, con=conOut)

			nWig <- nWig + nrow( thisWIG)
		}
		close( conOut)
		cat( "\nWrote WIG file: ", fileOut, "\nN_Wiggles: ", nWig)

		} # end of 'for strand...'

		return()
	}
	# end of all local functions...



	# we either start fresh, or add to an existing...
	if ( reload) {
		cat( "\nCreating new WIG files for: ", sampleID)
		# pre-erase any files we will make and all the temp files we will need...
		eraseAllOutputWigFiles( sampleID)
		eraseAllScratchFiles( sampleID)
		FileNameSet <- filein
	} else {
		cat( "\nAdding to existing WIG files for: ", sampleID)
		curfiles <- expandWIGtoScratchFiles( sampleID)
		FileNameSet <- c( curfiles, filein)
	}

	# try to gracefully handle a missing file, as in the case if zero reads ever aligned
	if ( ! file.exists( filein)) {
		totalWiggleCountPerSpecies <- tallyWiggleCountsBySpecies( allSeqWiggleCnts)
		compressWIGfromScratchFiles( sampleID)
		eraseAllScratchFiles( sampleID)
		setCurrentSpecies( saveSpecies)
		return( 0)
	}

	# open that BAM file, and kow the dictionary
	conBam <- bamReader( filein)
	refData <- getRefData( conBam)

	# we will make wiggles one buffer at a time
	allAligns <- 0
	hasMore <- TRUE
	repeat {
		if ( ! hasMore) break
		if (verbose) cat( "\nReadBAM..")
		chunk <- getNextChunk( conBam, n=readBufferSize, alignedOnly=TRUE)
		nNow <- size(chunk)
		if ( nNow < 1) break
		if ( nNow < readBufferSize) hasMore <- FALSE
		nbuf <- nbuf + 1
		nAlign <- nNow
		allAligns <- allAligns + nAlign
		cat( " N:", formatC( allAligns, big.mark=",", format="d"))

		# accumulate the bits we need
		wtStrs <- getTag( chunk, BAMTAG_READWEIGHT)
		weights <- as.numeric( wtStrs)
		weights[ is.na(weights)] <- 1
		refIDs <- refID( chunk)
		seqIDs <- refID2seqID( refIDs, refData=refData)
		strands <- ifelse( reverseStrand(chunk), "-", "+")
		starts <- position( chunk)
		readLens <- alignLength( chunk)

		# there MIGHT be soft clipping from local mode alignment
		readLens <- softClipAlignLength( readLens, cigarString(chunk))
		stops <- starts + readLens - 1

		# if mate pairs, then the second mate should get its strand swapped
		isMate2 <- secondInPair(chunk)
		doStrandSwap <- which( isMate2)
		if ( length( doStrandSwap)) {
			strands[ doStrandSwap] <- ifelse( strands[doStrandSwap] == "+", "-", "+")
		}

		# it's possible we have some zero length reads, mostly from splice ends
		# find and drop them now
		isZeroLen <- which( as.numeric(starts) > as.numeric(stops))
		if ( length(isZeroLen)) {
			starts <- starts[ -isZeroLen]
			stops <- stops[ -isZeroLen]
			weights <- weights[ -isZeroLen]
			strands <- strands[ -isZeroLen]
			readLens <- readLens[ -isZeroLen]
			refIDs <- refIDs[ -isZeroLen]
			seqIDs <- seqIDs[ -isZeroLen]
			nAlign <- length(starts)
		}

		# try to more accurately track the average length of reads
		# note that splices are broken into to shorter pieces
		thisReadLength <- mean( readLens)
		if (isSPLICEfile) thisReadLength <- thisReadLength * 2
		thisReadCount <- round( sum( weights))
		if ( RawReadCount == 0) {
			RawReadLength <- thisReadLength
			RawReadCount <- thisReadCount
		} else {
			RawReadLength <- ((RawReadLength * RawReadCount) + (thisReadLength * thisReadCount)) / 
						(RawReadCount + thisReadCount)
			RawReadCount <- RawReadCount + thisReadCount
		}
		# factorize to do each chromo strand as a group
		seqFac <- factor( refIDs)
		strandFac <- factor( strands)
		# remove what we're done with
		rm( refIDs, readLens)
		gc(); gc()

		# make the set of details we need
		bufDF <- data.frame( seqIDs, strands, as.numeric(starts), as.numeric(stops), 
					as.numeric(weights), stringsAsFactors=FALSE)
		colnames( bufDF) <- c( "SEQ_ID", "STRAND", "START", "STOP", "DEPTH")
		rownames( bufDF) <- 1:nAlign

		# dispatch each 
		cat( "  add wiggles..")
		cntsRaw <- multicore.tapply( 1:nrow(bufDF), INDEX=list( strandFac, seqFac), FUN=OneSeqOneStrand2wiggles)
		cntNames <- tapply( 1:nrow(bufDF), INDEX=list( strandFac, seqFac), FUN=function(x) bufDF$SEQ_ID[x[1]])
		names(cntsRaw) <- cntNames
		# the multicore 'may' throw an error, so catch and force all the counts to numeric
		notNull <- sapply( cntsRaw, FUN=is.null)
		cntsRaw <- unlist( sapply( cntsRaw, FUN=function(x) x[1]))
		cnts <- as.numeric( cntsRaw)
		if ( any( is.na(cnts))) {
			good <- which( ! is.na(cnts))
			#cat( "\n  Good Count returned: ", cntNames[good], "|", cntsRaw[good], "\n")
			bad <- which( is.na(cnts))
			#cat( "\n  Bad Count returned:  ", cntNames[bad], "|", cntsRaw[bad], "\n")
		}
		names(cnts) <- cntNames[notNull]

		# update the counts on the chromosomes we just saw
		whereCnt <- match( names( allSeqWiggleCnts), cntNames, nomatch=0)
		allSeqWiggleCnts[ whereCnt > 0] <- cnts[ whereCnt]

		# totals are based on everyone, not just this buffer
		totalWiggleCountPerSpecies <- tallyWiggleCountsBySpecies( allSeqWiggleCnts)
		totalWiggleCount <- sum( as.numeric( unlist( totalWiggleCountPerSpecies)))
		cat( "   N_Wiggles: ", formatC( totalWiggleCount, format="d", big.mark=","))

		if ( ! is.null( maxReads)) {
			if ( allAligns >= maxReads) break
		}
	}
	bamClose( conBam)

	if ( nbuf < 1) {
		warning( paste("Failed to read any alignments from file: ", filein))
		return( 0)
	}


	# after all buffers full, save one compressed WIG buffer object
	compressWIGfromScratchFiles( sampleID)

	# now free to erase those temp files.
	eraseAllScratchFiles( sampleID)

	cat( "\n")

	# reset now that we're done
	setCurrentSpecies( saveSpecies)
	return( nbuf)
}
