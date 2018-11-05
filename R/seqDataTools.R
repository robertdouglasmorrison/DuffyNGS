# seqDataTools.R


`detectFastqReadFormat` <- function( filein) {

	# get ready to read the fastq file in chuncks...
	filein <- allowCompressedFileName( filein)
	if ( ! file.exists( filein)) stop( paste("Can't find input file: ", filein))
	conIn <- openCompressedFile( filein, open="r")

	# there is often garbage at the front...
	# so step past it...
	tmp <- readLines( con=conIn, n=10000)

	# now use the last few
	nlines <- length(tmp)
	tmp <- tmp[ max( c(nlines-40+1, 1)):nlines]
	nlines <- length(tmp)
	close( conIn)

	readids <- tmp[ seq( 1, nlines, by=4)]
	readseqs <- tmp[ seq( 2, nlines, by=4)]
	scores <- tmp[ seq( 4, nlines, by=4)]

	# verify bases...
	mybases <- base::unlist( strsplit( readseqs, split=""))
	whobad <- which( !( mybases %in% c("A","C","G","T","N",".")))
	if ( length(whobad) > 0) {
		stop( paste( "Bad fastq file:  non-DNA bases detected:\nFile: ", filein, 
				"\nBad Bases: ", mybases[whobad]))
	}

	# see what separators are used in the read IDs
	readIDtype <- detectReadIDformat( readids)

	# detect the score encryption
	try33 <- phredScoreStringToInt( scores, "Phred33")
	try64 <- phredScoreStringToInt( scores, "Phred64")
	# none should be negative
	scoreType <- "Phred64"
	nTooSmall33 <- sum( try33 < 1)
	nTooSmall64 <- sum( try64 < 1)
	if ( nTooSmall33 < nTooSmall64) {
		scoreType <- "Phred33"
	}

	out <- list( "readIDtype"=readIDtype, "scoreType"=scoreType)
	return( out)
}


`detectReadIDformat` <- function( readids) {

	N <- length( readids)

	# see what separators are used in the read IDs
	colons <- base::unlist( gregexpr( ":", readids, fixed=T))
	underscores <- base::unlist( gregexpr( "_", readids, fixed=T))
	periods <- base::unlist( gregexpr( ".", readids, fixed=T))
	idealSolexa1.3 <- length(colons) / (4 * N)
	idealRosetta <- length(underscores) / (3 * N)
	idealSolexa1.8 <- length(colons) / (9 * N)
	idealUnusual1 <- length(periods) / (5 * N)
	
	bestGuess <- c( idealSolexa1.3, idealRosetta, idealSolexa1.8, idealUnusual1)
	names(bestGuess) <- c( "Solexa_1.3", "Rosetta", "Solexa_1.8", "Unusual1")
	best <- which.min( abs( bestGuess - 1))
	readIDtype <- names(bestGuess)[best]

	return( readIDtype)
}


`extractReadIDlaneNumber` <- function( txt, isFastq=TRUE, readIDtype="Rosetta") {

	# the read ID is " @<lane>_<tile>_<x>_<y> "
	out <- rep( 0, times=length( txt))

	if ( readIDtype == "Rosetta") {
		seps <- gregexpr( "_", txt[1], fixed=TRUE)[[1]]
		nsep <- length(seps) 
		sep <- seps[ nsep - 2]
		out <- as.integer( base::substr( txt, (sep-1), (sep-1)))
	} else if ( readIDtype == "Solexa_1.3") {
		# in derisi IDs there is a ":n:m:x:y" field...
		seps <- gregexpr( ":", txt[1], fixed=TRUE)[[1]]
		sep <- if (length(seps) > 0) seps[1] else 0
		if (length(seps) > 4) sep <- seps[ length(seps) - 3]
		out <- ifelse( sep > 0, as.integer( base::substr( txt, (sep+1), (sep+1))), 0)
	} else if ( readIDtype == "Solexa_1.8") {
		# in newest IDs there is a "instrumentID:runID:flowCellID:n:m:x:y a:b:c:" field...
		# in derisi IDs there is a ":n:m:x:y" field...
		seps <- gregexpr( ":", txt[1], fixed=TRUE)[[1]]
		out <- ifelse( length(seps) > 3, as.integer( base::substr( txt, (seps[3]+1), (seps[4]-1))), 0)
	} else if ( readIDtype == "Unusual1") {
		# in IDs there is a ".n.m.x.y.readnumber" field...
		seps <- gregexpr( ".", txt[1], fixed=TRUE)[[1]]
		out <- ifelse( length(seps) > 1, as.integer( base::substr( txt, (seps[1]+1), (seps[2]-1))), 0)
	} else {
		stop( paste( "extractReadIDlaneNumber:  unknowm readID type: ", readIDtype, 
				"Known:  Rosetta, Solexa_1.3, Solexa_1.8"))
	}

	if ( ! all( out > 0)) {
		cat("\nExtracting ReadID Lane Number failed:   readID sample=\n")
		print( head( txt))
	}

	return( out)
}


`extractReadIDterms` <- function( txt, isFastq=TRUE, readIDtype="Rosetta") {

	# the read ID is " @<lane>_<tile>_<x>_<y> " for older Rosetta style
	# in DeRisi and newer Solexa, the  ID is a "@name:lane:tile:x:y" style
	# in newest Solexa, the  ID is a "@name:str::lane:tile:x:y a:b:c:" style
	outLane <- outTile <- outX <- outY <- rep( 0, times=length( txt))

	strIn <- txt

	skipChar <- 0
	if ( readIDtype == "Rosetta") {
		if ( isFastq) skipChar <- 1
		splitMark <- "_"
		wantedTerms <- 1:4
		fixed <- TRUE
	} else if ( readIDtype == "Solexa_1.3") {
		splitMark <- ":"
		wantedTerms <- 2:5
		fixed <- TRUE
	} else if ( readIDtype == "Solexa_1.8") {
		splitMark <- ":| "
		wantedTerms <- 4:11
		fixed <- FALSE
	} else if ( readIDtype == "Unusual1") {
		splitMark <- "."
		wantedTerms <- 2:5
		fixed <- TRUE
	} else {
		stop( paste( "extractReadIDterms:  unknowm readID type: ", readIDtype, 
				"Known:  Rosetta, Solexa_1.3, Solexa_1.8"))
	}
	markPtrs <- gregexpr( splitMark, strIn[1], fixed=fixed)[[1]]
	nTermsEach <- length( markPtrs) + 1

	# solexa_1.8 does not always have a barcode
	hasBarcode <- FALSE
	if ( readIDtype == "Solexa_1.8") {
		hasBarcode <- TRUE
		if ( markPtrs[ length(markPtrs)] == base::nchar( strIn[1])) {
			nTermsEach <- nTermsEach - 1
			hasBarcode <- FALSE
		}
		# some datasets will have some pieces trimmed away already
		if( max( wantedTerms) > nTermsEach) {
			wantedTerms <- 4:nTermsEach
			hasBarcode <- FALSE
		}
	}
	# some 1.3 readIDs have colons in the prefix before the lane,tile,etc. part
	if ( readIDtype == "Solexa_1.3") {
		if ( length(markPtrs) > 4) {
			nn <- length(markPtrs)
			wantedTerms <- (nn-4+2) : (nn+1)
		}
	}
	# some Rosetta readIDs have text before the lane number
	if ( readIDtype == "Rosetta") {
		if ( length(markPtrs) > 3) {
			nn <- length(markPtrs)
			wantedTerms <- (nn-2) : (nn+1)
		}
	}


	allTerms <- strsplit( strIn, split=splitMark, fixed=fixed)
	mTerms <- matrix( base::unlist(allTerms), nrow=nTermsEach, ncol=length(strIn))
	if (skipChar > 0 && wantedTerms[1] == 1) {
		mTerms[1, ] <- base::substr( mTerms[1, ], skipChar+1, 10)
	}
	outLane <- as.integer( mTerms[ wantedTerms[1], ])
	if ( any( is.na(outLane))) {
		cat("\nFound non-integer Lane Numbers.    Tried ID field: ", wantedTerms[1])
		cat( "\nEncounterd:  ",  mTerms[ wantedTerms[1], head( which( is.na(outLane)))])
	}
	outTile <- as.integer( mTerms[ wantedTerms[2], ])
	if ( any( is.na(outTile))) {
		cat("\nFound non-integer Tile Numbers.    Tried ID field: ", wantedTerms[2])
		cat( "\nEncounterd:  ",  mTerms[ wantedTerms[2], head( which( is.na(outTile)))])
	}
	outX <- as.integer( mTerms[ wantedTerms[3], ])
	if ( any( is.na(outX))) {
		cat("\nFound non-integer X_position Values.    Tried ID field: ", wantedTerms[3])
		cat( "\nEncounterd:  ",  mTerms[ wantedTerms[3], head( which( is.na(outX)))])
	}
	# Y term may have junk...
	if ( readIDtype %in% c( "Rosetta", "Solexa_1.3", "Solexa_1.8")) {
		tmp <- mTerms[ wantedTerms[4], ] 
		tmp <- sub( "( |#|/).*$", "", tmp) 
		outY <- as.integer( tmp)
	} else {
		outY <- as.integer( mTerms[ wantedTerms[4], ])
	}
	if ( any( is.na(outY))) {
		cat("\nFound non-integer Y_position Values.    Tried ID field: ", wantedTerms[4])
		cat( "\nEncounterd:  ",  mTerms[ wantedTerms[4], head( which( is.na(outY)))])
	}

	isFiltered <- barcode <- NULL
	if ( readIDtype %in% c( "Solexa_1.8")) {
		if ( length( wantedTerms) >= 6) isFiltered <- mTerms[ wantedTerms[6], ]
		if ( length( wantedTerms) >= 8 && hasBarcode) barcode <- mTerms[ wantedTerms[8], ]
	}

	if ( any( is.na(outLane))) cat("\nBad Lane Numbers: ", mTerms[ wantedTerms[1], 
			head( which( is.na(outLane)))])
	
	if ( ! all( outLane > 0)) {
		cat("\nExtracting ReadID Terms failed:   readID sample=\n")
		print( head( txt))
	}

	# allow a combo of Lane+Tile
	outCombo <- paste( "L", outLane, "T", outTile, sep="")

	return( list( "Lane"=outLane, "Tile"=outTile, "X"=outX, "Y"=outY, "Filter"=isFiltered,
			"Barcode"=barcode, "Combo"=outCombo))
}



`partitionFastqByReadIDfield` <- function( filein, field=c( "Lane", "Combo", "Filter"), chunkSize=100000) {

	# partition one fastq file of multiple lanes of data into separate .fastq files for each lane.

	field <- match.arg( field)

	fileToUse <- allowCompressedFileName( filein)
	if ( ! file.exists( fileToUse)) stop( paste("Can't find input file: ", fileToUse))
	compressOutput <- ( regexpr( ".gz$", fileToUse) > 0)

	# get the type of readIDs
	ans <- detectFastqReadFormat( fileToUse)
	readIDtype <- ans$readIDtype
	cat( "\nFastq ReadID format:  ", readIDtype, "\n")

	conIn <- openCompressedFile( fileToUse, open="r")

	lineChunkSize <- chunkSize * 4
	nread <- 0
	fieldSet <- countSet <- fileSet <- vector()
	conSet <- list()
	partitionOverhead <- list( "fields"=fieldSet, "cons"=conSet, "counts"=countSet, "files"=fileSet)

	repeat {
		chunk <- readLines( conIn, n=lineChunkSize)
		if ( length( chunk) < 1) break

		nread <- nread + length( chunk)
		cat( "\rN_Reads: ", formatC( as.integer(nread/4), big.mark=","))
		partitionOverhead <- partitionFastq( chunk, field=field, readIDtype=readIDtype, 
				overhead=partitionOverhead, basefile=basename( fileToUse), 
				compressOutput=compressOutput)
		nParts <- length( partitionOverhead$counts)
		cat ( "  ", nParts, "  ", partitionOverhead$fields, "   ", formatC( as.integer( 
				partitionOverhead$counts), big.mark=","))
	}

	close( conIn)
	cat( "\nN_reads in original:  ", formatC( as.integer(nread/4), big.mark=","))
	for( i in 1:length( partitionOverhead$fields)) {
		thisField <- partitionOverhead$fields[ i]
		con <- partitionOverhead$cons[[ i]]
		cnt <- partitionOverhead$counts[ i]
		fil <- partitionOverhead$files[ i]
		close( con)
		cat( "\n",field, "=", thisField, "\tWrote file: ", fil, "\tN_Reads: ", formatC( as.integer(cnt),
				big.mark=","))
	}
	cat( "\n")
	return()
}


`partitionFastq` <- function( txt, field="Lane", readIDtype="Rosetta", overhead, basefile,
			compressOutput=FALSE) {

	# given a chunk of .fastq file text, use a field in the id to partition...

	# the id line is 1st, and it has the "lane" number in it...
	idsTxt <- txt[ seq( 1, length(txt), by=4)]
	idTerms <- extractReadIDterms( idsTxt, readIDtype=readIDtype)

	if ( ! (field %in% names(idTerms))) {
		cat( "\nInvalid ReadID field:  ", field, 
			"\nChoices are:   ", names( idTerms), "\n")
		stop()
	}

	# factor by lane number
	myFields <- idTerms[[ field]]
	fieldFac <- factor( myFields)
	thesefields <- levels( fieldFac)

	# any new lanes get their files, etc., set up...
	for ( i in 1:length(thesefields)) {
		thisfield <- thesefields[i]
		if ( thisfield %in% overhead$fields) next
		# this is a new lane, add it, open a new connection, etc.
		nnow <- length( overhead$fields) + 1
		overhead$fields[nnow] <- thisfield

		if ( regexpr( ".fastq", basefile, fixed=TRUE) > 0) {
			fileout <- sub( ".fastq", paste( ".", field, thisfield, ".fastq", sep=""), basefile, fixed=TRUE)
		} else if ( regexpr( ".fq", basefile, fixed=TRUE) > 0) {
			fileout <- sub( ".fq", paste( ".", field, thisfield, ".fq", sep=""), basefile, fixed=TRUE)
		} else {
			fileout <- paste( basefile, ".", field, thisfield, sep="")
		}
		if (compressOutput) {
			if ( regexpr( ".gz$", fileout) < 1) {
				fileout <- paste( fileout, "gz", sep=".")
			}
			con <- gzfile(fileout, open="wb")
		} else {
			con <- file(fileout, open="wt")
		}
		overhead$files[nnow] <- fileout
		overhead$cons[[nnow]] <- con
		overhead$counts[nnow] <- 0
	}

	# now partition those lines of text
	for ( field in thesefields) {
		whichfield <- base::match( field, overhead$fields)
		mycon <- overhead$cons[[ whichfield]]
		mycnt <- overhead$counts[ whichfield]

		who <- which( myFields == field)
		# turn these read locations into the actual 4 line pointers for each read
		starts <- ((who-1) * 4) + 1
		ends <- starts + 3
		lines <- base::unlist( base::mapply( FUN=`:`, starts, ends))

		# write those lines to that file
		writeLines( txt[ lines], con=mycon)

		# update the overhead
		overhead$counts[ whichfield] <- mycnt + length( starts)
	}

	return( overhead)
}


`removeOpticalDuplicates` <- function( fastq1, fastq2=NULL, readKey.len=32, max.XY.dist=250, max.edit.dist=3, 
					fout1=NULL, fout2=NULL, verbose=TRUE) {

	# search one or two fastq file(s) of FASTQ data, and drop optical duplicate reads
	# i.e. exact matches close together on the slide

	fastq1 <- allowCompressedFileName( fastq1)
	if ( ! file.exists( fastq1)) stop( paste("Can't find FASTQ file: ", fastq1))
	compressOutput <- grepl( ".gz$", fastq1)[1]
	if ( is.null(fastq2)) {
		doPairs <- FALSE
	} else {
		fastq2 <- allowCompressedFileName( fastq2)
		if ( ! file.exists( fastq2)) stop( paste("Can't find mate 2 FASTQ file: ", fastq2))
		doPairs <- TRUE
	}

	# create filenames for the output
	if ( is.null( fout1)) {
		fout1 <- fastq1
		if ( grepl( "\\.fastq",fout1)) fout1 <- sub( "\\.fastq", "DropOptDups.fastq", fout1)
		if ( grepl( "\\.fq",fout1)) fout1 <- sub( "\\.fq", "DropOD.fq", fout1)
		if (fout1 == fastq1) stop( "Unable to create new FASTQ filename.  Looked for 'fastq' or 'fq'.")
	}
	if ( doPairs && is.null( fout2)) {
		fout2 <- fastq2
		if ( grepl( "\\.fastq",fout2)) fout2 <- sub( "\\.fastq", "DropOptDups.fastq", fout2)
		if ( grepl( "\\.fq",fout2)) fout2 <- sub( "\\.fq", "DropOD.fq", fout2)
		if (fout2 == fastq2) stop( "Unable to create new FASTQ filename.  Looked for 'fastq' or 'fq'.")
	}

	# get the type of readIDs
	ans <- detectFastqReadFormat( fastq1)
	# get the type of readIDs
	ans <- detectFastqReadFormat( fastq1)
	readIDtype <- ans$readIDtype
	if (verbose) cat( "\nFastq ReadID format:  ", readIDtype, "\n")

	# open all the files we will read/write
	conIn1 <- openCompressedFile( fastq1, open="r")
	conOut1 <- openCompressedFile( fout1, open="w")
	if (doPairs) {
		conIn2 <- openCompressedFile( fastq2, open="r")
		conOut2 <- openCompressedFile( fout2, open="w")
	}

	# because of the way Illumina writes out the raw reads in geographic order along the tile, we only need a fixed
	# amount of look-ahead, based on how far away we want to search.  Empirically derived...
	readChunkSize <- ceiling( (max.XY.dist * 5)/100) * 100
	cat( "\nBufferSize: ", readChunkSize, "  XY.Dist: ", max.XY.dist, "  Edit.Dist: ", max.edit.dist)
	cat( "\n")
	lineChunkSize <- readChunkSize * 4

	# we will process one chunk at a time, but keep the "next" buffer in memory for faster look-ahead
	firstPass <- TRUE
	hasMoreData <- TRUE
	nRead <- totalDups <- 0
	doubleBuffered <- FALSE

	# accumulate a few stats
	ttlXYdist <- maxXYdist <- ttlLookAhead <- maxLookAhead <- 0

	# read in a buffer of reads and deduce their locations
	nBufferReads <- 0
	repeat {
		if ( ! hasMoreData) break

		# read the next buffer full
		chunk1 <- readLines( conIn1, n=lineChunkSize)
		if ( ! (Nlines1 <- length( chunk1))) hasMoreData <- FALSE
		if (doPairs && hasMoreData) {
			chunk2 <- readLines( conIn2, n=lineChunkSize)
			Nlines2 <- length( chunk2)
			if ( Nlines2 != Nlines1) stop( "Mater pair files not same size!")
		}
		nRead <- nRead + (Nlines1/4)
		nBufferReads <- nBufferReads + 1

		if (hasMoreData) {
			ptrs <- seq( 1, Nlines1, by=4)
			thisID1 <- chunk1[ ptrs]
			thisRead1 <- chunk1[ ptrs+1]
			thisQual1 <- chunk1[ ptrs+3]
			if (doPairs) {
				thisID2 <- chunk2[ ptrs]
				thisRead2 <- chunk2[ ptrs+1]
				thisQual2 <- chunk2[ ptrs+3]
			}
	
			# extract the geographic info
			idTerms <- extractReadIDterms( thisID1, readIDtype=readIDtype)
			thisTile <- idTerms$Tile
			thisX <- idTerms$X
			thisY <- idTerms$Y

			if (firstPass) {
				id1b <- thisID1
				read1b <- thisRead1
				qual1b <- thisQual1
				tile1b <- thisTile
				x1b <- thisX
				y1b <- thisY
				if (doPairs) {
					id2b <- thisID1
					read2b <- thisRead2
					qual2b <- thisQual2
				}
				firstPass <- FALSE
				doubleBuffered <-TRUE
				next
			}
		}

		# move buffer 2 up to the primary storage
		if ( ! exists( "id1b")) break
		id1a <- id1b
		read1a <- read1b
		qual1a <- qual1b
		tile1a <- tile1b
		x1a <- x1b
		y1a <- y1b
		if (doPairs) {
			id2a <- id2b
			read2a <- read2b
			qual2a <- qual2b
		}
		# and write to the second
		if (hasMoreData) {
			id1b <- thisID1
			read1b <- thisRead1
			qual1b <- thisQual1
			tile1b <- thisTile
			x1b <- thisX
			y1b <- thisY
			if (doPairs) {
				id2b <- thisID1
				read2b <- thisRead2
				qual2b <- thisQual2
			}
		} else {
			# no more new data
			if ( exists( "id1b")) {
				rm( id1b, read1b, qual1b, tile1b, x1b, y1b)
				if ( doPairs) rm( id2b, read2b, qual2b)
				doubleBuffered <- FALSE
			}
		}

		# OK, one or 2 buffers have all the read ID data we need to find Optical Dups.
		if (doubleBuffered) {
			myReads1 <- c( read1a, read1b)
			if (doPairs) myReads2 <- c( read2a, read2b)
			myTiles <- c( tile1a, tile1b)
			myX <- c( x1a, x1b)
			myY <- c( y1a, y1b)
		} else {
			myReads1 <- read1a
			if (doPairs) myReads2 <- read2a
			myTiles <- tile1a
			myX <- x1a
			myY <- y1a
		}

		# set up to trap all the duplicates in the first buffer, 
		# using a factoring of the first KK characters of the read
		bufferMax <- length( read1a)
		optDupsToDrop <- optDupXYdist <- optDupLookAhead <- vector()
		nDrops <- 0
		readKey <- substr( myReads1, 1, readKey.len)
		readFac <- factor( readKey)

		tapply( 1:length(myReads1), readFac, function(x) {
				# need at least 2
				if ((NX <- length(x)) < 2) return()
				# any entries in buffer 1?
				toTest <- which( x <= bufferMax)
				for ( k in toTest) {
					# if we have no others to compare to, quit
					if (k == NX) break
					x1 <- x[k]
					r1 <- myReads1[x1]
					xOthers <- x[ (k+1):NX]
					rOthers <- myReads1[ xOthers]
					# find how many edits apart the reads are
					editDist <- adist( r1, rOthers)
					goodByEdits <- which( editDist <= max.edit.dist)
					if ( ! length(goodByEdits)) next
					# check the XY distance
					tileOrigin <- myTiles[ x1]
					XOrigin <- myX[ x1]
					YOrigin <- myY[ x1]
					tileOthers <- myTiles[ xOthers]
					XOthers <- myX[ xOthers]
					YOthers <- myY[ xOthers]
					# weight the distance so a tile change is extreme
					dOthers <- sqrt( (XOthers-XOrigin)^2 + (XOthers-XOrigin)^2 + 
								((tileOthers-tileOrigin)*1000)^2)
					goodByXY <- which( dOthers <= max.XY.dist)
					if ( ! length(goodByXY)) next
					# to be a valid OD, must be both low edits and distance
					good <- intersect( goodByEdits, goodByXY)
					if ( ! length(good)) next
					# almost done...
					# if we have mate pairs, then the second read must also be just as close in edit space
					if (doPairs) {
						r2 <- myReads2[x1]
						r2Others <- myReads2[ xOthers[ good]]
						editDist2 <- adist( r2, r2Others)
						if ( all( editDist2 > max.edit.dist)) next
					}
					# this tested read is the same as at least 1 other read further ahead in the file
					nDrops <<- nDrops + 1
					optDupsToDrop[ nDrops] <<- x1
					optDupXYdist[ nDrops] <<- min( dOthers[good])
					optDupLookAhead[ nDrops] <<- min( xOthers[good]) - x1
				}
				return()
			})

		# at this point we have a set of reads to omit
		outID1 <- id1a
		outREAD1 <- read1a
		outQUAL1 <- qual1a
		if ( nDrops) {
			outID1 <- outID1[ -optDupsToDrop]
			outREAD1 <- outREAD1[ -optDupsToDrop]
			outQUAL1 <- outQUAL1[ -optDupsToDrop]
		}
		outTxt <- paste( outID1, outREAD1, "+", outQUAL1, sep="\n")
		writeLines( outTxt, con=conOut1, sep="\n")
		if (doPairs) {
			outID2 <- id2a
			outREAD2 <- read2a
			outQUAL2 <- qual2a
			if ( nDrops) {
				outID2 <- outID2[ -optDupsToDrop]
				outREAD2 <- outREAD2[ -optDupsToDrop]
				outQUAL2 <- outQUAL2[ -optDupsToDrop]
			}
			outTxt <- paste( outID2, outREAD2, "+", outQUAL2, sep="\n")
			writeLines( outTxt, con=conOut2, sep="\n")
		}

		totalDups <- totalDups + nDrops
		ttlXYdist <- ttlXYdist + sum( optDupXYdist,na.rm=T)
		maxXYdist <- round( max( maxXYdist, optDupXYdist))
		ttlLookAhead <- ttlLookAhead + sum( optDupLookAhead,na.rm=T)
		maxLookAhead <- max( maxLookAhead, optDupLookAhead)
		if ( nBufferReads %% 10 == 0) {
			pctOptDup <- as.percent( totalDups, big.value=nRead)
			avgXYdist <- round( ttlXYdist / totalDups, digits=2)
			avgLookAhead <- round( ttlLookAhead / totalDups, digits=0)
			cat( "\rN_Read:", nRead, " N_OptDup:", totalDups, " Pct:", pctOptDup, "  AvgXYdist:", avgXYdist, 
					" MaxXYdist:", maxXYdist, " AvgAhead:", avgLookAhead, 
					" MaxAhead:", maxLookAhead)
		}
		# continue reading and thinning until no more data
	}

	close( conIn1)
	close( conOut1)
	cat( "\nWrote file: ", fout1)
	if (doPairs) {
		close( conIn2)
		close( conOut2)
		cat( "\nWrote file: ", fout2)
	}
	cat( "\n")
	return( totalDups)
}


`opticalDuplicates` <- function( filein, chunkSize=100000, key.length=32, max.XY.dist=250, max.edit.dist=3, verbose=TRUE) {

	# search one fastq file of single lane/tile of raw read data for optical duplicates
	# i.e. exact matches close together in space

	fileToUse <- allowCompressedFileName( filein)
	if ( ! file.exists( fileToUse)) stop( paste("Can't find input file: ", fileToUse))
	compressOutput <- ( regexpr( ".gz$", fileToUse) > 0)

	# get the type of readIDs
	ans <- detectFastqReadFormat( fileToUse)
	readIDtype <- ans$readIDtype
	if (verbose) cat( "\nFastq ReadID format:  ", readIDtype, "\n")

	conIn <- openCompressedFile( fileToUse, open="r")

	lineChunkSize <- chunkSize * 4
	nread <- 0
	idSet <- readSet <- qualSet <- xSet <- ySet <- vector()
	nNow <- 0

	# read in all the reads and their locations
	repeat {
		chunk <- readLines( conIn, n=lineChunkSize)
		if ( (Nlines <- length( chunk)) < 1) break
		nread <- nread + Nlines
		if (verbose) cat( "\rN_Reads: ", formatC( as.integer(nread/4), big.mark=","))

		ptrs <- seq( 1, Nlines, by=4)
		thisID <- chunk[ ptrs]
		thisRead <- chunk[ ptrs+1]
		thisQual <- chunk[ ptrs+3]
		
		Nreads <- round( Nlines / 4)
		pOut <- (nNow+1) : (nNow+Nreads)
		nNow <- nNow + Nreads
		idSet[ pOut] <- thisID
		readSet[ pOut] <- thisRead
		qualSet[ pOut] <- thisQual
	
		idTerms <- extractReadIDterms( thisID, readIDtype=readIDtype)
		xSet[ pOut] <- idTerms$X
		ySet[ pOut] <- idTerms$Y
	}
	rm( thisID, thisRead, thisQual, idTerms)
	close( conIn)
	cat( "\nN_reads in original:  ", formatC( as.integer(nNow), big.mark=","))

	# now find those that are matches, within some tolerance
	# use a fixed width prefix for finding reads worth testing
	readKey <- substr( readSet, 1, key.length)
	readFac <- factor( readKey)
	cat( "\nN_ReadKeys: ", nlevels(readFac), "\nSearching..\n")
	xSet <- as.numeric( xSet)
	ySet <- as.numeric( ySet)

	# visit all the reads that may be duplicates
	nODs <- 0
	od1 <- od2 <- nEdits <- oDist <- vector()
	tapply( 1:nNow, readFac, function(x) {
			if ( (NN <- length(x)) < 2) return()
			if (verbose) cat( "\r", x[1], readKey[x[1]], NN)
			# find those pairs where the read strings are similar enough
			for ( i in 1:(NN-1)) {
				r1 <- x[i]
				x1 <- xSet[ r1]
				y1 <- ySet[ r1]
				dm <- adist( readSet[r1], readSet[x[(i+1):NN]])
				for ( j in (i+1):NN) {
					thisED <- dm[ 1, j-i]
					if ( thisED > max.edit.dist) next
					r2 <- x[j]
					x2 <- xSet[ r2]
					y2 <- ySet[ r2]
					thisOD <- sqrt( (x1-x2)^2 + (y1-y2)^2)
					if ( thisOD > max.XY.dist) next
					# we found a duplicate close in space!
					nODs <<- nODs + 1
					od1[nODs] <<- r1
					od2[nODs] <<- r2
					nEdits[nODs] <<- thisED
					oDist[nODs] <<- thisOD
					if (verbose) cat( "\rHit: ", x[1], nODs, i, j, thisED, thisOD)
				}
			}
		})
	
	# summarize what we found
	cat( "\nN_OpticalDuplicates: ", nODs, "\n")
	out <- data.frame( "ID1"=idSet[od1], "ID2"=idSet[od2], "Read1"=od1, "Read2"=od2, 
			"EditDist"=nEdits, "OpticalDistance"=oDist, stringsAsFactors=F)
	return( out)
}


`seq2fastq` <- function( filein, fileout, phred=TRUE, split="\t") {

	# turn a .seq file from Rosetta ( compressed or not) having Solexa-type quality scores
	# into a .fastq file with Phred-type quality scores ( compressed or not)

	filein <- allowCompressedFileName(filein)
	if ( ! file.exists( filein)) stop( paste("Can't find input file: ", filein))
	conIn <- openCompressedFile( filein, open="r")

	if ( regexpr( ".gz$", fileout) > 0) {
		conOut <- gzfile( fileout, open="w")
	} else {
		conOut <- file( fileout, open="w")
	}

	cat( "\nInput File: ", filein)

	chunkSize <- 100000
	nread <- 0

	repeat {
		chunk <- readLines( conIn, n=chunkSize)
		if ( length( chunk) < 1) break

		nread <- nread + length( chunk)
		newchunk <- illumina2Fastq( chunk, phred=phred, split=split)
		writeLines( newchunk, conOut)
		cat( "  ",formatC(nread, format="d", big.mark=","))
	}

	close( conIn)
	close( conOut)
	cat( "\nN_reads converted: ", nread, "\n")
}


`randomDNA2fastq` <- function( read.length=100, fileout, Nreads=1000000, mean.qual=36, sd.qual=2) {

	# generate random k-mers
	# into a .fastq file with Phred-type quality scores ( compressed or not)
	if ( regexpr( ".gz$", fileout) > 0) {
		conOut <- gzfile( fileout, open="w")
	} else {
		conOut <- file( fileout, open="w")
	}

	cat( "\nMaking ", Nreads, " random reads of length ", read.length, "  Mean.Q: ", mean.qual, "  SD.Q: ", sd.qual, "\n")

	chunkSize <- 100000
	nmade <- 0
	BASES <- c( "A", "C", "G", "T")
	PASTE <- base::paste

	repeat {
		if ( nmade >= Nreads) break
		nbase <- chunkSize * read.length
		baseVec <- sample( BASES, size=nbase, replace=T)
		qualVec <- round( rnorm( nbase, mean=mean.qual, sd=sd.qual)) + 33
		starts <- seq( 1, nbase, by=read.length)
		stops <- starts + read.length - 1
		N <- length(starts)
		readSeq <- readQual <- vector( length=N)
		sapply( 1:N, FUN=function(i) {
				readSeq[i] <<- PASTE( baseVec[ starts[i]:stops[i]], collapse="")
				readQual[i] <<- intToUtf8( qualVec[ starts[i]:stops[i]])
				return(NULL)
			})
		nmade <- nmade + N
		newchunk <- paste( "@", "Random", 1:N, "\n", readSeq, "\n+\n", readQual, sep="")
		writeLines( newchunk, conOut)
		cat( "\rN_Wrote:  ",formatC(nmade, format="d", big.mark=","))
	}
	close( conOut)
	cat( "\nN_reads created: ", nmade, "\n")
}


`randomDNA2fasta` <- function( read.length=100, fileout, Nreads=1000000) {

	# generate random k-mers
	# into a .fasta file
	conOut <- file( fileout, open="w")

	cat( "\nMaking ", Nreads, " random reads of length ", read.length, "\n")

	chunkSize <- 100000
	nmade <- 0
	BASES <- c( "A", "C", "G", "T")
	PASTE <- base::paste

	repeat {
		if ( nmade >= Nreads) break
		nbase <- chunkSize * read.length
		baseVec <- sample( BASES, size=nbase, replace=T)
		starts <- seq( 1, nbase, by=read.length)
		stops <- starts + read.length - 1
		N <- length(starts)
		readSeq <- vector( length=N)
		sapply( 1:N, FUN=function(i) {
				readSeq[i] <<- PASTE( baseVec[ starts[i]:stops[i]], collapse="")
				return(NULL)
			})
		nmade <- nmade + N
		newchunk <- paste( ">", "Random", 1:N, "\n", readSeq, sep="")
		writeLines( newchunk, conOut)
		cat( "\rN_Wrote:  ",formatC(nmade, format="d", big.mark=","))
	}
	close( conOut)
	cat( "\nN_reads created: ", nmade, "\n")
}


`illumina2Fastq` <- function( txt, phred=TRUE, split=":") {

	# given lines of text from a illumina .SEQ file of raw reads,
	# transform it into a .FASTQ file for input to an alignment program
	terms <- strsplit( txt, split=split, fixed=TRUE)

	# we should trap any bad records before we 'matrix' them...
	nExpect <- 7
	if ( split == "\t") nExpect <- 11
	nGot <- sapply( terms, length)
	if ( (nUse <- median( nGot)) > nExpect) nExpect <- nUse

	badOnes <- which( nGot < nExpect)
	if ( length( badOnes) > 0) {
		terms <- terms[ -badOnes]
		badtxt <- txt[ badOnes]
		cat( "\nCorrupt Lines:   skipping ", length(badOnes),"\n")
		print( badtxt)
	}

	flat <- base::unlist( terms)
	termsM <- matrix( flat, nrow=nExpect, ncol=length( terms))

	# older Illumina style was 'colon' separated...
	if ( split == ":") {
	
	# first field is id, ignore for now
	# next four are "lane / tile / X / Y"
	myID <- base::paste( termsM[2,], termsM[3,], termsM[4,], termsM[5,], sep="_")
	laneID <- termsM[ 2, ]

	# next is the RNA read string
	mySEQ <- termsM[ 6, ]

	# last is the quality scores, as space-separated "Solexa-type" integer scores
	# OR as a "Solexa-type" character string...
	myScores <- termsM[ 7, ]

	# try to detect which...
	hasSpaces <- regexpr( " ", myScores, fixed=TRUE)
	if ( all( hasSpaces > 0)) {
		# nothing to do....
	} else {
		fac <- factor( myScores)
		uniqueScores <- levels( fac)
		facPtrs <- tapply( 1:length(myScores), fac, FUN=NULL)
		cat( "  Solexa Conv Speedup=", formatC( (length(myScores)/length(uniqueScores)), digits=2, format="f"))
		uniqueSolexaStrs <- sapply( uniqueScores, FUN=solexaCharScoreToSolexaIntScore)
		myScores <- uniqueSolexaStrs[ facPtrs]
	}
	# convert to Phred style scores...
	if ( phred) {
		fac <- factor( myScores)
		uniqueScores <- levels( fac)
		facPtrs <- tapply( 1:length(myScores), fac, FUN=NULL)
		cat( "  Phred Conv Speedup=", formatC( (length(myScores)/length(uniqueScores)), digits=2, format="f"))
		uniquePhredStrs <- sapply( uniqueScores, FUN=solexaToPhred)
		myScores <- uniquePhredStrs[ facPtrs]
	}
	}

	if ( split == "\t") {  #if split is by tab, newer format with different spacing, etc...

	# next four are "lane / tile / X / Y"
	myIDa <- base::paste( termsM[1], termsM[2,], sep="_") 
	myIDb <- base::paste( termsM[3,], termsM[4,], termsM[5,], termsM[6,], sep=":")
	myIDc <- base::paste( termsM[7], termsM[8,], sep="/") 
	myID <- base::paste( myIDa, ":", myIDb, "#", myIDc, sep="")
	laneID <- termsM[ 3, ]

	# next is the RNA read string
	mySEQ <- termsM[ 9, ]

	# "Solexa-type" character string...
	myScores <- termsM[ 10, ]
	}


	# make the fastq format lines...
	out <- base::paste( "@", myID, "\n", mySEQ, "\n+\n", myScores, sep="")
	return( out)
}


`fastqPatternSearchDetails` <- function( filein, pattern="AGAGATCGGAAGATCTCGTATGCC") {

	# find all reads that contain a pattern string

	filein <- allowCompressedFileName( filein)
	if ( ! file.exists( filein)) stop( paste("Can't find input file: ", filein))

	ans <- detectFastqReadFormat( filein)
	readIDtype <- ans$readIDtype

	conIn <- openCompressedFile( filein, open="r")

	chunkSize <- 400000
	nreads <- nhits <- 0

	bigTileYes <- bigXYes <- bigYYes <- vector()
	bigTileNo <- bigXNo <- bigYNo <- vector()

	repeat {
		txt <- readLines( conIn, n=chunkSize)
		if ( length( txt) < 1) break

		nread <- round(length( txt)/4)
		nreads <- nreads + nread
		cat( "\nN_Reads: ", prettyNum( as.integer(nreads), big.mark=","))

		# extract all the lane,tile,X,Y info from the read IDs
		idsTxt <- txt[ seq( 1, length(txt), by=4)]
		idTerms <- extractReadIDterms( idsTxt, readIDtype=readIDtype)

		# turn the reads into a searchable set, and count how often we see the pattern
		readTxt <- txt[ seq( 2, length(txt), by=4)]
		subject <- DNAStringSet( readTxt)

		cnts <- vcountPattern( pattern, subject, max.mismatch=1)
		hits <- which( cnts > 0)
		nhit <- length( hits)
		nhits <- nhits + nhit
		cat( "  N_hits: ", nhit, "  Percent Hits: ", as.percent( nhit, big.value=nread))

		bigTileYes <- append( bigTileYes, idTerms$Tile[ hits])
		bigXYes <- append( bigXYes, idTerms$X[ hits])
		bigYYes <- append( bigYYes, idTerms$Y[ hits])

		noHits <- which( cnts == 0)
		bigTileNo <- append( bigTileNo, idTerms$Tile[ noHits])
		bigXNo <- append( bigXNo, idTerms$X[ noHits])
		bigYNo <- append( bigYNo, idTerms$Y[ noHits])
	}

	close( conIn)
	cat( "\nN_reads in file:  ", nreads)

	return( list( "TileYes"=bigTileYes, "XYes"=bigXYes, "YYes"=bigYYes, "TileNo"=bigTileNo, 
			"XNo"=bigXNo, "YNo"=bigYNo))
}


`SOLiD2fastq` <- function( csfasta.file, qual.file=sub( ".csfasta", "_QV.qual", csfasta.file), 
			 	fq.file=sub( ".csfasta", ".fq.gz", csfasta.file), comment.lines=3) {

	conIn <- file( csfasta.file, open="rt")
	cat( "\nReading colorspace fasta file: ", csfasta.file)
	conInQ <- file( qual.file, open="rt")
	cat( "\nReading colorspace quality file: ", qual.file)
	conOut <- gzfile( fq.file, open="w")

	# get past the comments...
	for(i in 1:comment.lines) {
		id1 <- readLines( conIn, n=1)
		id2 <- readLines( conInQ, n=1)
		if ( substr( id1, 1,1) != "#") break
	}

	chunkSize = 200000
	cat( "\nConverting..")
	nout <- 0
	repeat {
		txt <- readLines( conIn, n=chunkSize)
		txtQ <- readLines( conInQ, n=chunkSize)
		if ( length( txt) < 1) break

		who <- seq( 1,length(txt), by=2)
		id1 <- txt[ who]
		id2 <- txtQ[ who]
		seq <- txt[ who + 1]
		qualstr <- txtQ[ who + 1]

		# make the id look like an early 'Rosetta' style of Illumina ID
		id <- base::sub( ">", "@1_", id1)
		#seq <- convColorSpaceToACGT( seq)
		seq <- fastCS2DNA( seq)
		qualstr <- base::sapply( qualstr, solexaToPhred)

		writeLines( base::paste( id, seq, "+", qualstr, sep="\n"), con=conOut)
		if ( id1[1] != id2[1]) cat( "\rBad IDs: ", id1[1], id2[1])
		nout <- nout + length(who)
		cat( ".")
		if ( nout %% 1000000 == 0) cat( formatC(nout, format="d", big.mark=","))
	}
	close( conIn)
	close( conInQ)
	close( conOut)
	cat( "\nWrote file: ", fq.file, "\nN_Reads: ", formatC( nout, format="d", big.mark=","), "\n")
	return(NULL)
}


`convColorSpaceToACGT` <- function( txt) {

	A_JUMP <- c( "A","C","G","T","N")
	C_JUMP <- c( "C","A","T","G","N")
	G_JUMP <- c( "G","T","A","C","N")
	T_JUMP <- c( "T","G","C","A","N")

	v <- strsplit( txt, split="")
	ans <- base::sapply( v, function(x) {
		N <- length(x)
		cur <- x[1]
		x[ x == "."] <- '4'
		xin <- as.numeric( x[2:N]) + 1
		out <- vector()
		now <- ""
		for ( j in 1:(N-1)) {
			now <- if ( cur == 'A') A_JUMP[ xin[j]]
				else if ( cur == 'C') C_JUMP[ xin[j]]
				else if ( cur == 'G') G_JUMP[ xin[j]]
				else if ( cur == 'T') T_JUMP[ xin[j]]
				else 'N'
			out[j] <- now
			if ( now != "N") cur <- now
		}
		return( base::paste( out, collapse=""))
		})

	return(ans)
}


