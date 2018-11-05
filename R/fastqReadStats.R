# fastqReadStats.R


# spawn the fastq stats guy...
`dispatch.FastqReadStats` <- function( filein, sampleID, statsPath=".", chunkSize=100000, maxReads=NULL, 
				pause=0, baseOrder="ACGTN") {

	commandLine <- paste( "checkX11( width=7, height=8, xpos=10, ypos=10, bg='white'); ",
				" fastqReadStats( filein=\"", filein, 
				"\", sampleID=\"", sampleID, 
				"\", statsPath=\"", statsPath, "\", chunkSize=", as.integer(chunkSize), 
				", maxReads=", if (is.null(maxReads)) "NULL" else as.integer(maxReads), 
				", pause=", as.integer(pause),
				", baseOrder=\"", baseOrder, "\" )", sep="")

	spawnRcmd( commandLine, Rpackage="DuffyNGS", logFile=paste( sampleID, "fastqReadStats.log.txt", sep="."))
	return()
}


# assess read qualities

`fastqReadStats` <- function( filein, sampleID, statsPath="FastqReadStats", calcStats=TRUE, plotAllTiles=FALSE, 
			baseOrder="ACGTN", dropPatterns=NULL, chunkSize=100000, maxReads=NULL, 
			keepAllTiles=NULL, max.lanes.storage=3, max.tiles.storage=240, pause=0) {

	sampleStatsPath <- file.path( statsPath, sampleID)
	statsFile <- paste( sampleID, "rawReadStats.rda", sep=".")
	statsFile <- file.path( sampleStatsPath, statsFile)
	if ( ! file.exists( sampleStatsPath)) dir.create( sampleStatsPath, recursive=TRUE, showWarnings=TRUE)

	plotPath <- sampleStatsPath

	# only calculate if we need to
	if (calcStats || !file.exists( statsFile)) {
		filein <- allowCompressedFileName( filein)
		ans <- detectFastqReadFormat( filein)
		cat( "\n\nFastq file encoding:  ", filein)
		cat( "\n  Read ID type:       ", ans$readIDtype)
		cat( "\n  Score type:         ", ans$scoreType, "\n")

		calcQualityScores( filein, sampleID, statsFile=statsFile, plotPath=plotPath, 
				readIDtype=ans$readIDtype, scoreType=ans$scoreType, 
				baseOrder=baseOrder, dropPatterns=dropPatterns, chunkSize=chunkSize,
				maxReads=maxReads, keepAllTiles=keepAllTiles, 
				max.lanes.storage=max.lanes.storage, max.tiles.storage=max.tiles.storage,
				pause=pause)
	}

	ansScores <- plotQualityScores( asPNG=T, statsFile=statsFile, plotPath=plotPath, baseOrder=baseOrder) 
	cat( "\n\nFastq Score Quality Summary:\n");  print( ansScores)

	ansTiles <- plotQualityTileAverage( asPNG=T, statsFile=statsFile, plotPath=plotPath, baseOrder=baseOrder) 
	cat( "\n\nFastq Tile Quality Summary:\n");  print( ansTiles)

	plotQualityTileVariance( asPNG=T, statsFile=statsFile, plotPath=plotPath, baseOrder=baseOrder) 
	if ( plotAllTiles) {
		load( statsFile)
		tileNumbers <- as.integer( qualPlotTileIDs)
		plotQualityScores( asPNG=T, statsFile=statsFile, plotPath=plotPath, tileNumber=tileNumbers, 
				baseOrder=baseOrder) 
	}

	return( list( "ScoreSummary"=ansScores, "TileSummary"=ansTiles))
}



`getBaseOrder` <- function( baseOrder="ACGTN") {

	# we have a default ordering for base calls, but allow the user to permute it...
	baseOrder <- strsplit( toupper( baseOrder), split="", fixed=TRUE)[[1]]
	baseOrder <- unique.default( c( baseOrder, c("A","C","G","T","N")))
	where <- base::match( baseOrder, c("A","C","G","T","N"))
	return( where)
}


`calcQualityScores` <- function( filein, sampleID, statsFile="qualScoreData.Rd", plotPath=".", 
			readIDtype="Solexa_1.3", scoreType="Phred64", baseOrder="ACGTN", 
			dropPatterns=NULL, chunkSize=50000, maxReads=NULL, keepAllTiles=NULL, 
			max.lanes.storage=3, max.tiles.storage=240, pause=0) {

	# measure the base call quality scores ( and plot em') of a .fastq file

	# get ready to read the fastq file in chuncks...
	filein <- allowCompressedFileName( filein)
	if ( ! file.exists( filein)) stop( paste("Can't find input file: ", filein))
	conIn <- openCompressedFile( filein, open="r")
	on.exit( close( conIn))

	# set up the details for scores and counts
	scoreSums <- nScores <- scoreMeans <- NULL
	scoreReads <- NULL
	nBases <- 0
	MAX_LANES <- max.lanes.storage
	MAX_TILES <- max.tiles.storage
	LaneIDs <- LaneNAMES <- NULL
	TileIDs <- TileNAMES <- NULL
	scoreFile <- filein
	scoreSample <- sampleID
	countSums <- NULL
	countReps <- NULL
	nCounts <- nCBases <- nBadReads <- 0
	nBadReads <- matrix( 0, nrow=MAX_TILES, ncol=MAX_LANES)

	# set up storage for plotting details
	qualPlotFile <- qualPlotSampleID <- NULL
	qualPlotAvgScores <- NULL
	qualPlotAvgReadScores <- NULL
	qualPlotBases <- NULL
	qualPlotNscoresVec <- NULL
	qualPlotBaseCnts <- NULL
	qualPlotTileBaseCnts <- NULL
	qualPlotBaseReps <- NULL
	qualPlotTileBaseReps <- NULL
	qualPlotTileScores <- NULL
	qualPlotTileReadScores <- NULL
	qualPlotNcountsVec <- NULL
	qualPlotGrandAvgScore <- NULL
	qualPlotTotalReads <- NULL
	qualPlotLaneIDs <- qualPlotTileIDs <- NULL
	qualPlotReadType <- "Total"

	# lets capture the entire Qscore space
	allQscores <- NULL

	# older files have exactly 120 tiles, but newer formats have very large tile#s
	if ( is.null( keepAllTiles)) {
		keepAllTiles <- ( readIDtype %in% c( "Rosetta", "Solexa_1.0", "Solexa_1.3"))
	} else if ( ! is.logical( keepAllTiles)) {
		keepAllTiles <- FALSE
	}

	# local functions...


`addScores` <- function( scoresTxt, laneNumber, tileNumber, scoreType="Phred33") {

	# turn those strings into integers
	MAX_PHRED_SCORE <- getMaxPhredScore()
	m <- phredScoreStringToInt( scoresTxt, scoreType)

	# accumulate what we want to tally, binning by the lane & tile number...
	laneFac <- factor( laneNumber)	
	Nlane <- nlevels(laneFac)
	laneNames <- levels( laneFac)
	tileFac <- factor( tileNumber)
	Ntile <- nlevels(tileFac)
	tileNames <- levels( tileFac)
	laneTileCounts <- base::table( laneNumber, tileNumber)

	# set up the storage first time thru
	if ( nBases == 0) {
		nBases <<- ncol( m)
		scoreSums <<- array( 0, dim=c( nBases, MAX_TILES, MAX_LANES))
		scoreReads <<- array( 0, dim=c( MAX_PHRED_SCORE, MAX_TILES, MAX_LANES))
		nScores <<- matrix( 0, nrow=MAX_TILES, ncol=MAX_LANES)
		LaneIDs <<- laneNames
		TileIDs <<- tileNames
		newLanes <- newTiles <- NA
		if ( keepAllTiles) {
			TileIDs <<- 1:MAX_TILES
			TileNAMES <<- TileIDs
		}
	} else {
		# see if we have new tiles/lanes, and make sure they append...
		newLanes <- setdiff( laneNames, LaneIDs)
		if ( length( newLanes) > 0) LaneIDs <<- c( LaneIDs, newLanes)
		newTiles <- setdiff( tileNames, TileIDs)
		if ( length( newTiles) > 0) TileIDs <<- c( TileIDs, newTiles)
	}
	LaneNAMES <<- LaneIDs
	TileNAMES <<- TileIDs

	# as the number of lanes, tiles, grows, update the names
	blnk <- paste( "unused", 1:MAX_TILES, sep="_")
	if ( length(LaneIDs) < MAX_LANES) LaneNAMES <<- c( LaneIDs, blnk[ (length(LaneIDs)+1):MAX_LANES])
	if ( length(TileIDs) < MAX_TILES) TileNAMES <<- c( TileIDs, blnk[ (length(TileIDs)+1):MAX_TILES])
	# bail out early if we're about to exceed the memory allocation
	if ( length(LaneIDs) > MAX_LANES) {
		cat( "\nToo many Lanes detected.  Perhaps increase 'max.lanes.storage'...")
		return(NULL)
	}
	if ( length(TileIDs) > MAX_TILES) {
		cat( "\nToo many Tiles detected.  Perhaps increase 'max.tiles.storage'...")
		return(NULL)
	}
	bnames <- 1:nBases
	dimnames(scoreSums) <<- list( bnames, TileNAMES, LaneNAMES)
	snames <- 1:MAX_PHRED_SCORE
	dimnames(scoreReads) <<- list( snames, TileNAMES, LaneNAMES)
	dimnames(nScores) <<- list( TileNAMES, LaneNAMES)


	# 1: sum up the quality scores at each base position
	#    ( we will loop over all tiles per lane, so tiles first in the the list of factors!
	NCOLM <- ncol(m)
	baseSums <- tapply( 1:nrow(m), INDEX=list( tileFac, laneFac), FUN=function(x) {
			if ( length(x) < 1) return( NULL)
			out <- vector( length=NCOLM)
			for( k in 1:NCOLM) out[k] <- sum( m[ x, k])
			out
		})

	# now distribute them
	for( i in 1:Nlane) {
	    laneText <- as.character( laneNames[i])
	    lane <- match( laneNames[i], LaneIDs)
	    for( j in 1:Ntile) {
		tileText <- as.character( tileNames[j])
		tile <- match( tileNames[j], TileIDs)
		thisElement <- ((i-1)*Ntile + j)
		if ( is.null( baseSums[[ thisElement]] )) next
		scoreSums[ , tile, lane]  <<- scoreSums[ , tile, lane]  + baseSums[[ thisElement]]
		nScores[ tile, lane] <<- nScores[ tile, lane] + laneTileCounts[ laneText, tileText]
	}}

	# 2: sum up the mean quality score of each entire read
	#    ( we will loop over all tiles per lane, so tiles first in the the list of factors!
	readMeans <- tapply( 1:nrow(m), INDEX=list( tileFac, laneFac), FUN=function(x) {
			if ( length(x) < 1) return( NULL)
			out <- vector( length=length(x))
			for( k in 1:length(x)) out[k] <- mean.default( m[ x[k], ])
			out
		})

	# now distribute them
	zeros <- rep( 0, times=MAX_PHRED_SCORE)
	for( i in 1:Nlane) {
	    #laneText <- as.character( laneNames[i])
	    lane <- match( laneNames[i], LaneIDs)
	    for( j in 1:Ntile) {
		#tileText <- as.character( tileNames[j])
		tile <- match( tileNames[j], TileIDs)
		thisElement <- ((i-1)*Ntile + j)
		if ( is.null( readMeans[[ thisElement]] )) next
		# turn all these 'one mean score per read' into histogram-able counts
		# A: integerize, and clip
		tmp <- round( readMeans[[ thisElement]])
		tmp[ is.na(tmp)] <- 1
		tmp[ tmp < 1] <- 1
		tmp[ tmp > MAX_PHRED_SCORE] <- MAX_PHRED_SCORE
		# B: table, and use the names as the locations to update...
		tableOfScores <- base::table( tmp)
		tmpCnts <- zeros
		tmpCnts[ as.integer( names( tableOfScores))] <- tableOfScores
		names( tmpCnts) <- NULL
		scoreReads[ , tile, lane]  <<- scoreReads[ , tile, lane]  + tmpCnts
	}}

	return(m)
}


`addCounts` <- function( seqTxt, laneNumber, tileNumber) {

	# accumulate what we want to tally, binning by the lane number...
	laneFac <- factor( laneNumber)
	tileFac <- factor( tileNumber)

	# turn those strings into vector of single characters
	acgtList <- strsplit( seqTxt, split="", fixed=TRUE)
	alllens <- base::nchar( seqTxt)
	maxlen <- max( alllens)

	# accumulate the freq of each base
	acgt <- c( "A", "C", "G", "T", "N")
	if ( nCBases == 0) {
		nCBases <<- maxlen
		nCounts <<- matrix( 0, nrow=MAX_TILES, ncol=MAX_LANES)
		countSums <<- array( 0, dim=c( MAX_LANES, MAX_TILES, 5, nCBases), dimnames=list( NULL, NULL, acgt, 1:nCBases))
		countReps <<- array( 0, dim=c( MAX_LANES, MAX_TILES, nCBases), dimnames=list( NULL, NULL, 1:nCBases))
	}

	# as the number of lanes, tiles, grows, update the names
	bnames <- 1:nCBases
	dimnames(nCounts) <<- list( TileNAMES, LaneNAMES)
	dimnames(countSums) <<- list( LaneNAMES, TileNAMES, acgt, bnames)
	dimnames(countReps) <<- list( LaneNAMES, TileNAMES, bnames)

	# make all the base vectors equal length
	needPad <- which( alllens < maxlen)
	if ( length( needPad) > 0) {
		for( k in needPad) acgtList[[ k]] <- base::append( acgtList[[k]], rep( "", times=(maxlen-alllens[k])))
	}

	mBases <- matrix( base::unlist( acgtList), nrow=maxlen, ncol=length( seqTxt))

	newTallyBases <- function( m, lane, tile) {

		TABLE <- base::table
		BASES_CONST <- c("A","C","G","T","N")

		cnts <- apply( m, MARGIN=1, function(x) {
				tbl <- TABLE( x)
				out <- tbl[ BASES_CONST]
				out[ is.na(out)] <- 0
				if ( "." %in% names(tbl)) out[5] <- out[5] + tbl[[ "."]]
				out
			})

		countSums[ lane, tile, , ] <<- countSums[ lane, tile, , ] + cnts
		nBadReads[tile, lane] <<- nBadReads[tile, lane] + sum( cnts[ 5, ])
		reps <- apply( m, MARGIN=2, function(x) return( ( x == c( x[2:length(x)]," "))))
		# catch singularitiy...
		if ( length( dim(reps)) > 1) {
			repCnts <- apply( reps, MARGIN=1, sum)
		} else {
			repCnts <- ifelse( reps, 1, 0)
		}
		if ( length(repCnts) > 0) countReps[ lane, tile, ] <<- countReps[ lane, tile, ] + repCnts
		return()
	}


	# visit each set of lane/tile reads
	tapply( 1:ncol(mBases), INDEX=list( laneFac, tileFac), FUN=function(x) {
		if ( length(x) < 1) return()
		mylane <- laneNumber[x[1]]
		lane <- match( mylane, LaneIDs)
		mytile <- tileNumber[x[1]]
		tile <- match( mytile, TileIDs)
		tmp <- mBases[ , x, drop=FALSE]
		newTallyBases( tmp, lane=lane, tile=tile)
		nCounts[ tile, lane] <<- nCounts[ tile, lane] + ncol( tmp)
	})

	return()
}


`reportScores` <- function( verbose=TRUE) {

	nLaneScores <- apply( nScores, 2, sum)
	lanes <- which( nLaneScores > 0)
	if (verbose) {
		cat( "\nN_Lanes:            ", length( lanes))
		cat( "\nLane Names:         ", LaneIDs[ lanes])
		cat( "\nN_Scores per Lane:  ", as.integer( nLaneScores[lanes]))
		cat( "\nN_bases per read:   ", nBases, "\n")
	}

	avgScores <- matrix( nrow=nBases, ncol=length(lanes))
	rownames( avgScores) <- 1:nBases
	colnames( avgScores) <- LaneIDs[ lanes]
	for( i in 1:length(lanes)) {
		avgScores[ ,i] <- apply( scoreSums[ , ,lanes[i]], 1, sum) / nLaneScores[lanes[i]]
		if (verbose) cat( "\nAvg Scores per Base:   Lane =", LaneIDs[lanes[i]], "\n", round( avgScores[ ,i]))
	}

	nTileScores <- apply( nScores, 1, sum)
	tiles <- which( nTileScores > 0)
	if ( verbose) {
		cat( "\nN_Tiles:            ", length( tiles))
		cat( "\nTile Names:         ", TileIDs[tiles])
		cat( "\nN_Scores per Tile:  ", as.integer( nTileScores[tiles]))
	}

	avgTileScores <- array( 0, dim=c( nBases, length( tiles), length(lanes)), dimnames=list( 1:nBases,
				TileIDs[tiles], LaneIDs[lanes]))
	for( i in 1:length(lanes)) {
	    lane <- lanes[i]
	    thisNscores <- nScores[ , lane]
	    for( j in 1:length(tiles)) {
	    	tile <- tiles[j]
		avgTileScores[ ,j,i] <- scoreSums[ , tile ,lane] / thisNscores[ tile]
	}}

	GrandAvgScores <- apply( avgScores, MARGIN=2, FUN=mean.default)
	if (verbose) cat( "\n\nAvg All Base Scores per Lane: \n ", paste( "\tLane: ", LaneIDs[lanes], "\t", 
			formatC( GrandAvgScores, digits=2, format="f"), collapse="\n"))

	MAX_PHRED_SCORE <- getMaxPhredScore()
	avgReadScores <- matrix( nrow=MAX_PHRED_SCORE, ncol=length(lanes))
	rownames( avgReadScores) <- 1:MAX_PHRED_SCORE
	colnames( avgReadScores) <- LaneIDs[lanes]
	for( i in 1:length(lanes)) {
		avgReadScores[ ,i] <- apply( scoreReads[ , , lanes[i]], 1, sum) / nLaneScores[ lanes[i]] * 100
	}

	avgTileReadScores <- array( 0, dim=c(MAX_PHRED_SCORE, length( tiles), length(lanes)), dimnames=list(
					1:MAX_PHRED_SCORE, TileIDs[tiles], LaneIDs[lanes]))
	for( i in 1:length(lanes)) {
	    lane <- lanes[i]
	    thisNscores <- nScores[ , lane]
	    for( j in 1:length(tiles)) {
		tile <- tiles[j]
		if( thisNscores[ tile] > 0) avgTileReadScores[ ,j,i] <- scoreReads[ , tile, lane] / thisNscores[ tile] * 100
	}}

	# trap any invalid values
	avgScores[ ! is.finite(avgScores)] <- 0
	avgTileScores[ ! is.finite(avgTileScores)] <- 0
	avgReadScores[ ! is.finite(avgReadScores)] <- 0
	avgTileReadScores[ ! is.finite(avgTileReadScores)] <- 0
	nBases[ ! is.finite(nBases)] <- 0
	nScores[ ! is.finite(nScores)] <- 0
	GrandAvgScores[ ! is.finite(GrandAvgScores)] <- 0

	# stash what we need to plot for later
	qualPlotFile <<- scoreFile
	qualPlotSampleID <<- scoreSample
	qualPlotAvgScores <<- avgScores
	qualPlotTileScores <<- avgTileScores
	qualPlotAvgReadScores <<- avgReadScores
	qualPlotTileReadScores <<- avgTileReadScores
	qualPlotBases <<- nBases
	qualPlotNscoresVec <<- nScores
	qualPlotGrandAvgScore <<- GrandAvgScores
	qualPlotLaneIDs <<- LaneIDs
	qualPlotTileIDs <<- TileIDs
	qualPlotReadType <<- "Total"
	return()
}


`reportCounts` <- function( verbose=TRUE) {

	nLaneCounts <- apply( nCounts, 2, sum)
	lanes <- which( nLaneCounts > 0)
	mList <- mRepList <- list()
	for( i in 1:length( lanes)) {
		lane <- lanes[i]
		m <- 100 * countSums[ lane, , , ] / nLaneCounts[lane]
		mList[[i]] <- m
		m <- 100 * countReps[ lane, , ] / nLaneCounts[lane]
		mRepList[[i]] <- m
	}
	names( mList) <- names(mRepList) <- LaneIDs[lanes]

	if (verbose) {
		cat( "\n\nN_Bad Reads ( non-ACGT bases detected): \n ", paste( "\tLane: ", LaneIDs[lanes], "\t", 
			apply( as.matrix(nBadReads[,lanes]), 2, sum), "\t", 
			formatC( (100 * (apply( as.matrix(nBadReads[,lanes]), 2, sum)/nLaneCounts)[lanes]), 
			digits=2, format="f"), "%", collapse="\n"))
	}

	nTileCounts <- apply( nCounts, 1, sum)
	tiles <- which( nTileCounts > 0)
	mTileList <- mTileRepList <- vector( mode="list", length=(length(lanes)*length(tiles)))

	for( i in 1:length( lanes)) {
	    lane <- lanes[i]
	    for( j in 1:length( tiles)) {
		tile <- tiles[j]
		thisN <- (i-1) * length(tiles) + j
		m <- 100 * countSums[ lane, tile, , ] / nCounts[tile, lane]
		mTileList[[ thisN ]] <- m
		names( mTileList)[thisN] <- paste( "Tile", TileIDs[tile], "_Lane", LaneIDs[lane], sep="")
		m <- 100 * countReps[ lane, tile, ] / nCounts[tile, lane]
		mTileRepList[[ thisN ]] <- m
		names( mTileRepList)[thisN] <- paste( "Tile", TileIDs[tile], "_Lane", LaneIDs[lane], sep="")
	}}

	# trap any invalid values
	for (k in 1:length(mList)) { m <- mList[[k]]; m[ ! is.finite(m)] <- 0; mList[[k]] <- m }
	for (k in 1:length(mTileList)) { m <- mTileList[[k]]; m[ ! is.finite(m)] <- 0; mTileList[[k]] <- m }
	for (k in 1:length(mRepList)) { m <- mRepList[[k]]; m[ ! is.finite(m)] <- 0; mRepList[[k]] <- m }
	for (k in 1:length(mTileRepList)) { m <- mTileRepList[[k]]; m[ ! is.finite(m)] <- 0; mTileRepList[[k]] <- m }
	nCounts[ ! is.finite(nCounts)] <- 0

	# stash what we need to plot
	qualPlotBaseCnts <<- mList
	qualPlotTileBaseCnts <<- mTileList
	qualPlotBaseReps <<- mRepList
	qualPlotTileBaseReps <<- mTileRepList
	qualPlotNcountsVec <<- nCounts
	return()
}

# end of local functions...


	# 4 lines per read
	readChunkSize <- chunkSize * 4
	nread <- 0
	hasMore <- TRUE
	firstPass <- TRUE
	ansTiles <- ansScores <- NULL

	repeat {
		if ( ! hasMore) break
		if ( ! is.null( maxReads) && nread >= maxReads) break
		cat( "\nBuffer Read..")
		chunk <- readLines( conIn, n=readChunkSize)
		if ( length( chunk) < 1) break
		if ( length( chunk) < readChunkSize) hasMore <- FALSE

		# in a FASTQ file, it takes 4 lines for each read;  1 is ID, 2nd is the seq, 4th is scores
		ptrs <- seq.int( 1, length(chunk), by=4)
		idsTxt <- chunk[ ptrs]
		seqTxt <- chunk[ ptrs+1]
		scoresTxt <- chunk[ ptrs+3]

		# allow dropping of some sequences...
		if ( ! is.null(dropPatterns)) {
			cat( " dropping..")
			readSEQfac <- factor( seqTxt)
			drops <- vector()
			for( pattern in dropPatterns) {
				drops <- base::append( drops, grep( pattern, levels(readSEQfac)))
			}
			picks <- which( !( tapply( 1:length(seqTxt), INDEX=readSEQfac, FUN=NULL) %in% drops))
			idsTxt <- idsTxt[ picks ]
			seqTxt <- seqTxt[ picks ]
			scoresTxt <- scoresTxt[ picks ]
		}

		nread <- nread + length( idsTxt)

		cat( "  N_reads: ", prettyNum( as.integer(nread), big.mark=","))

		cat( "  parse readID..")
		ans <- extractReadIDterms( idsTxt, readIDtype=readIDtype)
		laneNumber <- ans$Lane
		tileNumber <- ans$Tile

		# some older formats now have big tile numbers...
		if (firstPass) {
			if ( max( tileNumber, na.rm=T) > MAX_TILES) keepAllTiles <- FALSE
			firstPass <- FALSE
		}

		cat( "  tabulate..")
		if ( scoreType == "SolexaChar") {
			scoresTxt <- solexaCharScoresToPhredScores( scoresTxt)
		}
		qscores <- addScores( scoresTxt, laneNumber, tileNumber, scoreType)
		if ( is.null( qscores)) {
			cat( "\nMaximum Lanes/Tiles exceeded...  Quiting early...")
			return( list( "ScoreSummary"=ansScores, "TileSummary"=ansTiles))
		}
		addCounts( seqTxt, laneNumber, tileNumber)

		# try doing the plots every time thru...
		cat( "  plot..")
		reportScores( verbose=FALSE)
		reportCounts( verbose=FALSE)

		#saveFile <- paste( "qualPlotData_", basename(filein), ".Rd", sep="")
		saveFile <- statsFile
		save( qualPlotFile, qualPlotSampleID, qualPlotAvgScores, qualPlotTileScores, 
			qualPlotBases, qualPlotBaseCnts, qualPlotTileBaseCnts, qualPlotAvgReadScores, qualPlotTileReadScores,
			qualPlotNscoresVec, qualPlotNcountsVec, qualPlotGrandAvgScore, qualPlotLaneIDs, qualPlotTileIDs,
			qualPlotBaseReps, qualPlotTileBaseReps, qualPlotReadType, file=saveFile)

		if ( ! capabilities( "png") ) {
			cat( "\n Can not plot...  saving plot data for non-cluster computer...")
			cat( paste( "\n In R, enter:  'plotQualityScores( statsFile=\"", saveFile,"\")' \n", sep=""))
			ansTiles <- ansScores <- NULL
		} else {
			plotQualityTileVariance( asPNG=FALSE, statsFile=saveFile, plotPath=plotPath, baseOrder=baseOrder) 
			plotQualityTileVariance( asPNG=TRUE, statsFile=saveFile, plotPath=plotPath, baseOrder=baseOrder) 
			if ( pause > 0) Sys.sleep(pause)

			ansTiles <- plotQualityTileAverage( asPNG=FALSE, statsFile=saveFile, plotPath=plotPath, baseOrder=baseOrder) 
			ansTiles <- plotQualityTileAverage( asPNG=TRUE, statsFile=saveFile, plotPath=plotPath, baseOrder=baseOrder) 
			if ( pause > 0) Sys.sleep(pause)

			ansReps <- plotQualityBaseReps( asPNG=FALSE, statsFile=saveFile, plotPath=plotPath)
			ansReps <- plotQualityBaseReps( asPNG=TRUE, statsFile=saveFile, plotPath=plotPath)
			if ( pause > 0) Sys.sleep(pause)

			ansScores <- plotQualityScores( asPNG=FALSE, statsFile=saveFile, plotPath=plotPath, baseOrder=baseOrder)
			ansScores <- plotQualityScores( asPNG=TRUE, statsFile=saveFile, plotPath=plotPath, baseOrder=baseOrder)
			if ( pause > 0) Sys.sleep(pause)

		}

	}

	return( list( "ScoreSummary"=ansScores, "TileSummary"=ansTiles))
}

