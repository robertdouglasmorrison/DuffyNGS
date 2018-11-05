# plotReadStats.R -- functions to visualize read stats from either FASTQ or BAM files


`readStats` <- function( filein, sampleID, ...) {

	if ( regexpr( "\\.bam$", filein) > 0) {
		ans <- bamReadStats( filein, sampleID, ...)
	} else {
		ans <- fastqReadStats( filein, sampleID, ...)
	}
	return( ans)
}



`plotQualityBaseReps` <- function( asPNG=FALSE, statsFile=NULL, plotPath=".", laneNumber=NULL, tileNumber=NULL) {

	# allow re-plot...
	if ( ! is.null( statsFile)) load( statsFile)

	fname <- basename( qualPlotFile)
	sampleID <- qualPlotSampleID

	# number of plots based on number of lanes of reads
	myNscoresVec <- apply( qualPlotNscoresVec, 2, sum)
	lanes <- which( myNscoresVec > 0)
	nLanes <- length( lanes)
	allLaneNames <- qualPlotLaneIDs[lanes]

	# allow focused plots...
	if ( !is.null( laneNumber)) {
		who <- base::match( laneNumber, allLaneNames, nomatch=0)
		nLanes <- sum( who > 0)
		lanes <- laneNumber[ who > 0]
	}

	for( i in 1:nLanes) {
	    lane <- lanes[i]
	    laneText <- allLaneNames[lane]

	    # perhaps do a focused by tile look
	    byTile <- FALSE
	    myTilescoresVec <- qualPlotNscoresVec[ , lane]
	    tiles <- which( myTilescoresVec > 0)
	    allTileNames <- qualPlotTileIDs[tiles]
	    nTiles <- tiles <- 1
	    if ( ! is.null( tileNumber)) {
		who <- base::match( tileNumber, allTileNames, nomatch=0)
		nTiles <- sum( who > 0)
		tiles <- tileNumber[ who > 0]
		byTile <- TRUE
	    }

	    for( j in 1:nTiles) {
		tile <- tiles[j]
	    	tileText <- allTileNames[tile]

		plotname <- paste( sampleID, "_Lane",laneText,"_readRepeats.png",sep="")
		if ( byTile) plotname <- paste( sampleID, "_Lane",laneText,"_Tile",tileText,"_readRepeats.png",sep="")

		if ( asPNG) {
			xlim <- c(0, qualPlotBases*1.2)
			cexMain <- par( "cex.main") * 1.35
			cexLab <- par( "cex.axis") * 1.5
			cexNames <- par( "cex.axis") * 1.2
		} else {
			xlim <- c(0, qualPlotBases*1.35)
			cexMain <- par( "cex.main")
			cexLab <- par( "cex.axis") * 1.15
			cexNames <- par( "cex.axis")
		}
		if ( asPNG) {
			png( file=file.path( plotPath, plotname), width=800, height=600)
			par( mfcol=c(1,1))
			par("mai"=c(1.6, 1.22, 1.92, 0.40))
		} else {
			par( mfcol=c(1,1))
			par("mai"=c(1.6, 1.22, 1.92, 0.40))
		}


		thisCounts <- qualPlotBaseReps[[ i]]
		# the counts are per tile, so aggreate
		myCounts <- apply( thisCounts, 2, sum)

		if (byTile) {
			thisID <- paste( "Tile",tileText,"_Lane",laneText, sep="")
			who <- base::match( thisID, names( qualPlotTileBaseReps), nomatch=0)
			if ( who[1] == 0) next
			thisCounts <- qualPlotTileBaseReps[[ who[1] ]]
			myCounts <- thisCounts
		}

		mainText <- paste( "Adjacent Base Repeat Distribution: \nFile: ", 
				fname, "\nLane:", laneText, "      ", qualPlotReadType, "Reads: ", 
				prettyNum( as.integer( myNscoresVec[lane]),big.mark=","))

		barplot( myCounts, main=mainText, legend=FALSE, 
			xlim=xlim, ylim=c(0,100), 
			xlab="Base Position    (Cycle Number)", 
			ylab="Frequency of Adjacent Base Repeats", cex.lab=cexLab, cex.main=cexMain, 
			cex.names=cexNames, cex.axis=cexNames)

		legend( "topright", paste( c("Mean:  ", "StdDev:"), formatC( c(mean(myCounts),sd(myCounts)), 
				format="f", digits=2), sep="  "), cex=cexNames)

		#if ( ( ! asPNG) && (i < nLanes)) locator(1)
		if ( asPNG) dev.off()
	}}

	par( mfcol=c(1,1))
	par("mai"=c(1.02, 0.82, 0.82, 0.40))

	avgRepeats <- sapply( qualPlotBaseReps, function(x) { y <- apply( x, 2, sum); return( mean.default( y))})
	sdRepeats <- sapply( qualPlotBaseReps, function(x) { y <- apply( x, 2, sum); return( sd( y))})
	varRepeats <- sdRepeats * sdRepeats
	out <- list( "LaneNumbers"=allLaneNames, "TotalReads"=myNscoresVec[lanes], "RepeatAverages"=avgRepeats, 
			"RepeatVariances"=varRepeats)
	return( out)
}


`plotQualityScores` <- function( asPNG=FALSE, statsFile=NULL, plotPath=".", laneNumber=NULL, tileNumber=NULL,
			baseOrder="ACGTN") {

	# allow re-plot...
	if ( ! is.null( statsFile)) load( statsFile)

	fname <- basename( qualPlotFile)
	sampleID <- qualPlotSampleID

	# number of plots based on number of lanes of reads
	myNscoresVec <- apply( qualPlotNscoresVec, 2, sum)
	lanes <- which( myNscoresVec > 0)
	nLanes <- length( lanes)
	allLaneNames <- qualPlotLaneIDs[lanes]

	MAX_PHRED_SCORE <- getMaxPhredScore()

	# allow focused plots...
	if ( !is.null( laneNumber)) {
		who <- base::match( laneNumber, allLaneNames, nomatch=0)
		nLanes <- sum( who > 0)
		lanes <- laneNumber[ who > 0]
	}

	for( i in 1:nLanes) {
	    lane <- lanes[i]
	    laneText <- allLaneNames[lane]

	    # perhaps do a focused by tile look
	    byTile <- FALSE
	    myTilescoresVec <- qualPlotNscoresVec[ , lane]
	    tiles <- which( myTilescoresVec > 0)
	    allTileNames <- qualPlotTileIDs[tiles]
	    nTiles <- tiles <- 1
	    if ( ! is.null( tileNumber)) {
		who <- base::match( tileNumber, allTileNames, nomatch=0)
		nTiles <- sum( who > 0)
		tiles <- tileNumber[ who > 0]
		byTile <- TRUE
	    }

	    for( j in 1:nTiles) {
		tile <- tiles[j]
		tileText <- allTileNames[tile]

		plotname <- paste( sampleID, "_Lane",laneText,"_readStats.png",sep="")
		if ( byTile) plotname <- paste( sampleID, "_Lane",laneText,"_Tile",tileText,"_readStats.png",sep="")

		if ( asPNG) {
			xlim <- c(0, qualPlotBases*1.25)
			cexMain <- par( "cex.main") * 1.35
			cexLab <- par( "cex.axis") * 1.5
			cexNames <- par( "cex.axis") * 1.2
		} else {
			xlim <- c(0, qualPlotBases*1.30)
			cexMain <- par( "cex.main")
			cexLab <- par( "cex.axis") * 1.15
			cexNames <- par( "cex.axis")
		}
		if ( asPNG) {
			png( file=file.path( plotPath, plotname), width=900, height=900)
			par( mfcol=c(3,1))
			par("mai"=c(0.82, 1.02, 0.82, 0.40))
		} else {
			par( mfcol=c(3,1))
			par("mai"=c(0.82, 1.02, 0.62, 0.40))
		}


		thisScores <- qualPlotAvgScores[ , i]
		thisReads <- qualPlotAvgReadScores[ ,  i]
		thisCounts <- qualPlotBaseCnts[[ i]]
		# the counts are per tile, so aggreate
		myCounts <- apply( thisCounts, 2:3, sum)

		if (byTile) {
			thisScores <- qualPlotTileScores[ , j, i]
			thisReads <- qualPlotTileReadScores[ , j, i]
			thisID <- paste( "Tile",tileText,"_Lane",laneText, sep="")
			who <- base::match( thisID, names( qualPlotTileBaseCnts), nomatch=0)
			if ( who[1] == 0) next
			thisCounts <- qualPlotTileBaseCnts[[ who[1] ]]
			myCounts <- thisCounts
		}

		names( thisScores) <- 1:length( thisScores)
		mainText <- paste( "Base Quality Scores: \nFile: ", fname, "\nLane:", laneText, 
				"      ", qualPlotReadType, "Reads: ", prettyNum( as.integer( myNscoresVec[lane]),big.mark=","))
		if ( byTile) mainText <- paste( "Base Call Quality Scores: \nFile: ", fname, "\nLane:", 
				laneText, "      Tile:", tileText, "      ", qualPlotReadType, "Reads: ", 
				prettyNum( as.integer(qualPlotNcountsVec[tile,lane]),big.mark=","))

		borderColor <- if ( length( thisScores) < 75) "black" else NA
		myYmax <- max( thisScores, min(40,MAX_PHRED_SCORE), na.rm=T)

		barplot( thisScores, main=mainText, col="lightblue", 
			xlim=xlim, ylim=c(0,myYmax), xlab="Base Position    (Cycle Number)", 
			ylab="Avg Q Score per Cycle",
			cex.lab=cexLab, cex.main=cexMain, cex.names=cexNames, border=borderColor)

		par("mai"=c(0.82, 1.02, 0.42, 0.40))

		ord <- getBaseOrder( baseOrder)
		barplot( myCounts[ ord, ], main="Base Call Nucleotides:", legend=FALSE, 
			col=c(2,4,"orange",3,"brown")[ord], xlim=xlim, ylim=c(0,100), 
			xlab="Base Position    (Cycle Number)", 
			ylab="Percent by Base Type", cex.lab=cexLab, cex.main=cexMain, cex.names=cexNames, border=borderColor)
		legend( "topright", rev( c( "A", "C", "G", "T", "N")[ord]), 
			fill=rev( c(2,4,"orange",3,"brown")[ord]), cex=1.25)

		yLim <- c( 0, max( thisReads)*1.2)
		barplot( thisReads, main="Read Quality Scores:", ylim=yLim,
			col="lightgreen", ylab="Percent of All Reads", xlab="Average Q Score per Read",
			cex.lab=cexLab, cex.main=cexMain, cex.names=cexNames, border=borderColor)
		legend( "topleft", paste( "Avg Q Score = ", formatC( mean( qualPlotAvgScores[ ,i]), 
			digits=2, format="f")), cex=1.3, bg="white")

		#if ( ( ! asPNG) && (i < nLanes)) locator(1)
		if ( asPNG) dev.off()
	}}

	par( mfcol=c(1,1))
	par("mai"=c(1.02, 0.82, 0.82, 0.40))

	avgScore <- apply( qualPlotAvgScores, 2, mean.default)
	sdScore <- apply( qualPlotAvgScores, 2, sd)
	varScore <- sdScore * sdScore
	baseFractions <- lapply( qualPlotBaseCnts, function(x) { 
				y <- apply( x, 2:3, sum); 
				return( apply( y, 1,mean.default))
			})
	out <- list( "LaneNumbers"=allLaneNames, "TotalReads"=myNscoresVec[lanes], "ScoreAverages"=avgScore, 
			"ScoreVariances"=varScore, "BaseComposition"=baseFractions)
	return( out)
}


`plotQualityTileVariance` <- function( asPNG=FALSE, statsFile=NULL, plotPath=".", baseOrder="ACGTN") {

	# allow re-plot...
	if ( ! is.null( statsFile)) load( statsFile)

	fname <- basename( qualPlotFile)
	sampleID <- qualPlotSampleID

	# number of plots based on number of lanes of reads
	myNscoresVec <- apply( qualPlotNscoresVec, 2, sum)
	lanes <- which( myNscoresVec > 0)
	nLanes <- length( lanes)
	allLaneNames <- qualPlotLaneIDs[lanes]

	Nbases <- qualPlotBases
	MAX_PHRED_SCORE <- getMaxPhredScore()

	for( i in 1:nLanes) {
		lane <- lanes[i]
		laneText <- allLaneNames[lane]

		# number of tiles based on how many have data
		myTilescoresVec <- qualPlotNscoresVec[ , lane]
		tiles <- which( myTilescoresVec > 0)
		allTileNames <- qualPlotTileIDs[tiles]
		Ntiles <- length( tiles)

		plotname <- paste( sampleID, "_Lane",laneText,"_tileVariance.png",sep="")

		if ( asPNG) {
			xlim <- c(0, qualPlotBases*1.25)
			cexMain <- par( "cex.main") * 1.35
			cexLab <- par( "cex.axis") * 1.5
			cexNames <- par( "cex.axis") * 1.2
		} else {
			xlim <- c(0, qualPlotBases*1.30)
			cexMain <- par( "cex.main")
			cexLab <- par( "cex.axis") * 1.15
			cexNames <- par( "cex.axis")
		}
		if ( asPNG) {
			png( file=file.path( plotPath, plotname), width=900, height=900)
			par( mfcol=c(3,1))
			par("mai"=c(0.82, 1.02, 0.82, 0.40))
		} else {
			par( mfcol=c(3,1))
			par("mai"=c(0.82, 1.02, 0.72, 0.40))
		}

		# make a variance call on the base type counts of each tile...
		allPcts <- array( 0, dim=c( Ntiles, 5, Nbases))
		for ( k in 1:length( qualPlotTileBaseCnts)) {
			thisName <- names( qualPlotTileBaseCnts)[k]
			if ( regexpr( paste( "_Lane",laneText,sep=""), thisName) < 1) next
			thisCounts <- qualPlotTileBaseCnts[[ k ]]
			if ( is.null( thisCounts)) next
			myLane <- sub( "(Tile[0-9]+_Lane)([0-9]+)", "\\2", thisName)
			if ( myLane != laneText) next
			myTile <- sub( "(Tile)([0-9]+)(_Lane.*)", "\\2", thisName)
			whereTile <- base::match( myTile, allTileNames, nomatch=0)
			if ( whereTile < 1) next
			if ( whereTile > Ntiles) next
			allPcts[ whereTile, , ] <- thisCounts
		}
		myVars <- matrix( 0, nrow=5, ncol=Nbases)
		for ( kType in 1:5 )
		for ( kLoc in 1:Nbases) {
			x <- allPcts[ , kType, kLoc]
			mysd <- sd( x)
			myVars[ kType, kLoc] <- mysd * mysd
		}
		colnames( myVars) <- 1:Nbases

		thisZeros <- rep(0, times=Nbases)
		names( thisZeros) <- 1:Nbases
		mainText <- paste( "Base Quality Score:   Tile Variance \nFile: ", fname, "\nLane:", 
				laneText, "      ", qualPlotReadType, "Reads: ", 
				prettyNum( as.integer( myNscoresVec[lane]),big.mark=","))

		myYmax <- max( qualPlotTileScores[ , ,i], min(40,MAX_PHRED_SCORE), na.rm=T)

		ats <- barplot( thisZeros, main=mainText,col="lightblue", 
			xlim=xlim, ylim=c(0,myYmax), xlab="Base Position     (Cycle Number)", 
			ylab="Avg Q Score per Cycle",
			cex.lab=cexLab, cex.main=cexMain, cex.names=cexNames)
		X <- ats
		for( k in 1:Ntiles) {
			myTile <- tiles[k]
			Y <- qualPlotTileScores[ , myTile, i]
			lines( X, Y, col="lightblue", lwd=2, lty=1)
		}

		par("mai"=c(0.82, 1.02, 0.42, 0.40))
		thisYlim <- c(0,100)
		bigY <- max( apply( myVars,2,sum))
		if ( !is.na(bigY) && bigY > 100) thisYlim <- c( 0, bigY*0.9)

		borderColor <- if ( ncol( myVars) < 75) "black" else NA

		ord <- getBaseOrder( baseOrder)
		barplot( myVars[ord, ], main="Base Call Nucleotides:    Tile Variance ", legend=FALSE, 
			col=c(2,4,"orange",3,"brown")[ord], xlim=xlim, ylim=thisYlim, 
			xlab="Base Position     (Cycle Number)", 
			ylab="Variance % by Base Type", cex.lab=cexLab, cex.main=cexMain, 
			cex.names=cexNames, border=borderColor)
		legend( "topright", rev( c( "A", "C", "G", "T", "N")[ord]), 
			fill=rev( c(2,4,"orange",3,"brown")[ord]), cex=1.25)

		yLim <- c( 0, max( qualPlotTileReadScores[ , , i])*1.2)
		thisZeros <- rep(0, times=MAX_PHRED_SCORE)
		names( thisZeros) <- 1:MAX_PHRED_SCORE
		ats <- barplot( thisZeros, main="Read Quality Score:    Tile Variance ", ylim=yLim,
			col="lightgreen", ylab="Percent of All Reads", xlab="Average Q Score per Read",
			cex.lab=cexLab, cex.main=cexMain, cex.names=cexNames)
		X <- ats
		for( k in 1:Ntiles) {
			myTile <- tiles[k]
			Y <- qualPlotTileReadScores[ , myTile, i]
			lines( X, Y, col="lightgreen", lwd=2, lty=1)
		}

		legend( "topleft", paste( "Avg Q Score = ", formatC( mean( qualPlotAvgScores[ ,i]), 
			digits=2, format="f")), cex=1.3, bg="white")

		if ( asPNG) dev.off()
	}

	par( mfcol=c(1,1))
	par("mai"=c(1.02, 0.82, 0.82, 0.40))
	return()
}


`plotQualityTileAverage` <- function( asPNG=FALSE, statsFile=NULL, plotPath=".", baseOrder="ACGTN") {

	# allow re-plot...
	if ( ! is.null( statsFile)) load( statsFile)

	fname <- basename( qualPlotFile)
	sampleID <- qualPlotSampleID

	# number of plots based on number of lanes of reads
	myNscoresVec <- apply( qualPlotNscoresVec, 2, sum)
	lanes <- which( myNscoresVec > 0)
	nLanes <- length( lanes)
	allLaneNames <- qualPlotLaneIDs[lanes]

	Nbases <- qualPlotBases
	MAX_PHRED_SCORE <- getMaxPhredScore()

	for( i in 1:nLanes) {
		lane <- lanes[i]
		laneText <- allLaneNames[lane]

		# number of tiles based on how many have data
		myTilescoresVec <- qualPlotNscoresVec[ , lane]
		tiles <- which( myTilescoresVec > 0)
		allTileNames <- qualPlotTileIDs[tiles]
		Ntiles <- length( tiles)
		MaxTile <- max( tiles)
		maxTileNames <- paste( "unused", 1:MaxTile, sep="_")
		maxTileNames[tiles] <- allTileNames

		plotname <- paste( sampleID, "_Lane",laneText,"_tileAverage.png",sep="")

		if ( asPNG) {
			xlim <- c(0, Ntiles*1.25)
			cexMain <- par( "cex.main") * 1.35
			cexLab <- par( "cex.axis") * 1.5
			cexNames <- par( "cex.axis") * 1.2
		} else {
			xlim <- c(0, Ntiles*1.25)
			cexMain <- par( "cex.main")
			cexLab <- par( "cex.axis") * 1.15
			cexNames <- par( "cex.axis")
		}
		if ( asPNG) {
			png( file=file.path( plotPath, plotname), width=900, height=900)
			par( mfcol=c(3,1))
			par("mai"=c(0.62, 1.02, 0.82, 0.40))
		} else {
			par( mfcol=c(3,1))
			par("mai"=c(0.62, 1.02, 0.72, 0.40))
		}

		# tweak for when there are very few tiles...
		if ( Ntiles < 10) {
			myWidth <- 0.8
			mySpace <- 0.2
			xlim <- c(0, (10)*1.35)
		} else {
			myWidth <- 1
			mySpace <- NULL
		}

		# make a average call on the base type counts of each tile...
		allPcts <- array( 0, dim=c( MaxTile, 5, Nbases))
		for ( k in 1:length( qualPlotTileBaseCnts)) {
			thisName <- names( qualPlotTileBaseCnts)[k]
			if ( regexpr( paste( "_Lane",laneText,sep=""), thisName) < 1) next
			thisCounts <- qualPlotTileBaseCnts[[ k ]]
			if ( is.null( thisCounts)) next
			myLane <- sub( "(Tile[0-9]+_Lane)([0-9]+)", "\\2", thisName)
			if ( myLane != laneText) next
			myTile <- sub( "(Tile)([0-9]+)(_Lane.*)", "\\2", thisName)
			whereTile <- base::match( myTile, maxTileNames, nomatch=0)
			if ( whereTile < 1) next
			if ( whereTile > MaxTile) next
			allPcts[ whereTile, , ] <- thisCounts
		}
		myAvgs <- matrix( 0, nrow=5, ncol=MaxTile)
		for ( kType in 1:5 )
		for ( kTile in 1:MaxTile) {
			x <- allPcts[ kTile, kType, ]
			myavg <- mean.default( x)
			myAvgs[ kType, kTile] <- myavg
		}
		colnames( myAvgs) <- maxTileNames

		# the average score for each tile
		thisScores <- rep(0, times=MaxTile)
		for ( k in 1:Ntiles) {
			tileNum <- tiles[k]
			#Y <- qualPlotTileScores[ , k, i]
			Y <- qualPlotTileScores[ , tileNum, i]
			thisScores[ tileNum] <- mean.default(Y)
		}
		names( thisScores) <- maxTileNames

		# the Hi-Seq puts out some very large Tile numbers
		if ( MaxTile < 130) {
			showTiles <- 1:MaxTile
		} else {
			showTiles <- tiles
		}

		mainText <- paste( "Base Quality Score per Tile: \nFile: ", fname, "\nLane:", laneText, 
				"      ", qualPlotReadType, "Reads: ", prettyNum( as.integer( myNscoresVec[lane]),big.mark=","))
		
		borderColor <- if ( length( thisScores) < 75) "black" else NA
		myYmax <- max( thisScores[showTiles], min(40,MAX_PHRED_SCORE), na.rm=T)

		ats <- barplot( thisScores[showTiles], main=mainText,col="lightblue", xlim=xlim, ylim=c(0,myYmax), xlab="Tile Number",
				width=myWidth, space=mySpace, ylab="Avg Q Score per Tile", cex.lab=cexLab, cex.main=cexMain, 
				cex.names=cexNames, border=borderColor)
		hasScrs <- which( thisScores > 0)
		if ( length( hasScrs) > 0) {
			spt_avg <- mean.default( thisScores[ hasScrs])
			spt_sd <- sd( thisScores[ hasScrs])
		} else {
			spt_avg <- spt_sd <- NA
		}

		# the ACGT base mix
		par("mai"=c(0.82, 1.02, 0.42, 0.40))
		thisYlim <- c(0,100)
		ord <- getBaseOrder( baseOrder)
		barplot( myAvgs[ord, showTiles, drop=FALSE], main="Base Call Nucleotides per Tile:", legend=FALSE, 
			col=c(2,4,"orange",3,"brown")[ord], xlim=xlim, 
			ylim=thisYlim, xlab="Tile Number", ylab="Percent by Base Type", cex.lab=cexLab, 
			width=myWidth, space=mySpace, cex.main=cexMain, cex.names=cexNames, border=borderColor)
		legend( "topright", rev( c( "A", "C", "G", "T", "N")[ord]), fill=rev( c(2,4,"orange",3,"brown")[ord]), cex=1.25)

		# the number of reads
		yLim <- c( 0, max( myTilescoresVec)*1.2)
		thisCounts <- myTilescoresVec[ 1:MaxTile]
		names( thisCounts) <- maxTileNames
		ats <- barplot( thisCounts[ showTiles], main="Number of Reads per Tile:", ylim=yLim, xlim=xlim,
				col="lightgreen", ylab="Read Count", xlab="Tile Number", width=myWidth, space=mySpace,
				cex.lab=cexLab, cex.main=cexMain, cex.names=cexNames, border=borderColor)

		hasCnts <- which( thisCounts > 0)
		rpt_avg <- Y <- mean.default( thisCounts[ hasCnts])
		rpt_sd <- sd( thisCounts[ hasCnts])
		lines( c(0,MaxTile*1.25), rep(Y,2), lwd=2, lty=3, col=2)
		legend( "topleft", paste( "Mean Reads/Tile = ", round(Y)), lty=3, lwd=3, col=2, cex=1.2, bg="white")

		if ( asPNG) dev.off()
	}

	par( mfcol=c(1,1))
	par("mai"=c(1.02, 0.82, 0.82, 0.40))

	out <- list( "ScorePerTile_Avg"=spt_avg, "ScorePerTile_StdDev"=spt_sd, "ReadsPerTile_Avg"=rpt_avg, 
			"ReadsPerTile_StdDev"=rpt_sd)
	return( out)
}

