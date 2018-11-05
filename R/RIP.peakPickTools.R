# RIP.peakPickTools.R


# find all peaks in an object with X, Y elements
# specialized for RIP-seq data...
`which.in` <- function( a, b) base::which( base::match( a, b, nomatch=0) > 0)


`findRIPpeaks` <- function( tbl, peak.type=c("auto","gaussian","gumbel","pulse"), 
			cutoff.medians=3, min.width=50, max.width=1000, plot.file=NULL, 
			roc.file=NULL, visualize=interactive(), controlPeaksFile=NULL) {
				
	peak.type <- match.arg( peak.type)

	# given a data object with X, Y...
	allX <- as.integer( tbl$x)
	allY <- as.numeric( tbl$y)
	# force the data to be continuous, with no gaps in X
	Xrange <- range( allX)
	NX <- diff( Xrange) + 1
	if ( NX != length( allX)) {
		cat( "  Expand sparse X..")
		newX <- Xrange[1] : Xrange[2]
		newY <- rep.int( 0, NX)
		where <- allX - Xrange[1] + 1
		newY[where] <- allY
		allX <- newX
		allY <- newY
		rm( newX, newY, where)
	}
	medianY <- median( allY[ allY > 0])
	minimumY <- min( allY)
	cutoffY <- medianY * cutoff.medians
	strand <- tbl$strand

	# the universe of possible peak tops is everyone higher than the cutoff
	whoPossiblePeaks <- which( allY > cutoffY)
	whoExtrema <- which.extrema( allY, canonical.width=min.width*2)
	whoPossiblePeaks <- intersect( whoPossiblePeaks, whoExtrema)
	showPossibles <- whoPossiblePeaks

	# we will visit them in order by height
	ord <- order( allY[ whoPossiblePeaks], decreasing=TRUE)
	whoPossiblePeaks <- whoPossiblePeaks[ord]

	outM <- outSD <- outH <- outV <- outR <- outStart <- outStop <- outFloor <- vector()
	outType <- outStat <- outDrift <- outTimes <- vector()
	nout <- 0

	# include a neighborhood about 4x bigger than the peak's expected footprint
	mid.width <- mean( c( min.width, max.width))
	extra <- as.integer( mid.width * 4)

	# repeatedly find the tallest peak, until all possible are gone
	useY <- allY
	myplotfile <- plot.file
	visualize <- ( visualize && (! is.null( plot.file)))
	if ( visualize) checkX11( width=10, height=7, bg="white")

	while ( length( whoPossiblePeaks) > 0) {

		who <- whoPossiblePeaks[1]
		centerX <- allX[ who]
		centerY <- useY[ who]
		#cat( "\nDebug: who,X,Y,extra: ", who, centerX, centerY, extra)
		doPlot <- ( visualize && (nout < 10))
		if (doPlot) {
			myplotfile <- sub( "png$", paste( strand, (nout+1), "png",sep="."), plot.file)
		} else {
			myplotfile <- NULL
		}

		# the region to look at is a neighborhood around the center
		myPts <- max( 1, (who-extra)) : min( length(allX), (who+extra))
		ans <- findOneBestRIPpeak( allX[myPts], useY[myPts], centerX, centerY, cutoff=cutoffY, 
					min.width=min.width, max.width=max.width, peak.type=peak.type, 
					extrema.pts=allX[showPossibles], plot.file=myplotfile,
					visualize=visualize)
		if ( is.null( ans)) {
			whoPossiblePeaks <- setdiff( whoPossiblePeaks, who)
			next
		}

		nout <- nout + 1
		outM[nout] <- myX <- ans$center
		outSD[nout] <- mySD <- ans$width
		outH[nout] <- maxY <- ans$height
		outFloor[nout] <- myFloor <- ans$floor
		outV[nout] <- ans$volume
		outR[nout] <- ans$residual
		outDrift[nout] <- myDrift <- centerX - myX
		outTimes[nout] <- (myFloor + maxY) / medianY
		fitX <- ans$x
		fitY <- ans$y
		outStart[nout] <- myStart <- min( fitX)
		outStop[nout] <- myStop <- max( fitX)
		outType[nout] <- ans$type

		outStat[nout] <- "good"
		keep <- TRUE
		if ( is.infinite(mySD) || abs(mySD) > max.width*2) {
			# this exception needs to be trapped
			outStat[nout] <- "badWidth"
			outSD[nout] <- newSD <- max.width*2 * sign(mySD)
			bigtail <- abs(newSD * 3)
			myLeft <- round( myX - bigtail)
			myRight <- round( myX + bigtail)
			if ( is.infinite( myStart) || myStart < myLeft) outStart[nout] <- myLeft
			if ( is.infinite( myStop) || myStop > myRight) outStop[nout] <- myRight
			if ( outStart[nout] >= outStop[nout]) keep <- FALSE
		}
		if ( abs(myDrift) > max.width*1) outStat[nout] <- "badDrift"
		if ( abs(myDrift) > max.width*4) keep <- FALSE
		if ( myX < outStart[nout] | myX > outStop[nout]) outStat[nout] <- "badMean"
		theseX <- sort( allX[myPts])
		if ( myX < theseX[1] | myX > theseX[length(theseX)]) keep <- FALSE
		if ( maxY < medianY/8) outStat[nout] <- "badHeight"
		if ( maxY < 0) keep <- FALSE
		if ( (maxY+myFloor) < medianY) outStat[nout] <- "badHeight"
		if ( ! keep) {
			nout <- nout - 1
			whoPossiblePeaks <- setdiff( whoPossiblePeaks, who)
			next
		}

		if ( nout %% 200 == 0) {
			cat( "\n", strand, nout, myX, mySD, maxY, 
					"N_ToTest:", length(whoPossiblePeaks), "   ")
		}

		if ( outStat[nout] != "good") {
			whoPossiblePeaks <- setdiff( whoPossiblePeaks, who)
			next
		}

		# when we get a good peak, remove that fitted peak from the data itself
		#whoCovered <- which( allX %in% fitX)
		whoCovered <- which.in( allX, fitX)
		if ( length( whoCovered) < 1) {
			whoPossiblePeaks <- setdiff( whoPossiblePeaks, who)
			next
		}
		#whoFrom <- which( fitX %in% allX)
		whoFrom <- which.in( fitX, allX)
		
		# let's try an average of what's there with our model
		was <- useY[whoCovered]
		now <- fitY[ whoFrom]
		blend <- (was + now) / 2
		useY[whoCovered] <- useY[whoCovered] - blend
		wentNeg <- which( useY[whoCovered] < medianY)
		useY[ whoCovered[ wentNeg]] <- medianY

		# remove the locations covered by this peak from the possible places
		whoPossiblePeaks <- setdiff( whoPossiblePeaks, c( who, whoCovered))

		# every so often, re-assess whos the tallest left to visit
		if ( nout %% 20 == 0) {
			ord <- order( useY[ whoPossiblePeaks], decreasing=TRUE)
			whoPossiblePeaks <- whoPossiblePeaks[ord]
		}
		showPossibles <- whoPossiblePeaks
	}

	# do a bit of rounding...
	outM <- round( outM, digits=2)
	outSD <- round( outSD, digits=2)
	outH <- round( outH, digits=2)
	outFloor <- round( outFloor, digits=2)
	outStart <- round( outStart, digits=0)
	outStop <- round( outStop, digits=0)
	outV <- round( outV, digits=0)
	outR <- round( outR, digits=0)
	outDrift <- round( outDrift, digits=2)
	outTimes <- round( outTimes, digits=2)


	out <- data.frame( "center"=outM, "width"=outSD, "height"=outH, "floor"=outFloor, 
			"start"=outStart, "stop"=outStop, "volume"=outV, "residual"=outR,
			"drift"=outDrift, "type"=outType, "status"=outStat, "times"=outTimes,
			stringsAsFactors=F)
	rownames(out) <- 1:nrow(out)


	# we can give them all scores
	scoreAns <- scoreRIPpeaks( out, min.width=min.width, max.width=max.width, cutoff=medianY, strand=strand)
	out <- cbind( out, "score"=as.numeric(scoreAns$score))

	# use control peaks to assign P-values to the scores...
	# if given a file, use that...   else draw randomly from this file
	needRandomPeaks <- TRUE
	if ( ! is.null( controlPeaksFile)) {
		if ( ! file.exists( controlPeaksFile)) {
			cat( "\n'Control Peaks File' not found: ", controlPeaksFile)
			cat( "\nDefaulting to random peaks...")
		} else {
			ctrlTbl <- read.csv( controlPeaksFile, as.is=T)
			ctrlTbl <- subset( ctrlTbl, Class == strand)
			if ( nrow(ctrlTbl) > 0) {
				ctrlScore <- ctrlTbl$Score
				needRandomPeaks <- FALSE
				# trim out any NA controls
				ctrlScore <- ctrlScore[ ! is.na(ctrlScore)]
				Nctrl <- length(ctrlScore)
				cat( "\nPre-defined control peaks..   N =", Nctrl, "\n")
			} else {
				cat( "\nControl Score column not in control peaks file: ", wantedCol)
			}
		}
	}

	# get the population of random peaks as a negative control
	if (needRandomPeaks) {
		randPks <- randomRIPpeaks( tbl, N=1000, whoExtrema=whoExtrema, cutoff.medians=cutoff.medians, 
					min.width=min.width, max.width=max.width, strand=strand)
		ctrlScore <- randPks$score
		Nctrl <- length(ctrlScore)
	}

	# first measure of a cutoff is the 95th % of the controls
	quantCut <- quantile( ctrlScore, 0.95, na.rm=T)

	# use ROC to get another measure of cutoff
	NctrlUse <- min( Nctrl, 1000)
	isGood <- which( out$status == "good")
	NgoodUse <- min( length(isGood), 250)
	visualizeROC <- ( visualize && !is.null(roc.file))
	#cat( "\n", strand, visualize, roc.file, visualizeROC)
	rocCut <- peakScoreROC( goodScores=out$score[ isGood[ 1:NgoodUse]],
				badScores=ctrlScore[1:NctrlUse], label="RIP Peak Scores",
				visualize=visualizeROC)
	useCut <- max( c( rocCut, quantCut))
	cat( "\n\nPeak Score Thresholds:", strand, "\n  ROC:  ", rocCut, "\n  95% = ", 
		quantCut, "\n  Using:", useCut, "\n")
	if ( visualizeROC) {
		roc.file <- sub( "txt$", paste( strand, "ROC.png", sep="."), roc.file)
		dev.print( png, roc.file, width=700, height=720, bg='white')
	}

	# pvalue is based on percent of randoms that were better
	pval <- sapply( 1:nrow(out), function(ip) {
			n <- sum( ctrlScore >= out$score[ip])
			return( n / Nctrl)
		})
	out$p.value <- pval
	out$status[ out$score < useCut] <- "failScore"

	# final order by P-value
	ord <- order( out$p.value, -(out$score))
	out <- out[ ord, ]
	rownames(out) <- 1:nrow(out)
	
	return( list("peaks"=out, "scoreDetails"=scoreAns))
}


# find the best peak in the interval of data, given a suggested extrema center

`findOneBestRIPpeak` <- function( x, y, centerX=NULL, centerY=NULL, cutoff=median(y), min.width=50, max.width=1000,
				peak.type=c("auto","gaussian","gumbel","pulse"), 
				plot.file=NULL, extrema.pts=NULL, fit.floor=TRUE, 
				visualize=interactive()) {

	peak.type <- match.arg( peak.type)
	if ( is.null(centerX)) centerX <- x[ which.max(y)]
	if ( ! is.finite(centerX) || centerX < 1) return(NULL)
	if ( is.null(centerY)) centerY <- max(y, na.rm=T)
	if ( ! is.finite( centerY)) return(NULL)
	medianY <- median(y[y > 0], na.rm=T)

	visualize <- ( visualize && !is.null(plot.file))
	if (visualize) {
		xlim <- quantile(x, c(0.2,0.8), na.rm=T)
		if ( xlim[1] < 1000) xlim[1] <- 1
		ylim <- (range(y) * c(0.85,1.15))
		plot( x, y, main="Automatic Peak Finding", type="l", col=1, lwd=2, xlim=xlim, ylim=ylim)
		if ( ! is.null( extrema.pts)) {
			canSee <- which( x %in% extrema.pts)
			points( x[canSee], y[canSee], col=2)
			lines( xlim+c(-100,100), rep( medianY,2), lty=3, col=1, lwd=1)
		}
	}

	# try to estimate the starting width suggestion by walking part way down the peak
	width.step <- round( min.width * 0.1)
	if ( width.step < 2) width.step <- 2
	widleft <- widright <- width.step 
	gotleft <- gotright <- FALSE
	tooBig <- max.width * 1
	mycutoff <- max( cutoff, centerY*0.5)
	#if ( any( is.na( c( centerX, centerY, medianY, mycutoff)))) return( NULL)

	# step outward till we have enough info to know a starting width guess
	repeat {
		areaXs <- (centerX - widleft) : centerX
		#myYleft <- y[ which( x %in% areaXs)]
		myYleft <- y[ which.in( x, areaXs)]
		if ( any( myYleft < mycutoff)) {
			gotleft <- TRUE
		} else {
			widleft <- widleft + width.step
			if ( widleft > tooBig) gotleft <- TRUE
		}
		areaXs <- centerX : (centerX + widright)
		#myYright <- y[ which( x %in% areaXs)]
		myYright <- y[ which.in( x, areaXs)]
		if ( any( myYright < mycutoff)) {
			gotright <- TRUE
		} else {
			widright <- widright + width.step
			if ( widright > tooBig) gotright <- TRUE
		}
		if ( gotleft && gotright) break
	}
	centerX <- round( centerX + (widright - widleft)/2)
	mywidth <- round( (widleft + widright)/2)
	lowbar <- min( c( myYleft, myYright))
	if ( visualize) rect( centerX-mywidth, lowbar, centerX+mywidth, centerY, border="brown",
				lwd=2, lty=3)

	# OK, we have a rough idea of 1 SD, gather the data in that rough neighborhood 
	neighborWidth <- as.integer( max.width * 2)
	areaXs <- (centerX - neighborWidth) : (centerX + neighborWidth)
	#use <- which( x %in% areaXs)
	use <- which.in( x, areaXs)
	myXs <- x[ use]
	myYs <- y[ use]

	# see if we can clean the far edges of the neighborhood to better estimate the true background
	nx <- length( myXs)
	if ( sum( myXs > (centerX+min.width*3)) > (nx * 0.33)) {
		myYs <- cleanTailRegion( myXs, myYs, centerX, mywidth, "right")
	}
	if ( sum( myXs < (centerX-min.width*3)) > (nx * 0.33)) {
		myYs <- cleanTailRegion( myXs, myYs, centerX, mywidth, "left")
	}
	if ( visualize) {
		lines( myXs, myYs, col=6, lwd=2)
		lines( x, y, col=1, lwd=2)
	}

	if ( peak.type %in% c( "gumbel", "auto")) {
		ansGumbel <- fit.gumbel( myXs, myYs, fit.floor=fit.floor, start.center=centerX, 
				start.width=mywidth, start.height=centerY, start.floor=medianY)
		residualGumbel <- sum( abs( ansGumbel$residual))
		ans2 <- fit.gumbel( myXs, myYs, fit.floor=fit.floor, start.center=centerX, 
				start.width=-mywidth, start.height=centerY, start.floor=medianY)
		residual2 <- sum( abs( ans2$residual))
		if ( residual2 < residualGumbel) {
			ansGumbel <- ans2
			residualGumbel <- residual2
		}
		ans <- ansGumbel
		model <- "gumbel"
	}
	if ( peak.type %in% c( "gaussian", "auto")) {
		ansGauss <- fit.gaussian( myXs, myYs, fit.floor=fit.floor, start.center=centerX, 
				start.width=mywidth, start.height=centerY, start.floor=medianY)
		residualGauss <- sum( abs( ansGauss$residual))
		ans <- ansGauss
		model <- "gaussian"
	}
	if ( peak.type %in% c( "pulse", "auto")) {
		ansPulse <- fit.pulse( myXs, myYs, fit.floor=fit.floor, start.center=centerX, 
				start.width=mywidth, start.height=centerY, start.floor=medianY)
		residualPulse <- sum( abs( ansPulse$residual))
		ans <- ansPulse
		model <- "pulse"
	}
	if ( peak.type == "auto") {
		model <- "gaussian"
		ans <- ansGauss
		residualBest <- residualGauss
		if ( residualGumbel < residualBest) {
			ans <- ansGumbel
			residualBest <- residualGumbel
			model <- "gumbel"
		}
		if ( residualPulse < residualBest) {
			ans <- ansPulse
			residualBest <- residualPulse
			model <- "pulse"
		}
	}
	validFit <- drawable <- TRUE

	if (visualize) {
		lines( myXs, ans$y, col=4, lwd=2)
		text( myXs[ which.max( ans$y)], max( c( ans$y, centerY)), model, pos=3, font=2)
	}

	# package up the results
	floor <- if ( fit.floor) ans$floor else 0
	if ( floor < 0) floor <- 0
	fitX <- as.integer( myXs)
	fitY <- ans$y - floor
	maxY <- max( fitY)
	if ( ! any( fitY > 0)) validFit <- FALSE

	# decide which points are 'really in' the curve
	keep <- which( fitY >= (maxY *.05))
	if (length( keep) < 3) drawable <- FALSE
	fitX <- fitX[keep]
	fitY <- fitY[keep]
	yIn <- myYs[keep]

	# did the peak completely leave the area?
	spreadX <- diff( range(x))
	if ( abs( ans$center - centerX) > spreadX) validFit <- FALSE

	volume <- sum( fitY)
	residual <- sum( abs( yIn - (fitY+floor)))

	if (visualize) {
		if (drawable) {
			lines( fitX, fitY+floor, col=3, lwd=3)
			lines( c( rep(fitX[1],2), rep(fitX[length(fitX)],2)), 
			 	(c(fitY[1],0,0,fitY[length(fitY)]) + floor), col=3, lwd=3)
		}
		legend( "topleft", c( "Extrema Found", "Raw Peak Guess", "Clip Neighbors", 
			"Best Fit Curve", "Final Called Peak", "Local Median"),
			col=c(2,"brown",6,4,3,1), lwd=c(NA,2,2,2,3,1), lty=c(NA,2,1,1,1,3), 
			pch=c(19,NA,NA,NA,NA,NA), bg="white", cex=1.1)
		dev.flush()
		dev.print( png, plot.file, width=840, height=600, bg='white')
	}

	if ( ! validFit) return( NULL)

	out <- list( "center"=ans$center, "width"=ans$width,"height"=maxY, "x"=fitX, "y"=fitY,
			"volume"=volume, "residual"= residual, "floor"=floor, "type"=model)
	return(out)
}


`which.extrema` <- function( y, canonical.width=50) {

	cat( "  Selecting extrema..")
	nTrueM <- nTrueP <- rep( 0, times=length(y))

	# given the typical width of a peak (half width at half height), 
	# decide how far out from a potential extrema we need to test the slope of...
	nTests <- round( canonical.width * 0.33)
	if ( nTests < 3) nTests <- 3
	nPass <- round( nTests * 0.9)

	for ( n in 1:nTests) {
		dy <- diff( y, lag=n)
		dYm <- c( rep(0,n), dy)
		dYp <- c( dy, rep(0,n))
		nTrueM <- nTrueM + as.integer( dYm >= 0)
		nTrueP <- nTrueP + as.integer( dYp <= 0)
	}
	isTrueM <- which( nTrueM >= nPass)
	isTrueP <- which( nTrueP >= nPass)

	out <- intersect( isTrueM, isTrueP)
	cat( "  N =", length( out), "\n")
	return( out)
}


`cleanTailRegion` <- function( x, y, centerX, width, flank=c("left","right")) {

	# peaks on the tails higher than median can get cleaned
	flank <- match.arg( flank)
	N <- length(y)
	nonpeakX <- which( x < (centerX-width) | x > (centerX+width))
	#medY <- median( y[nonpeakX]) * 1.15
	medY <- quantile( y[nonpeakX], 0.25, na.rm=T) * 1.15
	whoLow <- which( y <= medY)
	whoHigh <- which( y > medY)

	if ( flank == "right") {
		firstVictim <- round( N * 0.66)
		victimsLow <- intersect( whoLow, which( x > centerX))
		if ( length( victimsLow) > 0) firstVictim <- victimsLow[1]
		victimsFlank <- intersect( whoHigh, firstVictim:N)
	} else {
		lastVictim <- round( N * 0.33)
		victimsLow <- intersect( whoLow, which( x < centerX))
		if ( length( victimsLow) > 0) lastVictim <- victimsLow[ length( victimsLow)]
		victimsFlank <- intersect( whoHigh, 1:lastVictim)
	}

	newYs <- rep( medY, times=length(victimsFlank))
	out <- y
	out[ victimsFlank] <- newYs
	return( out)
}


`scoreRIPpeaks` <- function( tbl, min.width=50, max.width=1000, cutoff=1, strand="Plus") {

	# given a table of found peaks, assign a score to each

	# times higher than the ambient background
	absHeight <- (tbl$height + tbl$floor)
	hTerm <- (absHeight - cutoff) / absHeight
	isNeg <- which( tbl$height <= 0 | absHeight <= 0 | hTerm <= 0)
	hTerm[ isNeg] <- 0.1

	# using a 2-sided test with a mesa top
	upperWid <- max.width * 1
	lowerWid <- min.width * 1
	absWidth <- rep.int( 1, nrow(tbl))
	absWidth <- ifelse( abs(tbl$width) > upperWid, (abs( tbl$width) / upperWid), absWidth)
	absWidth <- ifelse( abs(tbl$width) < lowerWid, (abs( tbl$width) / lowerWid), absWidth)
	# way too skinny...
	absWidth[ absWidth < 0.1] <- 0.1
	tooSkinny <- which( absWidth < 0.7)
	if ( length( tooSkinny) > 0) absWidth[ tooSkinny] <- 0.7 / absWidth[ tooSkinny]
	absWidth[ absWidth < 1] <- 1
	wTerm <- sqrt( 1 / absWidth)

	# drift too far
	absDrift <- abs(tbl$drift) / (max.width * 1)
	absDrift[ absDrift < 1] <- 1
	dTerm <- sqrt( 1 / absDrift)

	# relative H vs W:  there is a typical shape...remember 'width' is the HWHH
	relHW <- tbl$height / abs(tbl$width)
	maxRelative <- 30
	relHW[ relHW > maxRelative] <- maxRelative
	rTerm <- sqrt( relHW / maxRelative)

	return( data.frame( "loc"=tbl$center, "strand"=rep(strand,nrow(tbl)), "ht"=hTerm, 
			"wd"=wTerm, "dr"=dTerm, "rel"=rTerm, "score"=(hTerm*wTerm*dTerm*rTerm),
			stringsAsFactors=F))
}



`randomRIPpeaks` <- function( tbl, N=100, whoExtrema=NULL, cutoff.medians=2, min.width=50, max.width=1000, strand="Plus") {
				
	peak.type <- "gaussian"
	cat( "\nRandom control peaks..\n")

	# given a data object with X, Y...
	allX <- as.integer( tbl$x)
	allY <- as.numeric( tbl$y)
	medianY <- median( allY[ allY > 0])
	cutoffY <- medianY * cutoff.medians
	strand <- tbl$strand

	# the universe of possible peak tops is everyone higher than the cutoff
	if ( is.null( whoExtrema)) {
		whoExtrema <- which.extrema( allY, canonical.width=min.width)
	}
	N <- min( N, length( whoExtrema))

	# take a random sample of extrema
	randomWhos <- sample( whoExtrema, size=N*2)

	outM <- outSD <- outH <- outV <- outR <- outStart <- outStop <- outFloor <- vector()
	outType <- outStat <- outDrift <- vector()
	nout <- 0

	# include a neighbor hood about 4x bigger than the peak's footprint
	mid.width <- mean( c( min.width, max.width))
	extra <- as.integer( mid.width * 4)

	for ( who in randomWhos) {

		centerX <- allX[ who]
		centerY <- allY[ who]

		# the region to look at is a neighborhood around the center
		myPts <- max( 1, (who-extra)) : min( length(allX), (who+extra))
		ans <- findOneBestRIPpeak( allX[myPts], allY[myPts], centerX, centerY, cutoff=cutoffY, 
					min.width=min.width, max.width=max.width, peak.type=peak.type)
		if ( is.null( ans)) next

		nout <- nout + 1
		if ( nout > N) break

		outM[nout] <- myX <- ans$center
		outSD[nout] <- mySD <- ans$width
		outH[nout] <- maxY <- ans$height
		outFloor[nout] <- myFloor <- ans$floor
		outV[nout] <- ans$volume
		outR[nout] <- ans$residual
		outDrift[nout] <- myDrift <- centerX - myX
		fitX <- ans$x
		fitY <- ans$y
		outStart[nout] <- min( fitX)
		outStop[nout] <- max( fitX)
		outType[nout] <- ans$type

		outStat[nout] <- "good"
		if ( abs(mySD) > max.width*1) outStat[nout] <- "badWidth"
		if ( abs(myDrift) > max.width*1) outStat[nout] <- "badDrift"
		if ( myX < outStart[nout] | myX > outStop[nout]) outStat[nout] <- "badMean"
		if ( maxY < medianY/4) outStat[nout] <- "badHeight"
		if ( (maxY+myFloor) < medianY) outStat[nout] <- "badHeight"

		if ( nout %% 200 == 0) cat( "\n", strand, nout, myX, mySD, maxY, outStat[nout])
	}

	# do a bit of rounding...
	outM <- round( outM, digits=2)
	outSD <- round( outSD, digits=2)
	outH <- round( outH, digits=2)
	outFloor <- round( outFloor, digits=2)
	outStart <- round( outStart, digits=0)
	outStop <- round( outStop, digits=0)
	outV <- round( outV, digits=0)
	outR <- round( outR, digits=0)
	outDrift <- round( outDrift, digits=2)

	out <- data.frame( "center"=outM, "width"=outSD, "height"=outH, "floor"=outFloor, 
			"start"=outStart, "stop"=outStop, "volume"=outV, "residual"=outR,
			"drift"=outDrift, "type"=outType, "status"=outStat, stringsAsFactors=F)
	rownames(out) <- 1:nrow(out)

	ans <- scoreRIPpeaks( out, min.width=min.width, max.width=max.width, cutoff=medianY, strand=strand)
	out <- cbind( out, "score"=as.numeric(ans$score))
	
	return(out)
}

