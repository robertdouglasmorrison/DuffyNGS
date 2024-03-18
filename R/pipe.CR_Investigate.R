# pipe.CR_Investigate.R


`dispatch.CR_Investigate` <- function( sampleID, annotationFile="Annotation.txt",
				optionsFile="Options.txt", mode="normal",
				blastIndex=getOptionValue( optionsFile, "blastIndex", notfound="nt"),
				doCR=TRUE, doBlast=TRUE, maxReads=500000, maxTime=1000, maxCycles=10, 
				ratePerCycle=1, maxCR=4000, pause=0,
				nIterations=1000, nBest=10, results.path=NULL, makePlots=TRUE, verbose=TRUE) {

	commandLine <- paste( "checkX11( width=8, height=6, xpos=30, ypos=30, bg='white'); ",
				" pipe.CR_Investigate( sampleID=\"", 
				sampleID, "\", annotationFile=\"", annotationFile, 
				"\", optionsFile=\"", optionsFile, 
				"\", blastIndex=\"", blastIndex,
				"\", mode=\"", mode, "\", doCR=", doCR, ", doBlast=", doBlast, 
				", maxReads=", as.integer(maxReads), ", maxTime=", as.integer(maxTime),
				", maxCycles=", as.integer(maxCycles), 
				", ratePerCycle=", as.numeric(ratePerCycle), 
				", maxCR=", as.integer(maxCR), 
				", pause=", as.integer(pause), 
				", nIterations=", as.integer(nIterations), ", nBest=", as.integer(nBest),
				", results.path=", if (is.null(results.path)) "NULL" else 
				paste( "\"", results.path, "\"", sep=""),
				", makePlots=", makePlots,
				", verbose=", verbose,
				" )", sep="")

	spawnRcmd( commandLine, Rpackage="DuffyNGS", logFile=paste( sampleID, 
			"CR_Investigate.log.txt", sep="."))

	return()
}



`pipe.CR_Investigate` <- function( sampleID, annotationFile="Annotation.txt",
		optionsFile="Options.txt", mode=c( "normal", "QuickQC", "genomic"),
		blastIndex=getOptionValue( optionsFile, "blastIndex", notfound="nt"),
		doCR=FALSE, doBlast=FALSE, maxReads=500000, maxTime=1000, maxCycles=10, 
		ratePerCycle=NULL, maxCR=4000, pause=0, 
		nIterations=1000, nBest=10, results.path=NULL, makePlots=TRUE,
		trim5=NULL, trim3=NULL, verbose=TRUE) {

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\n\nStarting 'Consensus Reads Investigate' for Sample:     ", sampleID, "\n\n")
	}
	gc()
	
	require( Biostrings)

	startT <- proc.time()

	mode <- match.arg( mode)

	# get the file of noHit reads
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=verbose)
	}
	infile <- paste( sampleID, "noHits", "fastq", sep=".") 
	infile <- file.path( results.path, "fastq", infile)
	if ( mode == "QuickQC") {
		infile <- file.path( results.path, "fastq", paste( sampleID, "QuickQC.noHits.fastq", sep="."))
	}
	if ( mode == "genomic") {
		infile <- file.path( results.path, "fastq", paste( sampleID, "not.genomic.fastq", sep="."))
	}
	nNoHits <- getFileLineCount( infile, sampleID) / 4

	# either reload or build the Unique Short Reads data structure
	usrFile <- file.path( results.path, "USR", paste( sampleID, "USR.rda", sep="."))
	if ( file.exists( usrFile)) {
		cat( "\nLoading existing USR dataset:  ", usrFile)
		load( usrFile, envir=.GlobalEnv)
	} else {
		if ( is.null( trim5)) trim5 <- as.numeric( getOptionValue( optionsFile, "trim5", notfound=0))
		if ( is.null( trim3)) trim3 <- as.numeric( getOptionValue( optionsFile, "trim3", notfound=0))
		ans <- USR_setup( infile, sampleID, results.path, trim5=trim5, trim3=trim3, 
				Nkeep=maxReads)
		doCR <- doBlast <- TRUE
	}

	# set up folder to hold results
	cr.path <- file.path( results.path, "CR", sampleID)
	if ( ! file.exists( cr.path)) dir.create( cr.path, recursive=TRUE, showWarnings=FALSE)
	crPngPath <- file.path( cr.path, paste( "CR", "pngPlots", sep="."))
	if ( ! file.exists( crPngPath)) dir.create( crPngPath, showWarnings=FALSE)

	# run the automatic CR growing tool
	crFile <- file.path( cr.path, paste( sampleID, "CR.rda", sep="."))
	if ( doCR || !file.exists(crFile)) {
		bestCR <- autoRunCR( nBest=nBest, nIterations=nIterations, maxTime=maxTime, 
					maxCycles=maxCycles, ratePerCycle=ratePerCycle, 
					maxCR=maxCR, pause=pause, makePlots=makePlots,
					contextFile=crFile, pngPath=crPngPath,
					label=sampleID)
		CRT_best <<- bestCR
		saveCRTcontext( crFile)
		doBlast <- TRUE
	} else {
		cat( "\nLoading existing CR dataset:  ", crFile)
		load( crFile, envir=.GlobalEnv)
		bestCR <- CRT_best
	}

	# file to hold the Blast results
	xmlOutFile <- file.path( cr.path, paste( sampleID, "blastOutput.xml", sep="."))

	# see what these maight be via Blast...
	if ( !file.exists(xmlOutFile)) doBlast <- TRUE
	blastResults <<- NULL
	if ( ! is.null( blastIndex)) {
		blastResults <<- CRblaster( sampleID, crIDs=NULL, nBest=nBest, xmlOutFile=xmlOutFile,
					doBlast=doBlast, optionsFile=optionsFile, blastIndex=blastIndex, 
					evalue=0.1, wordsize=11)
	}

	# turn the Blast results into a table form text file
	if ( ! is.null( blastResults)) {

		fileout <- file.path( results.path, "CR", sampleID, paste( sampleID, "CR.BlastSummary.txt", sep="."))
		tbl <- summarizeCRblastOutput( blastResults, fileout=fileout)

		# and into clickable HTML, with plots for each
		cat( "\nMaking HTLM table...")
		htmlFile <- sub( ".txt$", ".html", fileout)
		htmlTitle <- paste( "Consensus Reads:   'no hits' reads from sample:  ", sampleID)
		htmlPngPath <- paste( "CR", "pngPlots", sep=".")
		# try to make the consensus sequence be line wrapped and 'cut/paste'-able
		tbl$CONSENSUS_SEQUENCE <- as.text.Fasta( as.Fasta( tbl$CR_ID, tbl$CONSENSUS_SEQUENCE),
				line.width=50)
		tbl$CONSENSUS_SEQUENCE <- sub( "\n", " <br> ", tbl$CONSENSUS_SEQUENCE)
		tbl$TEXT_DESCRIPTION_OF_BEST_HITS <- gsub( " | ", " <br> ", tbl$TEXT_DESCRIPTION_OF_BEST_HITS, 
				fixed=TRUE)
		table2html( tbl, htmlFile, title=htmlTitle, linkColumnNames="CR_ID", linkPaths=htmlPngPath)

		crToPlot <- as.integer( sub( "^cr_", "", tbl$CR_ID))
		if (makePlots) {
		    for( i in crToPlot) {
			pngFile <- file.path( crPngPath, paste( "cr_",i,".png",sep=""))
			png( pngFile, width=1000, height=700, bg="white")
			plotOneCRT( i, label=sampleID)
			dev.off()
		    }
		}

		# when we did the blast, also make a small summary image
		pipe.CR_Summary( sampleID, optionsFile=optionsFile, results.path=results.path, max.show=20)
	}
	
	# done...
	USR_cleanup()
	CR_cleanup()

	gc()
	return()
}


`pipe.CR_Summary` <- function( sampleID, optionsFile="Options.txt", results.path=NULL, max.show=20, 
				nTermsBlastHit=4, useLog=NULL) {

	require( Biostrings)

	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	}
	cr.path <- file.path( results.path, "CR", sampleID)
	filein <- file.path( cr.path, paste( sampleID, "CR.BlastSummary.txt", sep="."))
	tbl <- read.delim( filein, as.is=T)
	N <- nrow( tbl)

	# see if we can deduce the typical junk
	dnaStrs <- tbl$CONSENSUS_SEQUENCE
	explainGroup <- rep.int( "", N)
	polySize <- 50
	for ( base in c("A","C","G","T","N")) {
		thisName <- paste( "Poly", base, sep="_")
		polyBase <- paste( rep.int( base, polySize), collapse="")
		simScores <- pairwiseAlignment( dnaStrs, polyBase, type="local-global", scoreOnly=T)
		goodHits <- which( simScores >= polySize)
		if ( length( goodHits)) explainGroup[ goodHits] <- thisName
	}

	fwdAdapter <- getOptionValue( optionsFile, "ForwardAdapter", notfound="", verbose=F)
	revAdapter <- getOptionValue( optionsFile, "ReverseAdapter", notfound="", verbose=F)
	if ( nchar(fwdAdapter) > 10) {
		simScores <- pairwiseAlignment( dnaStrs, fwdAdapter, type="local-global", scoreOnly=T)
		goodHits <- which( simScores >= nchar(fwdAdapter))
		if ( length( goodHits)) explainGroup[ goodHits] <- "ForwardAdapter"
		simScores <- pairwiseAlignment( dnaStrs, myReverseComplement(fwdAdapter), type="local-global", scoreOnly=T)
		goodHits <- which( simScores >= nchar(fwdAdapter))
		if ( length( goodHits)) explainGroup[ goodHits] <- "RevComp(FwdAdapter)"
	}
	if ( nchar(revAdapter) > 10) {
		simScores <- pairwiseAlignment( dnaStrs, revAdapter, type="local-global", scoreOnly=T)
		goodHits <- which( simScores >= nchar(revAdapter))
		if ( length( goodHits)) explainGroup[ goodHits] <- "ReverseAdapter"
		simScores <- pairwiseAlignment( dnaStrs, myReverseComplement(revAdapter), type="local-global", scoreOnly=T)
		goodHits <- which( simScores >= nchar(revAdapter))
		if ( length( goodHits)) explainGroup[ goodHits] <- "RevComp(RevAdapter)"
	}
	
	# for all the chunks not yet explained, grab the first K words from their free text
	freeText <- convertHypertext( tbl$BLAST_HIT)
	freeTerms <- strsplit( freeText, split=" +")
	freeText <- sapply( freeTerms, function(x) {
			n <- min( length(x), nTermsBlastHit); 
			return( paste( x[1:n], collapse=" "))
		})
	tooLong <- which( nchar( freeText) > 50)
	if ( length(tooLong)) freeText[ tooLong] <- paste( substr( freeText[tooLong], 1, 47), "...", sep="")
	#freeText <- sub( "(^[A-Za-z0-9.:]+ [A-Za-z0-9.:]+)( .+$)", "\\1", tbl$BLAST_HIT)
	notExplained <- which( explainGroup == "")
	explainGroup[notExplained] <- freeText[ notExplained]
	explainGroup[ explainGroup == ""] <- "No BLAST hits found"

	myPcts <- tbl$PCT_READS
	myCounts <- tbl$N_READS
	ans <- tapply( myPcts, factor(explainGroup), FUN=sum, na.rm=T)
	ans2 <- tapply( myCounts, factor(explainGroup), FUN=sum, na.rm=T)

	ord <- order( ans, decreasing=T)
	ans <- ans[ ord]
	ans2 <- ans2[ ord]

	# very small fractions can break log scale...
	pcts <- as.numeric(ans)
	drops <- which( pcts < 0.0005)
	if ( length(drops)) {
		ans <- ans[ -drops]
		pcts <- pcts[ -drops]
		ans2 <- ans2[ -drops]
	}

	# control how many get drawn, and how...
	Nshow <- length(ans)
	if ( Nshow > max.show) {
		ans <- ans[1:max.show]
		pcts <- pcts[1:max.show]
		Nshow <- max.show
	}
	if (is.null( useLog)) useLog <- ( max(pcts)/min(pcts) > 20)
	if ( useLog) {
		xLimits <- range( ans, na.rm=T) * c( 0.75,5)
		myLog="x"
		xLabel <- "Percent of No Hit Reads  (log scale)" 
		xMidpt <- logmean( xLimits)
	} else {
		xLimits <- c( 0, max(ans,na.rm=T)*2) 
		myLog=""
		xLabel <- "Percent of No Hit Reads" 
		xMidpt <- mean( xLimits)
	}

	# make that plot
	savMAI <- par( "mai")
	par( mai=c(1,1,0.8,0.4))
	mp <- barplot( pcts, main=paste( "NoHit 'Consensus Reads' Summary:    ", sampleID), 
			horiz=T, xlab=xLabel, xlim=xLimits, log=myLog, col=rainbow( Nshow, end=0.75))

	textFont <- 2
	textCex <- 1
	textShow <- paste( names(ans), "  =", ifelse( pcts > 0.05, as.percent(pcts,big.value=100), 
				as.percent(pcts,big.value=100,digits=4)), sep="")
	myX <- pcts
	myPos <- ifelse( pcts >= xMidpt,  2, 4)
	if ( Nshow > 10) textCex <- 1 - (Nshow-10)*0.02
	text( myX, mp, textShow, pos=myPos, font=textFont, cex=textCex)

	dev.flush()
	printPlot( file.path( cr.path, paste( sampleID, "CR.SummaryPlot.pdf", sep="."))
	par( mai=savMAI)

	# make a text version too
	out <- data.frame( "Blast.Explanation"=names(ans), "NoHit.Percentage"=pcts, "N_READS"=as.numeric(ans2), 
			stringsAsFactors=F)
	outFile <- file.path( cr.path, paste( sampleID, "CR.SummaryPercentages.txt", sep="."))
	write.table( out, outFile, sep="\t", quote=F, row.names=F)

	return( ans)
}


`pipe.Extract.CR_Summary` <- function( sampleIDset, annotationFile="Annotation.txt", optionsFile="Options.txt", 
				results.path=NULL, verbose=TRUE) {

	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=verbose)
	}

	out <- data.frame()
	for ( sampleID in sampleIDset) {

		cr.path <- file.path( results.path, "CR", sampleID)
		inFile <- file.path( cr.path, paste( sampleID, "CR.SummaryPercentages.txt", sep="."))
		if ( ! file.exists(inFile)) {
			cat( "\nWarning: CR Summary file not found: ", inFile)
			next
		}
		tbl <- read.delim( inFile, as.is=T)
		sml <- data.frame( "SampleID"=sampleID, tbl, stringsAsFactors=F)
		out <- rbind( out, sml)
	}

	# final order by No Hit peroentage
	ord <- order( out$NoHit.Percentage, decreasing=T)
	out <- out[ ord, ]
	if ( nrow(out)) rownames(out) <- 1:nrow(out)

	return( out)
}

