# pipe.ChIPpeaks.R

# run a sample's  Alignment results thruough the ChIP peak picker

`pipe.ChIPpeaks` <-
function( sampleID, annotationFile="Annotation.txt", optionsFile="Options.txt", speciesID=NULL,
		results.path=NULL, loadWIG=FALSE, storage.mode="normal", 
		doPeakSearch=TRUE, peak.type="auto", cutoff.medians=3, canonical.width=50, 
		use.multicore=TRUE, p.value=0.25, verbose=FALSE, 
		watermark.gene=sub( "_[BS][0-9]+$", "", sampleID), 
		controlPeaksFile=NULL, scaledReadCount=NULL, visualize=interactive())
{

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nStarting pipe 'ChIPpeaks' on Sample:     ", sampleID, "\n")
	}
	startTotalT <- proc.time()
	grandTotalReads <- 0

	# get annotation and option facts we need...
	origSampleID <- originalSamplePairID( sampleID, annotationFile)
	annT <- readAnnotationTable( annotationFile)
	optT <- readOptionsTable( optionsFile)
	target <- getAndSetTarget( optionsFile, verbose=T)

	if ( use.multicore) {
		nCores <- as.integer( getOptionValue( optT, "nCores", notfound="4"))
		multicore.setup( nCores)
	}

	allSpecies <- getCurrentTargetSpecies()
	if ( ! is.null(speciesID)) allSpecies <- speciesID

	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	if ( ! file.exists( results.path)) dir.create( results.path, recursive=TRUE, showWarnings=FALSE)


	# do each species all the way through
	for( speciesID in allSpecies) {
	
	cat( "\n\nExtracting results for Species:  ", speciesID,"\n")
	startT <- proc.time()
	gc()
	setCurrentSpecies( speciesID)
	speciesPrefix <- getCurrentSpeciesFilePrefix()

	# find and/or load our WiggleBin data structure
	fileOutWIG <- file.path( results.path, "wig", paste( sampleID, speciesPrefix, "WIG.rda", sep="."))

	if ( !file.exists( fileOutWIG) || loadWIG) {
		
		cat( "\nBuilding Wiggles data structure from alignments...")
		pipe.AlignToWig( sampleID, annotationFile=annotationFile, optionsFile=optionsFile,
				results.path=results.path, dataType="ChIP-seq", mode=storage.mode)

		#get that data back when done... object 'wiggles'
		load( file=fileOutWIG)
	} else {
		cat( "\nLoading pre-existing Wiggles data...")
		load( file=fileOutWIG)
	}
	#wiggles <<- wiggles

	fileOutPeaks <- paste( sampleID, speciesPrefix, "ChIPpeaks.txt", sep=".")
	pathOut <- file.path( results.path, "ChIPpeaks", sampleID)
	if ( ! file.exists( pathOut)) dir.create( pathOut, recursive=TRUE, showWarnings=FALSE)
	fileOutPeaks <- file.path( pathOut, fileOutPeaks)

	peaks <- calcWigChIPpeaks( wiggles, fileout=fileOutPeaks, doPeakSearch=doPeakSearch,
				peak.type=peak.type, 
				cutoff.medians=cutoff.medians, canonical.width=canonical.width,
				use.multicore=use.multicore, verbose=verbose, 
				controlPeaksFile=controlPeaksFile, 
				scaledReadCount=scaledReadCount, watermark.gene=watermark.gene,
				visualize=visualize)

	reads <- WIG_getTotalRawReads( wiggles)
	grandTotalReads <- grandTotalReads + reads$Unique + reads$Multi

	myTime <- elapsedProcTime( startT, proc.time(), N=(reads$Unique + reads$Multi))
	if (verbose) {
		cat( "\n\nFinished Species:  ", speciesID)
		cat( "\n\nSpecies Timing Stats: \n")
		print( myTime)
		gc()
	}

	}  # end of all species...

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\n\nFinished pipe.ChIPpeaks:  ", sampleID, "\tSpecies set:  ", allSpecies)
	}
	if ( length( allSpecies) > 1) {
		myTime <- elapsedProcTime( startTotalT, proc.time(), N=grandTotalReads)
		cat( "\n\nSample Timing Stats: \n")
		print( myTime)
	}


	# do all the post-pick operations
	ChIPpeaksAddDNA( sampleID, optionsFile=optionsFile, p.value=p.value, results.path=results.path)
	ChIPpeaksToExcel( sampleID, optionsFile=optionsFile, p.value=p.value, results.path=results.path)
	ChIPpeaksToHTML( sampleID, optionsFile=optionsFile, tailWidth=300, p.value=p.value, results.path=results.path,
				scaledReadCount=scaledReadCount, visualize=visualize)

	return( sampleID)
}


`dispatch.ChIPpeaksPipeline` <- function( sampleID, annotationFile="Annotation.txt",
					optionsFile="Options.txt", banner="", mode="normal",
					maxReads=NULL, pause=0, results.path=NULL) {

	commandLine <- paste( "X11( width=8, height=6, xpos=20, ypos=20, bg='white'); ",
				" pipe.ChIPpeaks( sampleID=\"", 
				sampleID, "\", annotationFile=\"", annotationFile, 
				"\", optionsFile=\"", optionsFile, 
				"\", banner=\"", banner, 
				"\", mode=\"", mode, "\"", 
				", maxReads=", if (is.null(maxReads)) "NULL" else as.integer(maxReads), 
				", pause=", as.integer(pause), 
				", results.path=", if (is.null(results.path)) "NULL" else 
					paste("\"",results.path,"\"",sep=""), 
				" )", sep="")

	spawnRcmd( commandLine, Rpackage="DuffyRNAseq", logFile=paste( sampleID, "ChIPpeaks.log.txt", sep="."))

	return()
}


`ChIPpeaksToHTML` <- function( sampleID, optionsFile="Options.txt", p.value=0.25, results.path=NULL,
				scaledReadCount=NULL, visualize=interactive(), ...) {

	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.null(results.path)) results.path <- getOptionValue( optionsFile, "results.path", notfound=".")
	wigfile <- file.path( results.path, "wig", paste( sampleID, prefix, "WIG.rda", sep="."))
	load( wigfile)
	chipfile <- file.path( results.path, "ChIPpeaks", sampleID, 
				paste( sampleID, prefix, "ChIPpeaks.txt", sep="."))
	tbl <- read.delim( chipfile, as.is=T)

	htmlfile <- file.path( results.path, "ChIPpeaks", sampleID, 
				paste( sampleID, prefix, "ChIPpeaks.html", sep="."))
	localPngPath <- "ChIP.plots"
	pngPath <- file.path( results.path, "ChIPpeaks", sampleID, localPngPath)
	if ( ! file.exists( pngPath)) dir.create( pngPath, recursive=T)


	myColumns <- c( "GENE_ID", "VPM", "P_Value", "Score", "Cheight", "Ccenter", "Ctimes", 
			"Fpvalue", "Rpvalue", "Cpvalue", "Fscore", "Rscore", "Cscore", "Ftype", "Rtype", "Ctype",
			"Fcenter", "Rcenter", "Fheight", "Rheight", "Fwidth", "Rwidth", "Cwidth",
			"Ffloor", "Rfloor", "Fstatus", "Rstatus", "Cstatus", "Model_Error", "Fstart", 
			"Rstart", "Fstop","Rstop", "Fvpm", "Rvpm")
	html <- tbl[ , myColumns]

	# trim to only plot the better P value peaks
	keep <- which( html$P_Value <= p.value)
	if ( length(keep) < 100) keep <- 1 : min(100,nrow(tbl))
	html <- html[ keep, ]
	rownames(html) <- 1:nrow(html)
	tbl <- tbl[ keep, ]

	html$GENE_ID <- paste( html$GENE_ID, as.integer(round(html$Ccenter)), sep="_")

	# any formatting...
	for( j in c( 2,5,6,17:20,24:25,30:33)) {
		html[ ,j] <- formatC( html[ ,j], format="d", digits=0)
	}
	for( j in c( 3,8:10)) {
		html[ ,j] <- formatC( html[ ,j], format="e", digits=1)
	}
	for( j in c( 4, 11:13)) {
		html[ ,j] <- formatC( html[ ,j], format="f", digits=3)
	}
	for( j in c( 7,21:23,29,34:35)) {
		html[ ,j] <- formatC( html[ ,j], format="f", digits=2)
	}

	# make a high level plot
	cat( "\nChromosome Plots...")
	seqIDset <- rownames( wiggles$Info$RawReads)
	extraText <- vector()
	extraText[1] <- paste( "<BR>",  "Overview Plots: &nbsp; ")
	if ( visualize) {
	    checkX11( width=10, height=7, bg="white")
	    for ( seqID in seqIDset) {
		plotWIGregion( wiggles, seqID, showDetect=F, label=sampleID)
		genomefile <- paste( "Overview", seqID, "png", sep=".")
		genomepath <- file.path( pngPath, genomefile)
		genomepathlocal <- file.path( localPngPath, genomefile)
		dev.print( png, genomepath, width=900, height=600, bg='white')
		extraText <- c( extraText, paste( " &nbsp;", as.link( genomepathlocal, seqID)))
		cat( "  ", seqID)
	    }

	    # the TF gene too...
	    geneID <- sub( "_.+", "", sampleID)
	    if ( geneID %in% getCurrentGeneMap()$GENE_ID) {
	    	plotWIGgene( wiggles, geneID, showDetect=T, label=sampleID, tailWidth=1200, minYmax=40)
		genefile <- paste( "Overview", geneID, "png", sep=".")
		genepath <- file.path( pngPath, genefile)
		genepathlocal <- file.path( localPngPath, genefile)
		dev.print( png, genepath, width=900, height=600, bg='white')
		extraText <- c( extraText,  " <BR>   ")
		extraText <- c( extraText, paste( " Transcription Factor 'Watermark' Gene Plot: &nbsp;", 
				as.link( genepathlocal, geneID)))
		cat( "  ", geneID)
	    }
	}

	table2html( html, htmlfile, title=paste( "Top ChIP Peaks:  &nbsp;  &nbsp;  ", sampleID), 
		extraHTMLtext=extraText, linkPaths=localPngPath)

	cat( "\n")
	if ( visualize) {
	    for ( i in 1:nrow(tbl)) {
		thisG <- html$GENE_ID[i]
		pngfile <- file.path( pngPath, paste( thisG, "png", sep="."))
		png( pngfile, width=900, height=600)
		plotChIPpeak( wiggles, tbl, chipRow=i, label=sampleID, scaledReadCount=scaledReadCount,
				...)
		dev.off()
		cat( "\r", i, thisG)
	    }
	} else {
		cat( "\nNo Plotting in batch mode...")
	}
	cat( "\n")
	return()
}


`ChIPpeaksToExcel` <- function( sampleID, optionsFile="Options.txt", p.value=0.25, results.path=NULL) {

	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.null(results.path)) results.path <- getOptionValue( optionsFile, "results.path", notfound=".")
	chipfile <- file.path( results.path, "ChIPpeaks", sampleID, 
				paste( sampleID, prefix, "ChIPpeaks.txt", sep="."))
	tbl <- read.delim( chipfile, as.is=T)

	csvfile <- sub( "txt$", "Excel.csv", chipfile)

	# trim to only keep the better P value peaks
	keep <- which( tbl$P_Value <= p.value)
	if ( length(keep) < 100) keep <- 1 : min(100,nrow(tbl))
	tbl <- tbl[ keep, ]

	excelterm <- paste( tbl$GENE_ID, round( tbl$Ccenter), sep="_")
	
	tbl2 <- data.frame( "GENE_LINK"=excelterm, tbl, stringsAsFactors=F)
	out <- addExcelLinksToTable( tbl2, linkColumnNames="GENE_LINK", linkPaths="ChIP.plots")

	# turn the DNA into fasta format...
	if ( "DNA" %in% colnames(out)) {
		todo <- which( out$DNA != "")
		newdna <- paste( ">", out$GENE_ID[todo], "\n", out$DNA[todo], sep="")
		out$DNA[todo] <- newdna
	}

	write.table( out, csvfile, sep=",", quote=T, qmethod="double", row.names=F)
	return( out)
}


`ChIPpeaksAddDNA` <- function( sampleID, optionsFile="Options.txt", p.value=0.25, results.path=NULL) {

	prefix <- getCurrentSpeciesFilePrefix()
	gmap <- getCurrentGeneMap()
	if ( is.null(results.path)) results.path <- getOptionValue( optionsFile, "results.path", notfound=".")
	chipfile <- file.path( results.path, "ChIPpeaks", sampleID, 
				paste( sampleID, prefix, "ChIPpeaks.txt", sep="."))
	tbl <- read.delim( chipfile, as.is=T)
	gene2seqPtr <- match( tbl$GENE_ID, gmap$GENE_ID, nomatch=0)
	seqID <- rep( NA, times=nrow(tbl))
	seqID[ gene2seqPtr > 0] <- gmap$SEQ_ID[ gene2seqPtr]

	myDNA <- rep( "", times=nrow(tbl))
	curSeqID <- ""
	genomicFastaFile <- getOptionValue( optionsFile, genomicFastaFile)
	if ( is.na( genomicFastaFile)) stop( "no genomic Fasta filename in 'options' file")

	# only visit rows that are decent
	goodPvalues <- which( tbl$P_Value <= p.value)
	for( i in goodPvalues) {
		thisSeqID <- seqID[i]
		if ( is.na( thisSeqID)) next
		if ( thisSeqID != curSeqID) {
			curSeqID <- thisSeqID
			dna <- getFastaSeqFromFilePath( genomicFastaFile, thisSeqID)
		}

		myfrom <- tbl$Fcenter[i]
		myto <- tbl$Rcenter[i]
		n <- myto - myfrom + 1
		if ( is.na(n)) next
		cntr <- as.integer( (myto + myfrom) / 2)
		if ( n < 100) {
			myfrom <- cntr - 50
			myto <- cntr + 50
		}
		if ( n > 200) {
			myfrom <- cntr - 100
			myto <- cntr + 100
		}
		myDNA[i] <- substr( dna, myfrom, myto)
	}
	tbl$DNA <- myDNA
	write.table( tbl, chipfile, sep="\t", quote=F, row.names=F)

	ndone <- sum( nchar( myDNA) > 0)
	cat( "\nAppended DNA fragments for ", ndone, " ChIP peaks.")
	return( NULL)
}


`plotMultiChIPpeaks` <- function( sampleIDset=NULL, annotationFile="Annotation.txt", results.path=NULL,
				optionsFile="Options.txt", speciesID=NULL, Ngenes=100, 
				label="", ...) {

	annT <- readAnnotationTable( annotationFile)
	if ( is.null( sampleIDset)) {
		sampleIDset <- annT$SampleID
	} else {
		annT <- subset( annT, SampleID %in% sampleIDset)
		sampleIDset <- annT$SampleID
	}

	optT <- readOptionsTable( optionsFile)
	if ( is.null(results.path)) results.path <- getOptionValue( optT, "results.path", notfound=".")
	target <- getAndSetTarget( optionsFile, verbose=T)
	if ( ! is.null( speciesID)) setCurrentSpecies( speciesID)
	speciesID <- getCurrentSpecies()
	prefix <- getCurrentSpeciesFilePrefix()

	# gather up all those wiggle tracks
	WIGlist <- loadMultipleWIGs( sampleIDset, speciesID, results.path)

	# gather up all those ChIP peak files
	prefix <- getCurrentSpeciesFilePrefix()
	PEAKlist <- vector( mode="list")
	allGenes <- vector()
	for ( i in 1:length(sampleIDset)) {
		sampleID <- sampleIDset[i]
		chipfile <- file.path( results.path, "ChIPpeaks", sampleID, 
				paste( sampleID, prefix, "ChIPpeaks.txt", sep="."))
		if ( file.exists( chipfile)) {
			tbl <- read.delim( chipfile, as.is=T)
			n <- min( Ngenes, nrow(tbl))
			PEAKlist[[i]] <- tbl[ 1:n, ]
			allGenes <- c( allGenes, tbl$GENE_ID[1:n])
		} else {
			cat( "\nChIP peak file not found: ", basename(chipfile), "   Skipped..")
		}
	}

	allGenes <- sort( unique( allGenes))
	outpath <- file.path( results.path, "ChIPpeaks", "Multi_Sample.plots")
	if( ! file.exists( outpath)) dir.create( outpath, recursive=T)

	cat( "\nPlotting Genes: ", length( allGenes), "\n")
	for ( ig in 1:length(allGenes)) {
		g <- allGenes[ig]
		thisfile <- paste( g, "png", sep=".")
		thisfile <- file.cleanSpecialCharactersFromFileName( thisfile)
		thisfile <- file.path( outpath, thisfile)

		plotMultiWIGgene( WIGlist, annT$Color, gene=g, ...)
		dev.print( png, thisfile, width=900, height=600, bg='white')
		cat( "\r", ig, g, "   ")
	}
	cat( "\nDone.\n")

	return( length(allGenes))
}


`ChIPpeakSeparationSummary` <- function( sampleIDset, optionsFile="Options.txt", p.value=0.25, results.path=NULL) {

	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.null(results.path)) results.path <- getOptionValue( optionsFile, "results.path", notfound=".")

	peakSep <- vector()
	for ( sampleID in sampleIDset) {
		chipfile <- file.path( results.path, "ChIPpeaks", sampleID, 
				paste( sampleID, prefix, "ChIPpeaks.txt", sep="."))
		tbl <- read.delim( chipfile, as.is=T)
		keep <- which( tbl$P_Value <= p.value)
		tbl <- tbl[ keep, ]
		cat( "\nAvg Peak Seperation: ", sampleID, nrow(tbl))
		if ( nrow(tbl) < 1) next

		thisSep <- tbl$Rcenter - tbl$Fcenter
		myMean <- mean( thisSep, na.rm=T)
		mySD <- sd( thisSep, na.rm=T)
		cat( "\tMean: ", myMean, "\tSD: ", mySD)
		peakSep <- c( peakSep, thisSep)
	}
	myMean <- mean( peakSep, na.rm=T)
	mySD <- sd( peakSep, na.rm=T)

	out <- list( "mean"=myMean, "sd"=mySD)
	return( out)
}

