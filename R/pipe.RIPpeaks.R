# pipe.RIPpeaks.R

# run a sample's  Alignment results thruough the RIP peak picker

`pipe.RIPpeaks` <-
function( sampleID, annotationFile="Annotation.txt", optionsFile="Options.txt", speciesID=NULL,
		results.path=NULL, loadWIG=FALSE, storage.mode="normal", 
		doPeakSearch=TRUE, peak.type="auto", cutoff.medians=3, min.width=50, max.width=1000,
		use.multicore=FALSE, p.value=0.075, verbose=FALSE, 
		controlPeaksFile=NULL, scaledReadCount=NULL, visualize=interactive())
{

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nStarting pipe 'RIPpeaks' on Sample:     ", sampleID, "\n")
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
				results.path=results.path, dataType="RIP-seq", mode=storage.mode)

		#get that data back when done... object 'wiggles'
		load( file=fileOutWIG)
	} else {
		cat( "\nLoading pre-existing Wiggles data...")
		load( file=fileOutWIG)
	}
	#wiggles <<- wiggles

	fileOutPeaks <- paste( sampleID, speciesPrefix, "RIPpeaks.txt", sep=".")
	pathOut <- file.path( results.path, "RIPpeaks", sampleID)
	if ( ! file.exists( pathOut)) dir.create( pathOut, recursive=TRUE, showWarnings=FALSE)
	fileOutPeaks <- file.path( pathOut, fileOutPeaks)

	peaks <- calcWigRIPpeaks( wiggles, fileout=fileOutPeaks, doPeakSearch=doPeakSearch,
				peak.type=peak.type, cutoff.medians=cutoff.medians, 
				min.width=min.width, max.width=max.width,
				use.multicore=use.multicore, verbose=verbose, 
				controlPeaksFile=controlPeaksFile, scaledReadCount=scaledReadCount, 
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
		cat( "\n\nFinished pipe.RIPpeaks:  ", sampleID, "\tSpecies set:  ", allSpecies)
	}
	if ( length( allSpecies) > 1) {
		myTime <- elapsedProcTime( startTotalT, proc.time(), N=grandTotalReads)
		cat( "\n\nSample Timing Stats: \n")
		print( myTime)
	}


	# do all the post-pick operations
	RIPpeaksAddDNA( sampleID, optionsFile=optionsFile, results.path=results.path, p.value=p.value)
	RIPpeaksToExcel( sampleID, optionsFile=optionsFile, results.path=results.path, p.value=p.value)
	RIPpeaksToHTML( sampleID, optionsFile=optionsFile, results.path=results.path, tailWidth=1000, p.value=p.value,
				scaledReadCount=scaledReadCount, visualize=visualize)

	return( sampleID)
}


`dispatch.RIPpeaksPipeline` <- function( sampleID, annotationFile="Annotation.txt",
					optionsFile="Options.txt", banner="", mode="normal",
					maxReads=NULL, pause=0, results.path=NULL) {

	commandLine <- paste( "X11( width=8, height=6, xpos=20, ypos=20, bg='white'); ",
				" pipe.RIPpeaks( sampleID=\"", 
				sampleID, "\", annotationFile=\"", annotationFile, 
				"\", optionsFile=\"", optionsFile, 
				"\", banner=\"", banner, 
				"\", mode=\"", mode, "\"", 
				", maxReads=", if (is.null(maxReads)) "NULL" else as.integer(maxReads), 
				", pause=", as.integer(pause), 
				", results.path=", if (is.null(results.path)) "NULL" else 
					paste("\"",results.path,"\"",sep=""), 
				" )", sep="")

	spawnRcmd( commandLine, Rpackage="DuffyRNAseq", logFile=paste( sampleID, "RIPpeaks.log.txt", sep="."))

	return()
}


`RIPpeaksToHTML` <- function( sampleID, optionsFile="Options.txt", results.path=NULL, p.value=0.05,
				scaledReadCount=NULL, visualize=interactive(), ...) {

	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.null( results.path)) results.path <- getOptionValue( optionsFile, "results.path", notfound=".")
	wigfile <- file.path( results.path, "wig", paste( sampleID, prefix, "WIG.rda", sep="."))
	load( wigfile)
	ripfile <- file.path( results.path, "RIPpeaks", sampleID, 
				paste( sampleID, prefix, "RIPpeaks.txt", sep="."))
	tbl <- read.delim( ripfile, as.is=T)

	htmlfile <- file.path( results.path, "RIPpeaks", sampleID, 
				paste( sampleID, prefix, "RIPpeaks.html", sep="."))
	localPngPath <- "RIP.plots"
	pngPath <- file.path( results.path, "RIPpeaks", sampleID, localPngPath)
	if ( ! file.exists( pngPath)) dir.create( pngPath, recursive=T)


	myColumns <- c( "GENE_ID", "VPM", "P_Value", "Score", "Height", "Center", 
			"Type", "Strand", "Start", "Stop", "Width", "Times", "Floor",
			"Model_Error", "status")
	html <- tbl[ , myColumns]

	# trim to only plot the better P value peaks
	keep <- which( html$P_Value <= p.value)
	if ( length(keep) < 100) keep <- 1 : min(100,nrow(tbl))
	html <- html[ keep, ]
	rownames(html) <- 1:nrow(html)
	tbl <- tbl[ keep, ]

	html$GENE_ID <- paste( html$GENE_ID, as.integer(round(html$Center)), sep="_")

	# any formatting...
	for( j in c( 2,5,6,9:13)) {
		html[ ,j] <- formatC( html[ ,j], format="d", digits=0)
	}
	for( j in c( 3)) {
		html[ ,j] <- formatC( html[ ,j], format="e", digits=1)
	}
	for( j in c( 4)) {
		html[ ,j] <- formatC( html[ ,j], format="f", digits=3)
	}
	#for( j in c( 7,21:23,29,34:35)) {
	#	html[ ,j] <- formatC( html[ ,j], format="f", digits=2)
	#}

	# make a high level plot
	cat( "\nChromosome Plots...")
	seqIDset <- rownames( wiggles$Info$RawReads)
	extraText <- vector()
	extraText[1] <- paste( "<BR>",  "Overview Plots: &nbsp; ")
	if ( visualize) {
	    checkX11( width=10, height=7, bg="white")
	    for ( seqID in seqIDset) {
		plotWIGregion( wiggles, seqID, showDetect=F, label=sampleID)
		dev.flush()
		genomefile <- paste( "Overview", seqID, "png", sep=".")
		genomepath <- file.path( pngPath, genomefile)
		genomepathlocal <- file.path( localPngPath, genomefile)
		dev.print( png, genomepath, width=900, height=600, bg='white')
		extraText <- c( extraText, paste( " &nbsp;", as.link( genomepathlocal, seqID)))
		cat( "  ", seqID)
	    }
	}

	table2html( html, htmlfile, title=paste( "Top RIP Peaks:  &nbsp;  &nbsp;  ", sampleID), 
		extraHTMLtext=extraText, linkPaths=localPngPath)

	cat( "\n")
	if ( visualize) {
	    for ( i in 1:nrow(tbl)) {
		thisG <- html$GENE_ID[i]
		pngfile <- file.path( pngPath, paste( thisG, "png", sep="."))
		png( pngfile, width=900, height=600)
		plotRIPpeak( wiggles, tbl, ripRow=i, label=sampleID, scaledReadCount=scaledReadCount,
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


`RIPpeaksToExcel` <- function( sampleID, optionsFile="Options.txt", results.path=NULL, p.value=0.05) {

	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.null( results.path)) results.path <- getOptionValue( optionsFile, "results.path", notfound=".")
	ripfile <- file.path( results.path, "RIPpeaks", sampleID, 
				paste( sampleID, prefix, "RIPpeaks.txt", sep="."))
	tbl <- read.delim( ripfile, as.is=T)

	csvfile <- sub( "txt$", "Excel.csv", ripfile)

	# trim to only keep the better P value peaks
	keep <- which( tbl$P_Value <= p.value)
	if ( length(keep) < 100) keep <- 1 : min(100,nrow(tbl))
	tbl <- tbl[ keep, ]

	excelterm <- paste( tbl$GENE_ID, round( tbl$Center), sep="_")
	
	tbl2 <- data.frame( "GENE_LINK"=excelterm, tbl, stringsAsFactors=F)
	out <- addExcelLinksToTable( tbl2, linkColumnNames="GENE_LINK", linkPaths="RIP.plots")

	# turn the DNA into fasta format...
	if ( "DNA" %in% colnames(out)) {
		todo <- which( out$DNA != "")
		newdna <- paste( ">", out$GENE_ID[todo], "\n", out$DNA[todo], sep="")
		out$DNA[todo] <- newdna
	}

	write.table( out, csvfile, sep=",", quote=T, qmethod="double", row.names=F)
	return( out)
}


`RIPpeaksAddDNA` <- function( sampleID, optionsFile="Options.txt", results.path=NULL, p.value=0.05) {

	prefix <- getCurrentSpeciesFilePrefix()
	gmap <- getCurrentGeneMap()
	if ( is.null( results.path)) results.path <- getOptionValue( optionsFile, "results.path", notfound=".")
	ripfile <- file.path( results.path, "RIPpeaks", sampleID, 
				paste( sampleID, prefix, "RIPpeaks.txt", sep="."))
	tbl <- read.delim( ripfile, as.is=T)
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

		myfrom <- tbl$Center[i] - 50
		myto <- tbl$Center[i] + 50
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
	write.table( tbl, ripfile, sep="\t", quote=F, row.names=F)

	ndone <- sum( nchar( myDNA) > 0)
	cat( "\nAppended DNA fragments for ", ndone, " RIP peaks.")
	return( NULL)
}


`plotMultiRIPpeaks` <- function( sampleIDset=NULL, annotationFile="Annotation.txt", 
				optionsFile="Options.txt", result.path=NULL, 
				speciesID=NULL, Ngenes=100, label="", ...) {

	annT <- readAnnotationTable( annotationFile)
	if ( is.null( sampleIDset)) {
		sampleIDset <- annT$SampleID
	} else {
		annT <- subset( annT, SampleID %in% sampleIDset)
		sampleIDset <- annT$SampleID
	}

	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) results.path <- getOptionValue( optT, "results.path", notfound=".")
	target <- getAndSetTarget( optionsFile, verbose=T)
	if ( ! is.null( speciesID)) setCurrentSpecies( speciesID)
	speciesID <- getCurrentSpecies()
	prefix <- getCurrentSpeciesFilePrefix()

	# gather up all those wiggle tracks
	WIGlist <- loadMultipleWIGs( sampleIDset, speciesID, results.path)

	# gather up all those RIP peak files
	prefix <- getCurrentSpeciesFilePrefix()
	PEAKlist <- vector( mode="list")
	allGenes <- vector()
	for ( i in 1:length(sampleIDset)) {
		sampleID <- sampleIDset[i]
		ripfile <- file.path( results.path, "RIPpeaks", sampleID, 
				paste( sampleID, prefix, "RIPpeaks.txt", sep="."))
		if ( file.exists( ripfile)) {
			tbl <- read.delim( ripfile, as.is=T)
			n <- min( Ngenes, nrow(tbl))
			PEAKlist[[i]] <- tbl[ 1:n, ]
			allGenes <- c( allGenes, tbl$GENE_ID[1:n])
		} else {
			cat( "\nRIP peak file not found: ", basename(ripfile), "   Skipped..")
		}
	}

	allGenes <- sort( unique( allGenes))
	outpath <- file.path( results.path, "RIPpeaks", "Multi_Sample.plots")
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

