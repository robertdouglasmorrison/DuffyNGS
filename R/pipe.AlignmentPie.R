# pipe.AlignmentPie.R

`pipe.AlignmentPie` <- function( sampleID, annotationFile="Annotation.txt", optionsFile="Options.txt",
				fastqFile=NULL, banner="", mode=c("normal","QuickQC"), results.path=NULL,
				totalRawReads=NULL, useUSR=TRUE) {

	# get the needed paths, etc.
	optT <- readOptionsTable( optionsFile)

	mode <- match.arg( mode)

	if ( is.null(results.path)) {
		resultsPath <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
		if ( mode == "QuickQC") {
			resultsPath <- "QuickQC"
		}
	} else {
		resultsPath <- results.path
	}

	# we may need to treat the paired end sampleIDs special...
	originalSampleID <- originalSamplePairID( sampleID)
	pairSampleIDs <- getSamplePairIDs(originalSampleID)

	statsPath <- file.path( resultsPath, "AlignStats", sampleID)
	if ( mode == "normal") {
		statsPath <- file.path( resultsPath, "AlignStats", originalSampleID)
	}
	if ( ! file.exists( statsPath)) dir.create( statsPath, recursive=T, showWarnings=F)
	USRpath <- file.path( resultsPath, "USR")
	fastqPath <- getOptionValue( optT, "fastqData.path", notfound=".")

	# allow for a possible paired end data
	if ( is.null( fastqFile)) {
		fqFileIn <- getAnnotationValue( annotationFile, key=originalSampleID, columnArg="Filename")
	} else {
		fqFileIn <- fastqFile
	}

	fqFileSet <- strsplit( fqFileIn, split=", *")[[1]]
	fqFile <- fqFileSet[ match( sampleID, pairSampleIDs, nomatch=1)]

	# if not given the fastq file explicitly, turn it to its full pathname
	if ( is.null( fastqFile)) {
		fqFile <- file.path( fastqPath, fqFile)
		fqFile <- allowCompressedFileName( fqFile)
	}

	# get the size of the original file...
	if ( is.null( totalRawReads)) {
		nLines <- getFileLineCount(fqFile, sampleID)
		if ( nLines < 1) return(NULL)
		nRawReads <- round( nLines / 4)
	} else {
		nRawReads <- totalRawReads
	}
	cat( "\n  Raw: ", nRawReads)

	if (useUSR) {
		USRfile <- file.path( USRpath, paste( sampleID, "USR.rda", sep="."))
		if ( !file.exists( USRfile)) {
			pipe.USR( sampleID, annotationFile, optionsFile, resultsPath=resultsPath, mode=mode)
		}
		load( USRfile)
	}

	# the names of some files and what they contain depend on mode...
	nNoHits <- nMoreK <- nMoreKribo <- 0
	nUnique <- nMulti <- 0
	nUniqueRibo <- nMultiRibo <- 0
	nUniqueSplice <- nMultiSplice <- 0
	nPolyN <- 0
	adapterCounts <- c( 0, 0)
	names(adapterCounts) <- c("ForwardAdapter", "ReverseAdapter")

	if ( mode == "normal") {

		# the files on NoHit or MoreK are fastqs -- 4 per read
		nohitsFile <- file.path( resultsPath, "fastq", paste( sampleID, "noHits.fastq", sep="."))
		nohitsFile <- allowCompressedFileName( nohitsFile)
		if ( file.exists( nohitsFile)) nNoHits <- round( getFileLineCount( nohitsFile, sampleID) / 4)
		cat( "  NoHits: ", nNoHits)
		morekFile <- file.path( resultsPath, "fastq", paste( sampleID, "moreK.fastq", sep="."))
		morekFile <- allowCompressedFileName( morekFile)
		if ( file.exists( morekFile)) nMoreK <- round( getFileLineCount( morekFile, sampleID) / 4)
		cat( "  MoreK: ", nMoreK)
		morekFile <- file.path( resultsPath, "fastq", paste( sampleID, "moreK.ribo.fastq", sep="."))
		morekFile <- allowCompressedFileName( morekFile)
		if ( file.exists( morekFile)) nMoreKribo <- round( getFileLineCount( morekFile, sampleID) / 4)
		cat( "  MoreKribo: ", nMoreKribo)
	
		# all reads are mixed together in the BAM files, so count them out
		bamFile <- file.path( resultsPath, "riboClear", paste( sampleID, "ribo.converted.bam", sep="."))
		if ( file.exists( bamFile)) {
			cat( "\n  Counting RiboCleared in BAM file")
			ans <- countUniqueAndMultiReads( bamFile)
			nUniqueRibo <- ans$Unique
			cat( "  RiboUnique: ", nUniqueRibo)
			nMultiRibo <- ans$Multi
			cat( "  RiboMulti: ", nMultiRibo)
		}
		bamFile <- file.path( resultsPath, "align", paste( sampleID, "genomic.bam", sep="."))
		if ( file.exists( bamFile)) {
			cat( "\n  Counting Genomics in BAM file")
			ans <- countUniqueAndMultiReads( bamFile)
			nUnique <- ans$Unique
			cat( "  GenomeUnique: ", nUnique)
			nMulti <- ans$Multi
			cat( "  GenomeMulti: ", nMulti)
		}
		bamFile <- file.path( resultsPath, "splicing", paste( sampleID, "splice.bam", sep="."))
		if ( file.exists( bamFile)) {
			cat( "\n  Counting Splices in BAM file")
			ans <- countUniqueAndMultiReads( bamFile)
			nUniqueSplice <- ans$Unique
			cat( "  SpliceUnique: ", nUniqueSplice)
			nMultiSplice <- ans$Multi
			cat( "  SpliceMulti: ", nMultiSplice)
		}

	} else if ( mode == "QuickQC") {

		nohitsFile <- file.path( resultsPath, "fastq", paste( sampleID, "QuickQC.noHits.fastq.gz", sep="."))
		if ( file.exists( nohitsFile)) nNoHits <- round( getFileLineCount( nohitsFile, sampleID) / 4)

		# uniques and multis are alignments...
		bamFile <- file.path( resultsPath, "riboClear", paste( sampleID, "QuickQC.ribo.converted.bam", sep="."))
		if ( file.exists( bamFile)) {
			cat( "\n  Counting RiboCleared in BAM file")
			ans <- countUniqueAndMultiReads( bamFile)
			nUniqueRibo <- ans$Unique
			cat( "  RiboUnique: ", nUniqueRibo)
			nMultiRibo <- ans$Multi
			cat( "  RiboMulti: ", nMultiRibo)
		}
		bamFile <- file.path( resultsPath, "align", paste( sampleID, "QuickQC.genomic.bam", sep="."))
		if ( file.exists( bamFile)) {
			cat( "\n  Counting Genomics in BAM file")
			ans <- countUniqueAndMultiReads( bamFile)
			nUnique <- ans$Unique
			cat( "  GenomeUnique: ", nUnique)
			nMulti <- ans$Multi
			cat( "  GenomeMulti: ", nMulti)
		}
		bamFile <- file.path( resultsPath, "splicing", paste( sampleID, "QuickQC.splice.bam", sep="."))
		if ( file.exists( bamFile)) {
			cat( "\n  Counting Splices in BAM file")
			ans <- countUniqueAndMultiReads( bamFile)
			nUniqueSplice <- ans$Unique
			cat( "  SpliceUnique: ", nUniqueSplice)
			nMultiSplice <- ans$Multi
			cat( "  SpliceMulti: ", nMultiSplice)
		}
	} else {
		stop ( "unknown 'AlignmentPie' mode...")
	}

	# add in USR details if given
	if ( exists( "USR_nDropPolyN")) {
		# adapter counts and PolyN are based on the USR counts, but these may be a subset
		scaleFactor <- 1.0 / USR_pctCoverage
		nPolyN <- USR_nDropPolyN * scaleFactor
		if ( exists( "USR_AdapterHits") && !is.null(USR_AdapterHits)) {
			adapterCounts <- USR_AdapterHits * scaleFactor
		}
		nNoHits <- nNoHits - (sum(adapterCounts) + nPolyN)
	}

	# now we can build the pie, based on Mode
	pieCounts <- c( adapterCounts, nPolyN, nNoHits, nUniqueRibo, nMultiRibo, nUnique, nUniqueSplice, 
			nMulti, nMultiSplice, (nMoreK + nMoreKribo))
	pieNames <- c( names(adapterCounts), "Poly_N", "Other_NoHits", "RiboClear_Unique", "RiboClear_MultiHit",
			"Unique_Genomic", "Unique_Splice", "MultiHit_Genomic", "MultiHit_Splice", "More_K_Hits")
	pieColors <- c("pink", "magenta", "brown", 2, "orange","goldenrod", 3, "seagreen", 5,"dodgerblue", 4) 
	names( pieCounts) <- pieNames

	if ( mode == "normal") {
		mainText <- paste( "Alignment Success Overview:      \n Sample = ", sampleID, "         ", banner, 
			"\n Total Raw Reads = ", formatC( round( nRawReads), format="d", big.mark=","))
		piefile <- file.path( statsPath, paste( sampleID, "Alignment.Pie.png", sep="."))
	} else if ( mode == "QuickQC") {
		mainText <- paste( "Alignment Success Overview:  \n Quick QC:        Sample = ", sampleID, "         ", 
			banner, "\n Total Raw Reads = ", formatC( round( nRawReads), format="d", big.mark=","))
		piefile <- file.path( statsPath, paste( sampleID, "QuickQC.Alignment.Pie.png", sep="."))
	}

	# don't draw any 'tiny' wedges
	keep <- which( pieCounts > (sum(pieCounts)*0.0010))
	pcts <- as.percent( pieCounts, big.value=nRawReads)
	useNames <- paste( pieNames, pcts, sep=" ")

	pie( pieCounts[keep], labels=useNames[keep], col=pieColors[keep], radius=0.9, clockwise=TRUE,
			cex.lab=1.1, init.angle=0, main=mainText)

	# draw it to disk too
	piefile <- file.path( statsPath, paste( sampleID, "Alignment.Pie.png", sep="."))
	png( piefile, width=800, height=600, bg="white")
	pie( pieCounts[keep], labels=useNames[keep], col=pieColors[keep], radius=0.9, clockwise=TRUE,
			cex.lab=1.1, init.angle=0, main=mainText)
	dev.off()

	ans <- ( list( "ReadCounts"=round(pieCounts), "Percents"=pcts))

	cat( "\n\nAlignment Pie:   \t", sampleID, "\n")
	print( ans)

	return( ans)
}
