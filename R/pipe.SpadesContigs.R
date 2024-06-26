# pipe.SPAdesContigs.R


`pipe.SpadesContigs` <- function( sampleID, fastqSource=NULL,
		annotationFile="Annotation.txt", optionsFile="Options.txt", 
		results.path=NULL, spades.path=dirname(Sys.which("spades.py")),
		spades.mode=c("isolate","rna","meta"), spades.args="", pairedEnd=NULL,
		kmerSizes=NULL, makePeptides=TRUE, folderName=NULL, keyword="Spades", verbose=TRUE) {

	require( Biostrings)
	if (version$major == "4" && as.numeric( version$minor) >= 4) require( pwalign)

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\n\nStarting 'SPAdesContigs' for Sample:     ", sampleID, "\n\n")
	}

	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=FALSE)
	}

	spades.mode <- match.arg( spades.mode)

	# spades wants a folder for writing all its results
	outpath <- file.path( results.path, "SpadesContigs", sampleID)
	if ( ! is.null( folderName)) {
		outpath <- file.path( results.path, "SpadesContigs", sampleID, folderName)
	}
	if ( ! file.exists( outpath)) dir.create( outpath, recursive=T)

	# what fastq files go into the tool?
	# default is the full raw fastq data before any alignments are done
	if ( is.null( fastqSource)) {
		fastqPath <- getOptionValue( optT, "fastqData.path", verbose=FALSE)
		fastqFile <- getAnnotationValue( annT, key=sampleID, columnArg="Filename", verbose=FALSE)
		fastqFile <- strsplit( fastqFile, split=", *")[[1]]
		fastqFile <- file.path( fastqPath, fastqFile)
		if (is.null( pairedEnd)) pairedEnd <- getAnnotationTrue( annT, key=sampleID, columnArg="PairedEnd", notfound=FALSE, verbose=FALSE)
		doPairedEnd <- (pairedEnd && length(fastqFile) == 2)
	} else { 
		# otherwise, some files from the 'fastq' results subfolder for this sample
		fastqPath <- file.path( results.path, "fastq")
		fastqFile <- paste( sampleID, fastqSource, "fastq.gz", sep=".")
		fastqFile <- file.path( fastqPath, fastqFile)
		doPairedEnd <- FALSE
		if (verbose) {
			cat( "\nUsing Fastq files:", fastqFile, sep="\n")
		}
	}
	# make sure those files exist
	fqExist <- file.exists( fastqFile)
	if ( ! all( fqExist)) {
		cat( "\nError:  some FASTQ files not found:\n")
		cat( fastqFile[ ! fqExist], "\n")
		return(NULL)
	}

	ans <- makeSpadesContigs( fastqFile, outpath=outpath, spades.path=spades.path, 
			spades.mode=spades.mode, kmerSizes=kmerSizes, 
			spades.args=spades.args, doPairedEnd=doPairedEnd, verbose=verbose)

	if ( makePeptides) {
		avgPepLen <- makeSpadesPeptides( sampleID, outpath=outpath, keyword=keyword, verbose=verbose)
		ans$AvgPeptide <- round(avgPepLen)
	}

	return( ans)
}


cleanupSpadesFiles <- function( path="results/SpadesContigs", verbose=T) {

	# delete any large unneeded files left behind by Velvet...
	folderPatterns <- c( "^misc$", "^tmp$", "^pipeline_state$", "^K[1-9]+$")
	filePatterns <- c( "^assembly_graph")

	totalBytes <- 0
	nfiles <- nfolders <- 0

	for (patt in folderPatterns) {
		files <- dir( path, pattern=patt, recursive=T, full.name=T, include.dir=T)
		if ( length( files) < 1) next
		# reorder, so all the regulare files come first
		ord <- order( file.info(files)$isdir, -nchar(files))
		files <- files[ ord]
		while ( length(files)) {
			f <- files[1]
			finfo <- file.info( f)
			bytes <- finfo$size
			nfiles <- nfiles + 1
			if (finfo$isdir) {
				# see if we delete it, or drill down first
				moreFiles <- dir( f, full.name=T, include.dir=T)
				if ( length(moreFiles)) {
					files <- c( moreFiles, files)
					next
				}
				nfolders <- nfolders + 1
				file.remove( f)
			}
			if (verbose) cat( "\r", nfiles, "  Size: ", bytes, "  File: ", f)
			totalBytes <- totalBytes + bytes
			file.delete( f)
			NF <- length( files)
			files <- if (NF > 1) files[ 2:NF] else character(0)
		}
		if (verbose) cat( "\n")
	}

	for (patt in filePatterns) {
		files <- dir( path, pattern=patt, recursive=T, full.name=T)
		if ( length( files) < 1) next
		for ( f in files) {
			bytes <- file.info( f)$size
			nfiles <- nfiles + 1
			if (verbose) cat( "\r", nfiles, "  Size: ", bytes, "  File: ", f)
			totalBytes <- totalBytes + bytes
			file.delete( f)
		}
		if (verbose) cat( "\n")
	}

	if (verbose) cat( "\nDeleted_Folders: ", nfolders, "\tDeleted_Files: ", nfiles, "\tDeleted Bytes: ", 
			prettyNum( totalBytes, big.mark=","), "\n")
}

