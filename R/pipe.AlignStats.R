# pipe.AlignStats.R

`pipe.AlignStats` <- function( sampleID, annotationFile="Annotation.txt",
				optionsFile="Options.txt", banner="", 
				mode=c("normal","QuickQC", "QuickQC.ribo","splicing", "QuickQC.splicing"), 
				chunkSize=500000, maxReads=NULL, pause=0, results.path=NULL,
				fastqFile=NULL, what=NULL, plot=TRUE) {

	mode <- match.arg( mode)

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	target <- getAndSetTarget( optionsFile, sampleID=sampleID, annotationFile=annotationFile)
	if ( is.null( results.path)) {
		resultsPath <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	} else {
		resultsPath <- results.path
	}

	statsPath <- file.path( resultsPath, "AlignStats", sampleID)
	if ( ! file.exists( statsPath)) dir.create( statsPath, recursive=TRUE)

	if ( is.null( fastqFile)) {
		rawFastq <- getRawFastqFileNames( sampleID, annotationFile, optionsFile, verbose=F)
		fqFileIn <- rawFastq$files
		asMatePairs <- rawFastq$asMatePairs
	} else {
		fqFileIn <- fastqFile
		asMatePairs <- FALSE
	}

	# get the size of the original file...
	if ( is.null( maxReads)) {
		nLines <- getFileLineCount(fqFileIn[1], sampleID)
		if ( nLines < 1) return(NULL)
		totalReads <- round( nLines / 4)
	} else {
		totalReads <- maxReads
	}

	# do all the riboClearing, then all the genomic...
	if ( mode %in% c( "normal", "QuickQC.ribo")) {

		alignPhase <- "RiboClear"
		fileSuffix <- "ribo.converted.bam"
		subfolder <- "riboClear"
		filein <- paste( sampleID, fileSuffix, sep=".")
		if ( mode != "normal") filein <- paste( sampleID, "QuickQC", fileSuffix, sep=".")
		filein <- file.path( resultsPath, subfolder, filein)
		if ( file.exists( filein)) {
			cat( "\n\nRiboCleared Reads:")
			useBanner <- paste( alignPhase, ":  ", banner, sep="")
			calcAlignStats( filein, sampleID, reload=TRUE, alignPhase=alignPhase, totalReads=totalReads, 
					banner=useBanner, statsPath=statsPath, chunkSize=chunkSize, 
					maxReads=maxReads, pause=pause, what=what, plot=plot)
		}
	}  # end of ribo stage

	if ( mode == "QuickQC.ribo") return(NULL)

	if ( mode %in% c( "normal", "QuickQC")) {
		alignPhase <- "Genomic"
		fileSuffix <- "genomic.bam"
		subfolder <- "align"
		filein <- paste( sampleID, fileSuffix, sep=".")
		if ( mode == "QuickQC") filein <- paste( sampleID, "QuickQC", fileSuffix, sep=".")
		filein <- file.path( resultsPath, subfolder, filein)
		if ( file.exists( filein)) {
			cat( "\n\nGenomic Reads:")
			useBanner <- paste( alignPhase, ":  ", banner, sep="")
			calcAlignStats( filein, sampleID, reload=TRUE, alignPhase=alignPhase, totalReads=totalReads, 
					banner=useBanner, statsPath=statsPath, chunkSize=chunkSize, 
					maxReads=maxReads, pause=pause, what=what, plot=plot)
		} else { 
			cat( "\nBAM file not found:  ", filein, "\nNo Genomic 'AlignStats' generated..")
		}
	}

	# there may be splices to do...
	if (TRUE) {
		cat( "\nBypassing 'Splice Align Stats' for now..  broken..")
	} else {
	    if ( mode %in% c( "normal", "splicing", "QuickQC.splicing")) {

		alignPhase <- "Splice"
		fileSuffix <- "splice.bam"
		subfolder <- "splicing"
		filein <- paste( sampleID, fileSuffix, sep=".")
		if ( mode == "QuickQC") filein <- paste( sampleID, "QuickQC", fileSuffix, sep=".")
		filein <- file.path( resultsPath, subfolder, filein)
		if ( file.exists( filein)) {
			cat( "\n\nSplice Reads:")
			useBanner <- paste( alignPhase, ":  ", banner, sep="")
			calcAlignStats( filein, sampleID, reload=TRUE, alignPhase=alignPhase, 
					totalReads=totalReads, banner=useBanner, 
					statsPath=statsPath, chunkSize=chunkSize, maxReads=maxReads,
					pause=pause, what=what, plot=plot)
		}
	    }
	}

	# only do the pie for 'normal' mode
	if ( mode != "normal") return(NULL)

	ans <- pipe.AlignmentPie( sampleID, annotationFile=annotationFile, optionsFile=optionsFile,
				banner=banner, mode=mode)

	return( ans)
}


`dispatch.AlignStats` <- function( sampleID, annotationFile="Annotation.txt",
					optionsFile="Options.txt", banner="", mode="normal",
					chunkSize=250000, maxReads=NULL, pause=0, results.path=NULL,
					fastqFile=NULL) {

	commandLine <- paste( "checkX11( width=9, height=7, xpos=20, ypos=20, bg='white'); ",
				" pipe.AlignStats( sampleID=\"", 
				sampleID, "\", annotationFile=\"", annotationFile, 
				"\", optionsFile=\"", optionsFile, 
				"\", banner=\"", banner, 
				"\", mode=\"", mode, "\", chunkSize=", as.integer(chunkSize),
				", maxReads=", if (is.null(maxReads)) "NULL" else as.integer(maxReads), 
				", pause=", as.integer(pause), 
				", results.path=", if (is.null(results.path)) "NULL" else 
					paste("\"",results.path,"\"",sep=""), 
				", fastqFile=", if (is.null(fastqFile)) "NULL" else 
					paste("\"",fastqFile,"\"",sep=""), 
				" )", sep="")

	if ( mode == "normal") {
		logFile=paste( sampleID, "alignStats.log.txt", sep=".")
	} else {
		logFile=paste( sampleID, mode, "alignStats.log.txt", sep=".")
	}
	spawnRcmd( commandLine, Rpackage="DuffyNGS", logFile=logFile)

	return()
}

