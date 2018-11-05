# pipe.AlignToWig.R


`pipe.AlignToWig` <- function( sampleID, annotationFile="Annotation.txt", optionsFile="Options.txt",
				results.path=NULL, dataType="RNA-seq", 
				mode=c("normal","QuickQC","spliceOnly")) {


	mode <- match.arg( mode)
	if (TRUE) {
		cat( verboseOutputDivider)
		cat( "\nStarting pipe 'AlignToWig' on Sample:     ", sampleID, "\nMode: ", mode, "\n")
	}

	optT <- readOptionsTable( optionsFile)
	target <- getAndSetTarget( optionsFile, sampleID=sampleID, annotationFile=annotationFile)
	readBufferSize <- as.numeric( getOptionValue( optT, "readBufferSize", notfound="1000000"))
	if ( multicore.currentCoreCount() < 2) {
		nCores <- as.integer( getOptionValue( optT, "nCores", notfound=4))
		multicore.setup( nCores)
	}

	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	readSense <- getReadSense( sampleID, annotationFile)
	reload <- TRUE

	# QuickQC mode has a different filename...
	if ( mode == "QuickQC") {
	    mySampleIDs <- getSamplePairIDs( sampleID, annotationFile)

	    for ( sample in mySampleIDs) {
		filein <- paste( sample, "QuickQC.splice.converted.bam", sep=".")
		filein <- file.path( results.path, "splicing", filein)
		readSense <- getReadSense( sample, annotationFile)
		if ( file.exists( filein)) {
			sortfile <- BAM.verifySorted( filein, index=FALSE, threads=1)
			if ( is.null(sortfile)) break
			alignToWig( sortfile, sampleID, readSense=readSense, reload=reload, 
				annotationFile=annotationFile, optionsFile=optionsFile, results.path=results.path, 
				readBufferSize=readBufferSize, dataType=dataType)
			reload <- FALSE
		}
	    }

	    for ( sample in mySampleIDs) {
		filein <- paste( sample, "QuickQC.genomic.bam", sep=".")
		filein <- file.path( results.path, "align", filein)
		readSense <- getReadSense( sample, annotationFile)
		sortfile <- BAM.verifySorted( filein, index=FALSE)
		if ( is.null(sortfile)) break
		alignToWig( sortfile, sampleID, readSense=readSense, reload=reload, 
				annotationFile=annotationFile, optionsFile=optionsFile, results.path=results.path, 
				readBufferSize=readBufferSize, dataType=dataType)
		reload <- FALSE
	    }
	    return()
	}

	# there is a chance we are doing strand specific paired end data
	# as of march 2016, no longer doing as separate _1, _2 files

	# now that the alignments got sorted, we don't need to load in small to big order...
	# do the genomic first, in case the splices are not done converting yet...
	if ( mode == "spliceOnly") reload <- FALSE
	if ( mode == "normal") {
	    	filein <- paste( sampleID, "genomic.bam", sep=".")
	    	filein <- file.path( results.path, "align", filein)
	    	readSense <- getReadSense( sampleID, annotationFile)
	    	sortfile <- BAM.verifySorted( filein, index=FALSE)
		if ( ! is.null(sortfile)) {
	    		nbuf <- alignToWig( sortfile, sampleID, readSense=readSense, reload=reload, 
					annotationFile=annotationFile, optionsFile=optionsFile, results.path=results.path, 
					readBufferSize=readBufferSize, dataType=dataType)
	    		if ( nbuf > 0) reload <- FALSE
		}
	}
	
	filein <- paste( sampleID, "splice.converted.bam", sep=".")
	filein <- file.path( results.path, "splicing", filein)
	readSense <- getReadSense( sampleID, annotationFile)
	if ( file.exists( filein)) {
		sortfile <- BAM.verifySorted( filein, index=FALSE, threads=1)
		if ( ! is.null(sortfile)) {
			nbuf <- alignToWig( sortfile, sampleID, readSense=readSense, reload=reload, 
					annotationFile=annotationFile, optionsFile=optionsFile, results.path=results.path, 
					readBufferSize=readBufferSize, dataType=dataType)
			if ( nbuf > 0) reload <- FALSE
		}
	}

	return()
}
