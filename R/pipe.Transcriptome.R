# pipe.Transcriptome.R

# turn the Wiggle Track pipeups into gene expression results

`pipe.Transcriptome` <- function( sampleID, annotationFile="Annotation.txt", optionsFile="Options.txt", 
				speciesID=NULL, results.path=NULL, dataType=NULL,
				altGeneMap=NULL, altGeneMapLabel=NULL, loadWIG=FALSE, verbose=TRUE,
				mode=c("normal","QuickQC"), exonsOnly=NULL)
{

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nStarting pipe 'Transcriptome' on Sample:     ", sampleID, "\t", altGeneMapLabel, "\n")
	}
	mode <- match.arg( mode)
	startTotalT <- proc.time()
	grandTotalReads <- 0

	optT <- readOptionsTable( optionsFile)
	target <- getAndSetTarget( optionsFile, sampleID=sampleID, annotationFile=annotationFile)
	allSpecies <- getCurrentTargetSpecies()
	if ( ! is.null(speciesID)) allSpecies <- speciesID

	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	if ( ! file.exists( results.path)) dir.create( results.path, showWarnings=FALSE)

	# other annotation facts we need...
	origSampleID <- originalSamplePairID( sampleID, annotationFile)
	sampleID <- origSampleID
	annT <- readAnnotationTable( annotationFile)

	if ( is.null(dataType)) dataType <- getAnnotationValue( annT, origSampleID, "DataType", notfound="RNA-seq")

	if ( mode == "normal" && dataType == "ChIP-seq") {
		pipe.ChIPpeaks( origSampleID, annotationFile, optionsFile, speciesID=speciesID, 
				results.path=results.path, loadWIG=loadWIG, verbose=FALSE)
		return()
	}

	strandSpecific <- getAnnotationTrue( annT, origSampleID, "StrandSpecific", notfound=FALSE)
	useBothStrands <- ! strandSpecific
	keepIntergenics <- getAnnotationTrue( annT, origSampleID, "KeepIntergenics", notfound=TRUE)
	if ( is.null( exonsOnly)) {
		exonsOnly <- getAnnotationTrue( annT, origSampleID, "ExonsOnly", notfound=FALSE)
	} else {
		cat( "\nExplicit argument:  \t", sampleID, ": \t", exonsOnly, sep="")
	}

	# during 'altGeneMap' runs, skip over species that are not covered...
	doingAltGeneMap <- FALSE
	if ( ! is.null( altGeneMap)) {
		doingAltGeneMap <- TRUE
		altSpeciesSet <- unique.default( getSpeciesFromSeqID( altGeneMap$SEQ_ID))
	}

	needLoadWIG <- loadWIG

	# do each species all the way through
	for( speciesID in allSpecies) {
	
	if ( doingAltGeneMap && (!( speciesID %in% altSpeciesSet))) next

	cat( "\n\nExtracting results for Species:  ", speciesID,"\n")
	startT <- proc.time()
	gc()
	setCurrentSpecies( speciesID)
	speciesPrefix <- getCurrentSpeciesFilePrefix()

	# find and/or load our WiggleBin data structure
	fileOutWIG <- file.path( results.path, "wig", paste( sampleID, speciesPrefix, "WIG.rda", sep="."))

	if ( !file.exists( fileOutWIG) || needLoadWIG) {
		
		cat( "\nLoading Wiggle Track data structure from alignments...")
		pipe.AlignToWig( sampleID, annotationFile=annotationFile, optionsFile=optionsFile,
				results.path=results.path, dataType=dataType, mode=mode)

		#get that data back when done... object 'wiggles'
		load( file=fileOutWIG)

		if ( mode == "normal") {
			cat( "\nCalculating Wiggle Track Metrics for Gene and Strand Specificity...\n")
			binMetricsFile <- paste( sampleID, speciesPrefix, "Metrics.txt", sep=".")
			binMetricsPath <- file.path( results.path, "summary")
			if ( ! file.exists( binMetricsPath)) dir.create( binMetricsPath, recursive=T)
			binMetricsFile <- file.path( binMetricsPath, binMetricsFile)
			calcWIGmetrics( wiggles, asDataFrame=FALSE, logFile=binMetricsFile)
		}
		needLoadWIG <- FALSE

	} else {
		cat( "\nLoading pre-existing Wiggles data...")
		load( file=fileOutWIG)
	}

	# we could be given a non-standard geneMap to use...
	if ( is.null( altGeneMap)) {

		# regular...
		fileOutTrans <- paste( sampleID, speciesPrefix, "Transcript.txt", sep=".")
		pathOut <- file.path( results.path, "transcript")
		if ( ! file.exists( pathOut)) dir.create( pathOut, showWarnings=FALSE)
		fileOutTrans <- file.path( pathOut, fileOutTrans)
		trans <- calcWigTranscriptome( wiggles, useBothStrands=useBothStrands, 
					keepIntergenics=keepIntergenics, exonsOnly=exonsOnly,
					fileout=fileOutTrans)

	} else {
		# alternate...
		gmap <- altGeneMap
		if ( is.character(gmap)) {
			gmap <- read.delim( file=altGeneMap, as.is=TRUE)
		}
		if ( ! all( c("GENE_ID", "SEQ_ID") %in% colnames(gmap))) 
			stop( paste( "pipe.Transcriptome: invalid alternate geneMap",
			"does not have required GENE_ID, SEQ_ID columns"))
		if ( is.null( altGeneMapLabel) || base::nchar( altGeneMapLabel) < 1) 
				stop( "pipe.Transcriptome: missing file name text term 'altGeneMapLabel' ")

		# only do it for the right species...
		if ( getSpeciesFromSeqID( gmap$SEQ_ID[1]) == speciesID) {

			cat( "\n\nDoing alternate gene map transcriptome:  ", altGeneMapLabel, "\n")
			fileOutTrans <- paste( sampleID, speciesPrefix, altGeneMapLabel, "Transcript.txt", sep=".")

			# use a sub-folder for alternate transcripts...
			pathOut <- file.path( results.path, "transcript", altGeneMapLabel)
			if ( ! file.exists( pathOut)) dir.create( pathOut, showWarnings=FALSE)
			fileOutTrans <- file.path( pathOut, fileOutTrans)

			trans <- calcWigTranscriptome( wiggles, geneMap=gmap, useBothStrands=useBothStrands, 
					keepIntergenics=TRUE, exonsOnly=FALSE, fileout=fileOutTrans)
		} else {
			cat( "\nThis species not in Alternate Gene Map... Skipping. ")
		}
	}

	reads <- WIG_getTotalRawReads( wiggles)
	grandTotalReads <- grandTotalReads + reads$Unique + reads$Multi

	myTime <- elapsedProcTime( startT, proc.time(), N=(reads$Unique + reads$Multi))
	if (verbose) {
		cat( "\n\nFinished Species:  ", speciesID)
		cat( "\n\nSpecies Timing Stats: \n")
		print( myTime)
		gc()
	}

	}  # end of all speceies...

	myTime <- elapsedProcTime( startTotalT, proc.time(), N=grandTotalReads)
	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\n\nFinished pipe.Transcriptome:  ", sampleID, "\tSpecies set:  ", allSpecies)
		cat( "\n\nSample Timing Stats: \n")
		print( myTime)
	}

	return( sampleID)
}



`dispatch.TranscriptPlusHTML` <- function( sampleID, annotationFile="Annotation.txt",
					optionsFile="Options.txt", banner="", mode="normal",
					maxReads=NULL, pause=0, results.path=NULL) {

	commandLine <- paste( "checkX11( width=8, height=6, xpos=20, ypos=20, bg='white'); ",
				" pipe.TranscriptPlusHTML( sampleID=\"", 
				sampleID, "\", annotationFile=\"", annotationFile, 
				"\", optionsFile=\"", optionsFile, 
				"\", banner=\"", banner, 
				"\", mode=\"", mode, "\"", 
				", maxReads=", if (is.null(maxReads)) "NULL" else as.integer(maxReads), 
				", pause=", as.integer(pause), 
				", results.path=", if (is.null(results.path)) "NULL" else 
					paste("\"",results.path,"\"",sep=""), 
				" )", sep="")

	spawnRcmd( commandLine, Rpackage="DuffyNGS", logFile=paste( sampleID, "transcript.log.txt", sep="."))

	return()
}


`pipe.TranscriptPlusHTML` <- function( sampleID, annotationFile="Annotation.txt",
				optionsFile="Options.txt", banner="", mode="normal", loadWIG=FALSE,
				maxReads=NULL, pause=0, results.path=NULL, triggerWarning=NULL, ...) {

	# get a bit of info for setting up
	optT <- readOptionsTable( optionsFile)
	target <- getAndSetTarget( optionsFile, sampleID=sampleID, annotationFile=annotationFile)
	allSpecies <- getCurrentTargetSpecies()

	# do the transcriptomes
	pipe.Transcriptome( sampleID, annotationFile, optionsFile, speciesID=NULL,
			results.path=results.path, loadWIG=loadWIG, mode=mode)
	
	# do we want a strand specific check?
	isStranded <- getAnnotationTrue( annotationFile, key=sampleID, columnArg="StrandSpecific", verbose=F)
	if ( isStranded && is.null(triggerWarning)) {
		strandAns <- pipe.StrandVerify( sampleID, annotationFile=annotationFile, optionsFile=optionsFile,
						speciesID=allSpecies[1], results.path=results.path)
		# don't treat these 2 P-values exactly as equals...
		#worstPval <- max( strandAns$NonGene_Pvalue, strandAns$StrandValue_Pvalue, na.rm=T)
		ng.Pval <- strandAns$NonGene_Pvalue
		if (ng.Pval > 0.25) triggerWarning <- "Non-Genes among highest expressers.  Read Sense may be wrong!"
		strandCC.Pval <- strandAns$StrandValue_Pvalue
		if (strandCC.Pval > 0.05) triggerWarning <- "Read Strand Verification failed.  Read Sense may be wrong!"
	}

	# now make some gene plots
	nGenes <- 25
	for( s in allSpecies) {
		pipe.TranscriptToHTML( sampleID, annotationsFile, optionsFile,
				speciesID=s, results.path=results.path, 
				pause=pause, N=nGenes, label=banner, triggerWarning=triggerWarning, ...)
	}
	return()
}
