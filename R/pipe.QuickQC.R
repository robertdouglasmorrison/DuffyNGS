# pipe.QuickQC.R

# do a quick fastq stats, genomic align, align stats, and transcript

`pipe.QuickQC` <- function( sampleID, annotationFile="Annotation.txt",
		optionsFile="Options.txt", results.path="QuickQC", banner="QuickQC", 
		chunkSize=1000000, maxReads=5000000, pause=0, verbose=T, nUSRkeep=1000000,
		mode=c("all", "alignOnly", "pieOnly")) {

	sink( file=paste( sampleID, "QuickQC.log.txt", sep="."), split=TRUE)
	on.exit( sink())

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\n\nStarting 'Quick QC pipeline' for Sample:     ", sampleID, "\n\n")
	}
	gc()

	optTbl <- readOptionsTable( optionsFile)
	annTbl <- readAnnotationTable( annotationFile)
	target <- getAndSetTarget( optionsFile, sampleID=sampleID, annotationFile=annotationFile)
	nCores <- as.integer( getOptionValue( optTbl, "nCores", notfound="4"))
	multicore.setup(nCores)

	if ( ! file.exists( results.path)) dir.create( results.path, recursive=TRUE, showWarnings=FALSE)

	# verify we see the file
	sampleSet <- annTbl$SampleID
	samplePtr <- base::match( sampleID, sampleSet, nomatch=0)
	if ( samplePtr < 1) {
		cat( "\n\nSampleID not in annotaion file:   given: ", sampleID, "  choices: ", sampleSet)
		return()
	}
	fastq.path <- getOptionValue( optTbl, "fastqData.path", notfound=".")
	statsPath <- file.path( results.path, "FastqReadStats")
	fastqFile <- getAnnotationValue( annTbl, key=sampleID, columnArg="Filename")

	# this may be paired end reads...
	pairedEnd <- getAnnotationTrue( annTbl, key=sampleID, columnArg="PairedEnd", notfound=FALSE)
	if (pairedEnd) {
		samplePairIDs <- getSamplePairIDs( sampleID, annotationFile)
		fastqFileSet <- strsplit( fastqFile, split=", *")[[1]]
		# only do half the wanted reads from each pair
		maxReads <- round( maxReads/2)
	} else {
		samplePairIDs <- sampleID
		fastqFileSet <- fastqFile
		# if comma sep list, just use the first for QuickQC...
		fastqFileSet <- strsplit( fastqFile, split=", *")[[1]][1]
	}

	dataType <- getAnnotationValue( annTbl, key=sampleID, columnArg="DataType", notfound="RNA-seq")
	readBufferSize <- as.numeric( getOptionValue( optTbl, "readBufferSize", notfound="1000000"))

	# set up sub folders
	folderSet <- c( "align", "riboClear", "fastq", "splicing", "USR", "CR", "transcript", "wig")
	if ( dataType == "ChIP-seq") folderSet <- c( "align", "ChIPpeaks", "fastq", "USR", "CR", "transcript", "wig")
	if ( dataType == "RIP-seq") folderSet <- c( "align", "riboClear", "RIPpeaks", "fastq", "splicing", "USR", 
							"CR", "transcript", "wig")
	for ( folder in folderSet) {
		subfolder <- file.path( results.path, folder)
		if ( ! file.exists( subfolder)) dir.create( subfolder, recursive=TRUE, showWarnings=FALSE)
	}

	# cludge to allow partial skipping of some steps
	mode <- match.arg( mode)

	for ( ipair in 1:length(samplePairIDs)) {

	sample <- samplePairIDs[ ipair]
	fastqFile <- fastqFileSet[ ipair]
	fastqFile <- file.path( fastq.path, fastqFile)
	filein <- allowCompressedFileName( fastqFile)
	genomicFile <- file.path( results.path, "align", paste( sample, "QuickQC.genomic.bam", sep="."))
	riboFile <- file.path( results.path, "riboClear", paste( sample, "QuickQC.ribo.bam", sep="."))
	riboFileOut <- file.path( results.path, "riboClear", paste( sample, "QuickQC.ribo.converted.bam", sep="."))
	notRiboFile <- file.path( results.path, "fastq", paste( sample, "QuickQC.notRibo.fastq.gz", sep="."))
	notGenomicFile <- file.path( results.path, "fastq", paste( sample, "QuickQC.notGenomic.fastq.gz", sep="."))
	spliceFile <- file.path( results.path, "splicing", paste( sample, "QuickQC.splice.bam", sep="."))
	spliceFileOut <- file.path( results.path, "splicing", paste( sample, "QuickQC.splice.converted.bam", sep="."))
	nohitFile <- file.path( results.path, "fastq", paste( sample, "QuickQC.noHits.fastq.gz", sep="."))

	# delete any pre-existing files
	file.delete( c( riboFile, riboFileOut))
	sortedFile <- sub( "bam$", "sorted.bam", genomicFile)
	file.delete( c( genomicFile, sortedFile))
	file.delete( c( notRiboFile, notGenomicFile))
	sortedFile <- sub( "bam$", "sorted.bam", spliceFileOut)
	file.delete( c( spliceFile, spliceFileOut, sortedFile))
	file.delete( nohitFile)
	USRfile <- file.path( results.path, "USR", paste( sample, "USR.rda", sep="."))
	file.delete( USRfile)

	if (pairedEnd) cat( "\n\nRunning QuickQC on Paired Sample part:  ", sample, "\n")

	cat( verboseOutputDivider)
	cat( "\n\nSample: ", sample, "\n")
	cat( "\nFastq file: ", filein)

	# show the user the type of reads/scores
	ans <- detectFastqReadFormat( filein)
	cat( "\n\nFastq file encoding:  ")
	cat( "\n  Read ID type:       ", ans$readIDtype)
	cat( "\n  Score type:         ", ans$scoreType, "\n\n")


	if ( mode == "all") {
		# start the Fastq Stats as a separate process..
		cat( verboseOutputDivider, "\nSpawning 'FastqReadStats' child process...\n")
		dispatch.FastqReadStats( filein=filein, sampleID=sample, statsPath=statsPath, chunkSize=chunkSize,
					maxReads=maxReads, pause=pause)
		Sys.sleep(2)
	}

	nRawReads <- NULL

	if ( mode != "pieOnly") {

	# set up the files, indexes , etc. for alignments...
	cat( verboseOutputDivider, "\nSetting up for 'Bowtie2' alignments...\n")

	indexPath <- getOptionValue( optTbl, "bowtie2Index.path", notfound=".")
	riboIndex <- getOptionValue( optTbl, "RiboIndex", notfound="Pf.RiboClear_idx")
	riboPolicy <- getOptionValue( optTbl, "RiboAlignmentPolicy", notfound=" ")
	riboMap <- getOptionValue( optTbl, "RiboMap", notfound="Pf.RiboMap.txt")
	spliceIndex <- getOptionValue( optTbl, "SpliceIndex", notfound="Pf.SpliceClear_idx")
	splicePolicy <- getOptionValue( optTbl, "SpliceAlignmentPolicy", notfound=" ")
	spliceMapPrefix <- getOptionValue( optTbl, "SpliceMapPrefix", notfound="spliceMap")

	# do a ribo clear step
	if ( dataType %in% c("RNA-seq","RIP-seq") && riboIndex != "") {

		if ( nchar(riboMap) < 4) {
			cat( "\nError:  option 'RiboMap' is not a valid filename")
			stop( "pipe.QuickQC is aborting..")
		}
		riboMap <- file.path( indexPath, riboMap)
		if ( ! file.exists(riboMap)) {
			cat( "\nError:  option 'RiboMap' file not found.  Tried: ", riboMap)
			stop( "pipe.QuickQC is aborting..")
		}

		cat( verboseOutputDivider, "\nPre-Alignment against Ribo Clearing Index...\n")
		ans <- fastqToBAM( filein, riboFile, sampleID=sample, k=5,
				optionsFile=optionsFile, annotationFile=annotationFile, noHitsFile=notRiboFile,
				alignIndex=riboIndex, alignPolicy=riboPolicy, maxReads=maxReads, verbose=T)
		nRawReads <- ans$RawReads
		nRibos <-  ans$UniqueReads + ans$MultiReads
		nNotRibo <- ans$NoHitReads
		quickFileLineCountRecord( filein, sample, lineCount=(nRawReads*4), readCount=nRawReads)
		quickFileLineCountRecord( notRiboFile, sample, lineCount=(nNotRibo*4), readCount=nNotRibo)
		quickFileLineCountRecord( riboFile, sample, lineCount=nRibos)
		notRiboNeedsDelete <- TRUE
		# if no aligments, delete that BAM file
		if ( nRibos < 1) file.delete( riboFile)
	} else {
		# no ribo clear...
		notRiboFile <- filein
		notRiboNeedsDelete <- FALSE
	}

	# do a genomic step
	cat( verboseOutputDivider, "\nMain Alignment against Genomic Index(s)...\n")
	genomicIndex <- getOptionValue( optTbl, "GenomicIndex", notfound="Pf.Genomic_idx")
	genomicPolicy <- getOptionValue( optTbl, "GenomicAlignmentPolicy", notfound=" ")
	ans <- fastqToBAM( notRiboFile, genomicFile, sampleID=sample, k=5,
			optionsFile=optionsFile, annotationFile=annotationFile, noHitsFile=notGenomicFile,
			alignIndex=genomicIndex, alignPolicy=genomicPolicy, maxReads=maxReads, verbose=T)
	nNotGenomic <- ans$NoHitReads
	quickFileLineCountRecord( notGenomicFile, sample, lineCount=(nNotGenomic*4), readCount=nNotGenomic)
	if ( riboIndex == "") {
		nRawReads <- ans$RawReads
		quickFileLineCountRecord( filein, sample, lineCount=(nRawReads*4), readCount=nRawReads)
	}
	if (notRiboNeedsDelete) file.delete( notRiboFile)
	nGenomic <-  ans$UniqueReads + ans$MultiReads
	quickFileLineCountRecord( genomicFile, sample, lineCount=nGenomic)

	# now that the genomic is done, we can convert the ribos, genomic, for doing align stats...
	if ( dataType %in% c("RNA-seq","RIP-seq") && riboIndex != "") {
		ans <- riboConvert( riboFile, riboFileOut, genomefile=genomicFile, sampleID=sampleID, 
				riboMapFile=riboMap, rawReadCount=nRawReads, 
				readBufferSize=readBufferSize, verbose=verbose)
		quickFileLineCountRecord( riboFileOut, sampleID, lineCount=ans$Alignments)
		cat("\n")
		cat( ans$textSummary)
		if (chunkSize > 0) {
		    cat( verboseOutputDivider, "\nSpawning 'AlignStats' child process on Ribo Cleared reads...\n")
		    dispatch.AlignStats( sampleID=sample, annotationFile=annotationFile,
				optionsFile=optionsFile, mode="QuickQC.ribo", banner=banner, chunkSize=chunkSize,
				maxReads=maxReads, pause=pause, results.path=results.path, fastqFile=filein)
		    Sys.sleep(2)
		}
		file.delete( riboFile)
	}

	# now convert the genomic to append the GeneIDs
	if (verbose) cat( "\n\nAppending GeneID terms to genomic alignments...\n")
	ans <- genomicConvert( genomicFile, sampleID=sampleID, readBufferSize=readBufferSize, 
				rawReadCount=nRawReads, verbose=verbose)
	quickFileLineCountRecord( genomicFile, sampleID, lineCount=ans$Alignments)

	# spawn off the Align stats process
	if (chunkSize > 0) {
	    cat( verboseOutputDivider, "\nSpawning 'AlignStats' child process on Genomic reads...\n")
	    dispatch.AlignStats( sampleID=sample, annotationFile=annotationFile,
			optionsFile=optionsFile, mode="QuickQC", banner=banner, chunkSize=chunkSize,
			maxReads=maxReads, pause=pause, results.path=results.path, fastqFile=filein)
	    Sys.sleep(2)
	}

	# do a splice alignment next
	if ( dataType %in% c("RNA-seq","RIP-seq") && spliceIndex != "") {

		cat( verboseOutputDivider, "\nAlignment against Splice Index...\n")
		ans <- fastqToBAM( notGenomicFile, spliceFile, sampleID=sample, k=5,
				optionsFile=optionsFile, annotationFile=annotationFile, noHitsFile=nohitFile,
				alignIndex=spliceIndex, alignPolicy=splicePolicy, maxReads=maxReads, verbose=T)
		nSplice <-  ans$UniqueReads + ans$MultiReads
		nNoHit <- ans$NoHitReads
		quickFileLineCountRecord( nohitFile, sample, lineCount=(nNoHit*4), readCount=nNoHit)
		quickFileLineCountRecord( spliceFile, sample, lineCount=nSplice)
		notGenomicNeedsDelete <- TRUE
		# if no aligments, delete that BAM file
		if ( nSplice < 1) {
			file.delete( spliceFile)
		} else {

			# convert and spawn the stats
			if (verbose) cat( "\n\nConverting Splices back to genomic alignments...\n")
			ans <- spliceConvert( spliceFile, spliceFileOut, genomefile=genomicFile, sampleID=sampleID, 
					spliceMapPath=indexPath, spliceMapPrefix=spliceMapPrefix, rawReadCount=nRawReads, 
					readBufferSize=readBufferSize/10, verbose=verbose)
			quickFileLineCountRecord( spliceFileOut, sampleID, lineCount=ans$Alignments)
			cat("\n")
			cat( ans$textSummary)
			if (chunkSize > 0) {
		    	cat( verboseOutputDivider, "\nSpawning 'AlignStats' child process on Splice reads...\n")
		    	dispatch.AlignStats( sampleID=sample, annotationFile=annotationFile,
					optionsFile=optionsFile, mode="QuickQC.splicing", banner=banner, chunkSize=chunkSize,
					maxReads=maxReads, pause=pause, results.path=results.path, fastqFile=filein)
		    	Sys.sleep(2)
			}
			# do not delete the unconverted splices...
			#file.delete( spliceFile)
		}
	} else {
		# no splicing
		file.rename( notGenomicFile, nohitFile)
		notGenomicNeedsDelete <- FALSE
	}
	if (notGenomicNeedsDelete) file.delete( notGenomicFile)


	}  # end of 'mode!=pieOnly'

	if ( mode %in% c( "all", "pieOnly")) {
	    USRpath <- file.path( results.path, "USR")
	    USRfile <- file.path( USRpath, paste( sample, "USR.rda", sep="."))
	    if ( ! file.exists( USRfile)) {
		# look at the noHits
		cat( verboseOutputDivider, "\nMining No-Hits...\n")
		# build the Unique Short Reads data structure
		pipe.USR( sample, annotationFile=annotationFile, optionsFile=optionsFile,
			mode="QuickQC", resultsPath=results.path, nUSRkeep=nUSRkeep)
	    }
	}

	if ( mode == "all") {

	    # now we can look at consensus reads
	    blastIndex <- getOptionValue( optTbl, "blastIndex", notfound="nt")
	    dispatch.CR_Investigate( sample, annotationFile, optionsFile, mode="QuickQC", 
		blastIndex=blastIndex, doCR=TRUE, doBlast=TRUE, pause=pause,
		maxReads=500000, maxTime=300, maxCycles=8, nIterations=500, nBest=12, 
		ratePerCycle=1, results.path=results.path, verbose=TRUE) 
	    Sys.sleep(2)

	}  # end of 'mode==all'


	# final step, the Alignment Pie
	checkX11( width=8, height=6, xpos=10, ypos=10, bg="white")
	if ( is.null( nRawReads)) nRawReads <- maxReads

	cat( verboseOutputDivider, "\nPie Chart of 'Quick Alignment' results...\n")
	trueMaxReads <- min( maxReads, nRawReads, na.rm=T)
	ans <- pipe.AlignmentPie( sample, annotationFile=annotationFile, optionsFile=optionsFile,
				mode="QuickQC", results.path=results.path, fastqFile=filein,
				totalRawReads=trueMaxReads)

	cat( "\n\nQuick QC:  Final Reads Analysis:  ", sample,  "\n")
	print( ans)

	USR_cleanup()
	gc()

	}  # end of pairs...

	# spawn off a quick transcriptome
	if ( mode %in% c( "all", "alignOnly")) {
		cat( verboseOutputDivider, "\nWiggle Pileups and Transcriptome...\n")
		pipe.TranscriptPlusHTML( sampleID=sampleID, annotationFile=annotationFile,
				optionsFile=optionsFile, mode="QuickQC", 
				loadWIG=TRUE, banner=banner, 
				maxReads=maxReads, pause=pause, results.path=results.path)
	}

	# and maybe some peak calling...
	if ( dataType == "RIP-seq") {
		pipe.RIPpeaks( sampleID=sampleID, annotationFile=annotationFile,
				optionsFile=optionsFile, results.path=results.path)
	}

	return( ans)
}

