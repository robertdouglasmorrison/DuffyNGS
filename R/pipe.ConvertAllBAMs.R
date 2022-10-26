# pipe.ConvertAllBAMs.R

# do any modifications to BAM files after Bowtie created them

`pipe.ConvertAllBAMs` <- function( sampleID,  annotationFile="Annotation.txt", optionsFile="Options.txt",
			dataType=c("RNA-seq", "DNA-seq", "ChIP-seq", "RIP-seq"), rawReadCount=NULL, verbose=TRUE) {

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nConverting BAM files for Sample:     ", sampleID, "\n\n")
	}
	gc()
	dataType <- match.arg( dataType)

	# convert/append any BAM file details
	optT <- readOptionsTable( optionsFile)
	results.path <- getOptionValue( optT, "results.path", notfound=".")
	riboFile <- file.path( results.path, "riboClear", paste( sampleID, "ribo.bam", sep="."))
	genomicFile <- file.path( results.path, "align", paste( sampleID, "genomic.bam", sep="."))
	spliceFile <- file.path( results.path, "splicing", paste( sampleID, "splice.bam", sep="."))

	nSpawn <- 0
	watchRibo <- FALSE
	watchRiboFile <- paste( sampleID, "convertingRiboBAM.done", sep=".")
	riboLogFile <- paste( sampleID, "convertRiboBAM.log.txt", sep=".")
	file.delete( watchRiboFile)
	if ( dataType == "RNA-seq" && file.exists( riboFile)) {
		watchRibo <- TRUE
		nSpawn <- nSpawn + 1
		dispatch.ConvertRiboBAM( sampleID, annotationFile, optionsFile, dataType, rawReadCount)
		Sys.sleep(2)
	}

	watchSplice <- FALSE
	watchSpliceFile <- paste( sampleID, "convertingSpliceBAM.done", sep=".")
	spliceLogFile <- paste( sampleID, "convertSpliceBAM.log.txt", sep=".")
	file.delete( watchSpliceFile)
	if ( dataType %in% c( "RNA-seq", "RIP-seq") && file.exists( spliceFile)) {
		watchSplice <- TRUE
		nSpawn <- nSpawn + 1
		dispatch.ConvertSpliceBAM( sampleID, annotationFile, optionsFile, dataType, rawReadCount)
		Sys.sleep(2)
	}

	watchGenomic <- FALSE
	watchGenomicFile <- paste( sampleID, "convertingGenomeBAM.done", sep=".")
	genomicLogFile <- paste( sampleID, "convertGenomicBAM.log.txt", sep=".")
	file.delete( watchGenomicFile)
	if ( file.exists( genomicFile)) {
		watchGenomic <- TRUE
		nSpawn <- nSpawn + 1
		dispatch.ConvertGenomicBAM( sampleID, annotationFile, optionsFile, dataType, rawReadCount)
		Sys.sleep(2)
	}


	if ( any( c( watchRibo, watchGenomic, watchSplice))) {
		cat( "\n\nWaiting for", nSpawn, "spawned jobs to complete")
		timeOutCount <- 0
		repeat {
			Sys.sleep(60)
			cat( ".")
			if ( watchRibo) {
				if ( file.exists(watchRiboFile)) {
					cat( "\n\nRibo Conversion done.\n")
					watchRibo <- FALSE
					file.delete( watchRiboFile)
					ans <- readLines( riboLogFile)
					writeLines( ans)
					file.delete( riboLogFile)
					cat( "\n")
				}
			}
			if ( watchSplice) {
				if ( file.exists(watchSpliceFile)) {
					cat( "\n\nSplice Conversion done.\n")
					watchSplice <- FALSE
					file.delete( watchSpliceFile)
					ans <- readLines( spliceLogFile)
					writeLines( ans)
					file.delete( spliceLogFile)
					cat( "\n")
				}
			}
			if ( watchGenomic) {
				if ( file.exists(watchGenomicFile)) {
					cat( "\n\nGenomic Conversion done.\n")
					watchGenomic <- FALSE
					file.delete( watchGenomicFile)
					ans <- readLines( genomicLogFile)
					writeLines( ans)
					file.delete( genomicLogFile)
					cat( "\n")
				}
			}
			if ( ! any( c( watchRibo, watchGenomic, watchSplice))) break

			# there is a tiny chance that one could die...  bale once genomic is done
			timeOutCount <- timeOutCount + 1
			if ( timeOutCount > 360 && !watchGenomic) {
				cat( "\nGenome done converting..  Timed out on some others..  Continuing on..")
				break
			}
		} 
	}

	cat( "\n\nFinished 'ConvertAllBAMs' on Sample:     ", sampleID, "\n\n")
	gc()
	return()
}


`dispatch.ConvertRiboBAM` <- function( sampleID,  annotationFile="Annotation.txt", optionsFile="Options.txt",
			dataType=c("RNA-seq", "DNA-seq", "ChIP-seq", "RIP-seq"), rawReadCount=NULL, verbose=TRUE) {

	dataType <- match.arg( dataType)

	commandLine <- paste( " pipe.ConvertRiboBAM( sampleID=\"", 
				sampleID, "\", annotationFile=\"", annotationFile, 
				"\", optionsFile=\"", optionsFile, 
				"\", dataType=\"", dataType, 
				"\", rawReadCount=", if (is.null(rawReadCount)) "NULL" else as.integer(rawReadCount), 
				" )", sep="")

	logFile=paste( sampleID, "convertRiboBAM.log.txt", sep=".")
	spawnRcmd( commandLine, Rpackage="DuffyNGS", logFile=logFile, verbose=FALSE)
	return()
}


`dispatch.ConvertGenomicBAM` <- function( sampleID,  annotationFile="Annotation.txt", optionsFile="Options.txt",
			dataType=c("RNA-seq", "DNA-seq", "ChIP-seq", "RIP-seq"), rawReadCount=NULL, verbose=TRUE) {

	dataType <- match.arg( dataType)

	commandLine <- paste( " pipe.ConvertGenomicBAM( sampleID=\"", 
				sampleID, "\", annotationFile=\"", annotationFile, 
				"\", optionsFile=\"", optionsFile, 
				"\", dataType=\"", dataType, 
				"\", rawReadCount=", if (is.null(rawReadCount)) "NULL" else as.integer(rawReadCount), 
				" )", sep="")

	logFile=paste( sampleID, "convertGenomicBAM.log.txt", sep=".")
	spawnRcmd( commandLine, Rpackage="DuffyNGS", logFile=logFile, verbose=FALSE)
	return()
}


`dispatch.ConvertSpliceBAM` <- function( sampleID,  annotationFile="Annotation.txt", optionsFile="Options.txt",
			dataType=c("RNA-seq", "DNA-seq", "ChIP-seq", "RIP-seq"), rawReadCount=NULL, verbose=TRUE) {

	dataType <- match.arg( dataType)

	commandLine <- paste( " pipe.ConvertSpliceBAM( sampleID=\"", 
				sampleID, "\", annotationFile=\"", annotationFile, 
				"\", optionsFile=\"", optionsFile, 
				"\", dataType=\"", dataType, 
				"\", rawReadCount=", if (is.null(rawReadCount)) "NULL" else as.integer(rawReadCount), 
				" )", sep="")

	logFile=paste( sampleID, "convertSpliceBAM.log.txt", sep=".")
	spawnRcmd( commandLine, Rpackage="DuffyNGS", logFile=logFile, verbose=FALSE)
	return()
}


`pipe.ConvertRiboBAM` <- function( sampleID,  annotationFile="Annotation.txt", optionsFile="Options.txt",
			dataType=c("RNA-seq", "DNA-seq", "ChIP-seq", "RIP-seq"), rawReadCount=NULL, verbose=TRUE) {

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nConverting Ribo BAM file for Sample:     ", sampleID, "\n\n")
	}
	gc()
	dataType <- match.arg( dataType)

	# convert/append any BAM file details
	optT <- readOptionsTable( optionsFile)
	target <- getAndSetTarget( optionsFile, sampleID=sampleID, annotationFile=annotationFile)
	results.path <- getOptionValue( optT, "results.path", notfound=".")
	bowtie2Par.defaults( optionsFile)
	readBufferSize <- as.integer( getOptionValue( optT, "readBufferSize", notfound="1000000"))
	genomicFile <- file.path( results.path, "align", paste( sampleID, "genomic.bam", sep="."))
	riboFile <- file.path( results.path, "riboClear", paste( sampleID, "ribo.bam", sep="."))
	riboFileOut <- file.path( results.path, "riboClear", paste( sampleID, "ribo.converted.bam", sep="."))
	indexPath <- getOptionValue( optT, "bowtie2Index.path")
	readSense <- getReadSense( sampleID, annotationFile)
	riboMap <- getOptionValue( optT, "RiboMap")
	if ( nchar( riboMap) < 3) {
		cat( "\nError:  Required option 'RiboMap' has no valid filename")
		stop( "Failed to convert ribo alignments")
	}
	riboMap <- file.path( indexPath, riboMap)

	if ( dataType %in% c( "RNA-seq", "RIP-seq") && file.exists( riboFile)) {
		if (verbose) cat( "\n\nConverting riboClear hits back to genomic positions...\n")

		ans <- riboConvert( riboFile, riboFileOut, genomefile=genomicFile, sampleID=sampleID, 
				riboMapFile=riboMap, rawReadCount=rawReadCount, 
				readBufferSize=readBufferSize, verbose=verbose)

		quickFileLineCountRecord( riboFileOut, sampleID, lineCount=ans$Alignments)

		sum.path <- file.path( results.path, "summary")
		if ( ! file.exists( sum.path)) dir.create( sum.path, recursive=T)
		sumfile <- file.path( sum.path, paste( sampleID, "ribo.Summary.txt", sep="."))
		con <- file( sumfile, open="at")
		writeLines( ans$textSummary, con=con, sep="")
		close( con)
	}

	# set the 'we are done' flag
	watchFile <- paste( sampleID, "convertingRiboBAM.done", sep=".")
	writeLines( "Done", watchFile)
	return()
}


`pipe.ConvertGenomicBAM` <- function( sampleID,  annotationFile="Annotation.txt", optionsFile="Options.txt",
			dataType=c("RNA-seq", "DNA-seq", "ChIP-seq", "RIP-seq"), rawReadCount=NULL, verbose=TRUE) {

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nConverting Genomic BAM file for Sample:     ", sampleID, "\n\n")
	}
	gc()
	dataType <- match.arg( dataType)

	# convert/append any BAM file details
	optT <- readOptionsTable( optionsFile)
	target <- getAndSetTarget( optionsFile, sampleID=sampleID, annotationFile=annotationFile)
	results.path <- getOptionValue( optT, "results.path", notfound=".")
	bowtie2Par.defaults( optionsFile)
	readBufferSize <- as.integer( getOptionValue( optT, "readBufferSize", notfound="1000000"))
	genomicFile <- file.path( results.path, "align", paste( sampleID, "genomic.bam", sep="."))

	if ( file.exists( genomicFile)) {
		if (verbose) cat( "\n\nAppending GeneID terms to genomic alignments...\n")
		ans <- genomicConvert( genomicFile, sampleID=sampleID, readBufferSize=readBufferSize, 
				rawReadCount=rawReadCount, verbose=verbose)

		quickFileLineCountRecord( genomicFile, sampleID, lineCount=ans$Alignments)

		sum.path <- file.path( results.path, "summary")
		if ( ! file.exists( sum.path)) dir.create( sum.path, recursive=T)
		sumfile <- file.path( sum.path, paste( sampleID, "genomic.Summary.txt", sep="."))
		con <- file( sumfile, open="at")
		writeLines( ans$textSummary, con=con, sep="")
		close( con)
	}

	# set the 'we are done' flag
	watchFile <- paste( sampleID, "convertingGenomeBAM.done", sep=".")
	writeLines( "Done", watchFile)
	return()
}


`pipe.ConvertSpliceBAM` <- function( sampleID,  annotationFile="Annotation.txt", optionsFile="Options.txt",
			dataType=c("RNA-seq", "DNA-seq", "ChIP-seq", "RIP-seq"), rawReadCount=NULL, verbose=TRUE) {

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nConverting Splice BAM file for Sample:     ", sampleID, "\n\n")
	}
	gc()
	dataType <- match.arg( dataType)

	# convert/append any BAM file details
	optT <- readOptionsTable( optionsFile)
	target <- getAndSetTarget( optionsFile, sampleID=sampleID, annotationFile=annotationFile)
	results.path <- getOptionValue( optT, "results.path", notfound=".")
	bowtie2Par.defaults( optionsFile)
	readBufferSize <- as.integer( getOptionValue( optT, "readBufferSize", notfound="1000000"))
	readBufferSize <- readBufferSize / 10
	#if ( readBufferSize > 200000) readBufferSize <- 200000
	genomicFile <- file.path( results.path, "align", paste( sampleID, "genomic.bam", sep="."))
	spliceFile <- file.path( results.path, "splicing", paste( sampleID, "splice.bam", sep="."))
	spliceFileOut <- file.path( results.path, "splicing", paste( sampleID, "splice.converted.bam", sep="."))
	indexPath <- getOptionValue( optT, "bowtie2Index.path")
	readSense <- getReadSense( sampleID, annotationFile)
	spliceMapPrefix <- getOptionValue( optT, "SpliceMapPrefix", notfound="spliceMap")

	if ( dataType %in% c( "RNA-seq", "RIP-seq") && file.exists( spliceFile)) {
		if (verbose) cat( "\n\nConverting splices back to genomic fragments...\n")

		ans <- spliceConvert( spliceFile, spliceFileOut, genomefile=genomicFile, 
					spliceMapPath=indexPath, spliceMapPrefix=spliceMapPrefix, 
					readSense=readSense, readBufferSize=readBufferSize, sampleID=sampleID,
					rawReadCount=rawReadCount)

		quickFileLineCountRecord( spliceFileOut, sampleID, lineCount=ans$Alignments)

		sum.path <- file.path( results.path, "summary")
		if ( ! file.exists( sum.path)) dir.create( sum.path, recursive=T)
		sumfile <- file.path( sum.path, paste( sampleID, "splice.Summary.txt", sep="."))
		con <- file( sumfile, open="at")
		writeLines( ans$textSummary, con=con, sep="")
		close( con)

		# turn the splice details into a file for each species
		tbl <- ans$spliceDetailSummary
		if ( !is.null(tbl) && !is.na(tbl) && nrow( tbl)) {
		    speciesSet <- sort( unique.default( tbl$SPECIES_ID))
		    speciesIDcolumn <- match( "SPECIES_ID", colnames(tbl), nomatch=0)
		    for ( i in 1:length( speciesSet)) {
			thisSpecies <- speciesSet[i]
			if ( is.na( thisSpecies) || thisSpecies == "") next
			thisPrefix <- getCurrentTargetFilePrefix()[ match(thisSpecies,getCurrentTargetSpecies())]
			sml <- subset( tbl, SPECIES_ID == thisSpecies)
			if (nrow(sml) < 1) next

			# prep the subset of splices in this one organism
			if ( speciesIDcolumn > 0) sml <- sml[ , -speciesIDcolumn]
			ttlReads <- sum( sml$N_READS, na.rm=T)
			pctReads <- sml$N_READS * 100 / ttlReads
			sml$PCT_READS <- format( pctReads, digits=3)
			sml$KEY_ID <- paste( sml$GENE_ID, sml$SPLICE_ID, sep="::")

			outfile <- paste( sampleID, "splice", thisPrefix, "Summary.txt", sep=".")
			outfile <- file.path( results.path, "splicing", outfile)
			write.table( sml, outfile, sep="\t", quote=FALSE, row.names=FALSE)
			cat( "\nWrote Splice Summary file: ", basename(outfile), "\tN_SpliceJunctions: ", nrow(sml))
		    }
		}
	}

	# set the 'we are done' flag
	watchFile <- paste( sampleID, "convertingSpliceBAM.done", sep=".")
	writeLines( "Done", watchFile)
	return()
}
