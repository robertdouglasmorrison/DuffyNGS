# pipe.Kmerize.R

# turn raw FASTQ reads into a 1-D table of Kmer counts


`pipe.Kmerize` <- function( sampleID=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt", 
				kmer.size=37, doCutadapt=TRUE, cutadaptProgram=Sys.which("cutadapt"), forceMatePairs=NULL, 
				buffer.size=1000000, maxReads=NULL, verbose=TRUE) {

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nStarting  'Kmerize Pipeline' on Sample:     ", sampleID,
			"\nStart Date/Time:   \t", date(), "\n")
	}

	results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)

	# file(s) to process comes from annotation and options...
	rawFastq <- getRawFastqFileNames( sampleID, annotationFile, optionsFile, verbose=FALSE)
	inputFastqFiles <- rawFastq$files
	asMatePairs <- rawFastq$asMatePairs

	# if we are doing the CutAdapt pass, check for those files, and call it if we need.
	if (doCutadapt) {
		# see if the files we need are already there
		cat( "\n")
		if( grepl( "\\.fastq", inputFastqFiles[1])) expectedFiles <- sub( "\\.fastq", ".trimmed.fastq", inputFastqFiles)
		if( grepl( "\\.fq", inputFastqFiles[1])) expectedFiles <- sub( "\\.fq", ".trimmed.fastq", inputFastqFiles)
		if ( all( file.exists( expectedFiles))) {
			cat( "\nUsing existing 'cutadapt' results..")
			filesToDo <- expectedFiles
		} else {
			# is the data paired end or not
			if (asMatePairs || (!is.null(forceMatePairs) && forceMatePairs)) {
				if (verbose) cat( "\nRunning 'CutAdapt' on mate pair files..")
				filesToDo <- cutadapt( file1=inputFastqFiles[1], file2=inputFastqFiles[2], 
                        				cutadaptProgram=cutadaptProgram, min.length=kmer.size)
			} else {
				filesToDo <- vector()
				for ( i in 1:length(inputFastqFiles)) {
					if (verbose) cat( "\nRunning 'CutAdapt' on single file: ", basename(inputFastqFiles[i]))
					filesToDo[i] <- cutadapt( file1=inputFastqFiles[i], file2=NULL, cutadaptProgram=cutadaptProgram, 
								min.length=kmer.size)
				}
			}
		}
	} else { 
		filesToDo <- inputFastqFiles
	}

	# set up to write out the results as they grow
	kmer.path <- file.path( results.path, "Kmers")
	if ( ! file.exists( kmer.path)) dir.create( kmer.path, recursive=T)
	outfile <- file.path( kmer.path, paste( sampleID, "Kmer", kmer.size, "Table.rda", sep="."))
	# removed any previous results
	file.delete( outfile)

	gc()
	startT <- proc.time()

	nReadsIn <- nDistinctKmers <- nTotalKmers <- 0
	bigKmerTable <- smlTable <- vector()
	nFiles <- length( filesToDo)
	kmer.size <- as.integer( kmer.size)

	cat( "\nN_Files to Kmerize: ", nFiles, "\n")
	for ( i in 1:nFiles) {
		ans <- kmerizeOneFastqFile( filesToDo[i], kmer.size=kmer.size, buffer.size=buffer.size, maxReads=maxReads)
		smlTable <- ans$Kmer.Table
		nReadsNow <- ans$nReadsIn
		if ( i == 1) {
			bigKmerTable <- smlTable
		} else {
			cat( "  merge..")
			bigKmerTable <- mergeTables( bigKmerTable, smlTable)
		}
		#cat( "  revCompClean..")
		#bigKmerTable <- kmerRevCompCleanup( bigKmerTable)

		cat( "  save..")
		save( bigKmerTable, file=outfile)

		nReadsIn <- nReadsIn + nReadsNow
		nDistinctKmers <- length( bigKmerTable)
		nTotalKmers <- sum( bigKmerTable)
		cat( "\nFile: ", i, basename(filesToDo[i]), "  N_Distinct: ", nDistinctKmers, "  N_Total: ", nTotalKmers, "\n")
		if ( ! is.null(maxReads) && nReadsIn >= maxReads) {
			cat( "\nReached 'maxReads' count. Stopping early..")
			break
		}
	}

	# now do the cleanup...
	cat( verboseOutputDivider)
	cat( "\n\nFinished 'Kmerize Pipeline' on Sample:     ", sampleID, "\n")
	cat( "\nTiming Stats: \n")
	myTime <- elapsedProcTime( startT, proc.time(), N=nTotalKmers, what="Kmer")
	print( myTime)

	out <- list( "nReadsIn"=nReadsIn, "nDistinctKmers"=nDistinctKmers, "nTotalKemrs"=nTotalKmers)
	return( out)
}


`pipe.KmerTable` <- function( sampleIDset, annotationFile="Annotation.txt", optionsFile="Options.txt", 
				kmer.size=37) {

	results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	kmer.path <- file.path( results.path, "Kmers")

	fileSet <- file.path( kmer.path, paste( sampleIDset, "Kmer", kmer.size, "Table.rda", sep="."))
	found <- file.exists( fileSet)
	if ( any( !found)) {
		cat( "\nSome Kmer files not found: ", fileSet[!found])
		fileSet <- fileSet[ found]
		sampleIDset <- sampleIDset[ found]
	}
	NS <- length(fileSet)
	if ( NS < 1) return(NULL)

	# we need to build up one giant vector of Kmers
	allKmers <- vector()
	cat( "\nFind union of all Kmers..\n")
	for ( f in fileSet) {
		cat( "  load: ", basename(f))
		load(f)
		myKmers <- names(bigKmerTable)
		cat( "  N=", length(myKmers))
		if ( ! length(allKmers)) {
			allKmers <- myKmers
		} else { 
			cat( "  combine..")
			allKmers <- union( allKmers, myKmers)
		}
	}
	cat( "  re-sort..")
	ord <- base::order( allKmers)
	allKmers <- allKmers[ord]
	cat( "  N_Kmers = ", NK <- length(allKmers))

	kmerTbl <- matrix( 0, nrow=NK, ncol=NS)
	colnames(kmerTbl) <- sampleIDset
	rownames(kmerTbl) <- allKmers

	cat( "\nReloading to fill table..\n")
	for ( i in 1:NS) {
		f <- fileSet[i]
		cat( "  load: ", basename(f))
		load(f)
		myKmers <- names(bigKmerTable)
		myCounts <- as.numeric( bigKmerTable)
		cat( "  lookup..")
		wh <- match( myKmers, allKmers)
		v <- rep.int(0,NK)
		cat( "  store..")
		v[ wh] <- myCounts
		kmerTbl[ , i] <- v
	}
	cat( "\nDone.\n")
	return( kmerTbl)
}


`kmerizeOneFastqFile` <- function( filein, kmer.size=37, buffer.size=1000000, maxReads=NULL) {

	# get ready to read the fastq file in chuncks...
	filein <- allowCompressedFileName( filein)
	if ( ! file.exists( filein)) stop( paste("Can't find input file: ", filein))
	conIn <- openCompressedFile( filein, open="r")
	on.exit( close( conIn))

	# 4 lines per read
	chunkSize <- buffer.size
	readChunkSize <- chunkSize * 4
	nread <- 0
	hasMore <- TRUE
	firstPass <- TRUE
	bigTable <- vector()

	repeat {
		if ( ! hasMore) break
		if ( ! is.null(maxReads) && nread >= maxReads) break
		cat( "\nBuffer Read..")
		chunk <- readLines( conIn, n=readChunkSize)
		if ( length( chunk) < 1) break
		if ( length( chunk) < readChunkSize) hasMore <- FALSE

		# in a FASTQ file, it takes 4 lines for each read;  1 is ID, 2nd is the seq, 4th is scores
		# all we want is the raw reads
		ptrs <- seq.int( 2, length(chunk), by=4)
		seqTxt <- chunk[ ptrs]
		nread <- nread + length( seqTxt)
		cat( "  N_reads: ", prettyNum( as.integer(nread), big.mark=","))

		smlTable <- kmerizeOneChunk( seqTxt, kmer.size=kmer.size)

		if ( ! length(bigTable)) {
			bigTable <- smlTable
		} else {
			cat( "  merge..")
			bigTable <- mergeTables( bigTable, smlTable)
		}
		cat( "  N_Kmer:", length( bigTable))
	}

	#cat( "  revCompClean..")
	#bigTable <- kmerRevCompCleanup( bigTable)

	return( list( "nReadsIn"=nread, "Kmer.Table"=bigTable))
}


`kmerizeOneChunk` <- function( seqs, kmer.size=37) {

	# make Kmers of every raw read in this chunk
	# set up to get all N-mers for all the sequences in this chunk
	SAPPLY <- base::sapply
	SUBSTR <- base::substr
	TABLE <- base::table

	sizeM1 <- kmer.size - 1

	# high chance of having duplicates, so only do each one once
	seqTbl <- TABLE( seqs)
	seqCounts <- as.numeric( seqTbl)
	seqs <- names( seqTbl)
	nSeq <- length( seqs)

	kmersOut <- vector()
	nOut <- 0

	# we can save a bit of time by doing all reads of the same length at one time
	seqLens <- base::nchar( seqs)
	lenFac <- factor( seqLens)
	cat( "  ReadLen=")
	tapply( 1:nSeq, lenFac, function(x) {
			
		# given all the pointers to seqs of the same length
		lenNow <- seqLens[ x[1]]
		cat( ".", lenNow, sep="")
		veryLastStart <- lenNow - sizeM1
		if ( veryLastStart < 1) return(NULL)
		fromSet <- 1:veryLastStart
		toSet <- fromSet + sizeM1
		nSubstrs <- length( fromSet)
		for ( j in x) {
			thisCount <- seqCounts[j]
			kmer <- SUBSTR( rep.int( seqs[j], nSubstrs), fromSet, toSet)
			# we only want to keep valid DNA chunks, so discard any 'N's
			isNNN <- grep( "N", kmer, fixed=T)
			if ( length(isNNN)) kmer <- kmer[ -isNNN]
			if ( ! length(kmer)) next
			# account for the multiplicity of this read
			if ( thisCount > 1) kner <- rep.int( kmer, thisCount)
			n <- length(kmer)
			now <- (nOut+1) : (nOut+n)
			kmersOut[now] <<- kmer
			nOut <<- nOut + n
		}
		return(NULL)
	})

	return( TABLE( kmersOut))
}

`kmerRevCompCleanup` <- function( kmerTable) {

	# the Kmers may be from both strands, and we only want to keep one form or each
	kmers <- names( kmerTable)
	cnts <- as.numeric( kmerTable)
	N <- length(kmerTable)

	# see what the RevComp over every Kmer is, and where it lies
	cat( "  make RevComps..")
	rcKmer <- vapply( kmers, FUN=myReverseComplement, FUN.VALUE="ACGT", USE.NAMES=F)
	cat( "  locate..")
	where <- match( rcKmer, kmers, nomatch=0)

	# we will make a new 1-D table that has just the first form of each Kmer
	# we can combine those that have both forms
	cat( "  join..")
	out <- rep.int( 0, N)
	hasRC <- which( where > 0)
	if ( length(hasRC)) {
		myLoc <- (1:N)[hasRC]
		rcLoc <- where[hasRC]
		firstLoc <- pmin( myLoc, rcLoc)
		out[firstLoc] <- cnts[myLoc] + cnts[rcLoc]
		names(out)[firstLoc] <- kmers[firstLoc]
	}

	# then those without their RevComp just stay as is
	noRC <- which( where == 0)
	if ( length(noRC)) {
		out[noRC] <- cnts[noRC]
		names(out)[noRC] <- kmers[noRC]
	}

	# lastly thown away the empty slots
	drops <- which( out == 0)
	if ( length( drops)) out <- out[ -drops]
	Nout <- length(out)
	cat( "  Pct Reduction: ", round( (N-Nout) * 100 / N), "%")

	# done
	return(out)
}

