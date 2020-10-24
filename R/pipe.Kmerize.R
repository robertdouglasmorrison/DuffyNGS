# pipe.Kmerize.R

# turn raw FASTQ reads into a 1-D table of Kmer counts


`pipe.Kmerize` <- function( sampleID=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt", 
				kmer.size=33, doCutadapt=TRUE, cutadaptProgram=Sys.which("cutadapt"), forceMatePairs=NULL, 
				buffer.size=1000000, maxReads=NULL, min.count=2, verbose=TRUE) {

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nStarting  'Kmerize Pipeline' on Sample:     ", sampleID,
			"\nStart Date/Time:   \t", date(), 
			"\nKmer Size:         \t", kmer.size, "\n")
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
		if( grepl( "\\.fq", inputFastqFiles[1])) expectedFiles <- sub( "\\.fq", ".trimmed.fq", inputFastqFiles)
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
	bigKmerTable <- vector()
	nFiles <- length( filesToDo)
	kmer.size <- as.integer( kmer.size)

	cat( "\nN_Files to Kmerize: ", nFiles, "\n")
	for ( i in 1:nFiles) {
		ans <- kmerizeOneFastqFile( filesToDo[i], kmer.size=kmer.size, buffer.size=buffer.size, 
					maxReads=maxReads, min.count=min.count)
		#smlTable <- ans$Kmer.Table
		nReadsNow <- ans$nReadsIn
		if ( i == 1) {
			bigKmerTable <- ans$Kmer.Table
		} else {
			cat( "  merge..")
			bigKmerTable <- mergeTables( bigKmerTable, ans$Kmer.Table)
	
			cat( "  redo RevComp check..")
			bigKmerTable <- kmerRevComp.1Dtable( bigKmerTable)
		}

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
	rm( bigKmerTable, ans)
	gc()

	cat( verboseOutputDivider)
	cat( "\n\nFinished 'Kmerize Pipeline' on Sample:     ", sampleID, "\n")
	cat( "\nTiming Stats: \n")
	myTime <- elapsedProcTime( startT, proc.time(), N=nTotalKmers, what="Kmer")
	print( myTime)

	out <- list( "nReadsIn"=nReadsIn, "nDistinctKmers"=nDistinctKmers, "nTotalKemrs"=nTotalKmers)
	return( out)
}


`pipe.KmerTable` <- function( sampleIDset, annotationFile="Annotation.txt", optionsFile="Options.txt", 
				kmer.size=33, min.count=5) {

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
	cat( "\nRe-sorting.. ")
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
	cat( "\nDone loading.\n")
	
	cat( "\nChecking for low coverage Kmers to drop..")
	rowMax <- apply( kmerTbl, 1, max)
	drops <- which( rowMax < min.count)
	if ( length(drops)) {
		kmerTbl <- kmerTbl[ -drops, ]
		cat( "  Removed", length(drops), "Kmers for max count <", min.count)
	}
	
	return( kmerTbl)
}


`pipe.KmerCompare` <- function( kmerTbl, sampleIDset, groupSet) {

	# just just the samples asked for
	where <- match( sampleIDset, colnames(kmerTbl))
	if ( any( is.na(where))) stop( "Some SampleIDs not in Kmer table")
	kmerTbl <- kmerTbl[ , where]
	NC <- ncol(kmerTbl)
	NR <- nrow(kmerTbl)
	
	# get our group factors
	grpFac <- factor( groupSet)
	grpLvls <- levels( grpFac)
	nGrps <- nlevels( grpFac)
	if ( nGrps != 2) stop( "Expected exactly 2 Sample Groups")
	cat( "\nBreakdown by Group:\n")
	print( table( groupSet))
	
	# convert everything to normalized unit
	cat( "\nConverting to LKPB (Log2 Kmers Per Billion)..")
	kSums <- apply( kmerTbl, 2, sum)
	lkpbTbl <- kmerTbl
	for ( i in 1:NC) {
		myCnts <- kmerTbl[ ,i]
		myKPB <- myCnts * 1e9 / kSums[i]
		myLKPB <- round( log2( myKPB + 1), digits=3)
		lkpbTbl[ , i] <- myLKPB
	}
	
	# we are now ready to do tests per Kmer
	cat( "\nRunning Linear Models on ", NR, " Kmers..\n")
	est <- pval <- vector( length=NR)
	avg <- matrix( 0, 2, NR)
	rownames(avg) <- paste( "Avg", grpLvls, sep="_")
	for ( i in 1:NR) {
		v <- lkpbTbl[ i, ]
		ans <- lm( v ~ grpFac)
		ans2 <- summary(ans)
		cf <- coef(ans2)
		est[i] <- cf[ 2,1]
		pval[i] <- cf[ 2,4]
		avg[ ,i] <- tapply( v, grpFac, mean)
		if ( i %% 1000 == 0) cat( "\r", i, rownames(kmerTbl)[i], est[i], pval[i])
	}
	
	out <- data.frame( "Kmer"=rownames(kmerTbl), t(avg), "Difference"=est, "P.Value"=pval,
			stringsAsFactors=F)
	ord <- diffExpressRankOrder( out$Difference, out$P.Value)
	out <- out[ ord, ]
	rownames(out) <- 1:NR
	
	return(out)
}


`pipe.KmerAlignToGenome` <- function( kmers, optionsFile="Options.txt", quiet=TRUE) {

	# map some Kmers onto the reference genome
	N <- length( kmers)
	
	# make a temp FASTA file
	kmerIDs <- paste( "Kmer", 1:N, sep="_")
	kmerFastaFile <- "Temp.Kmers.BowtieInput.fasta"
	writeFasta( as.Fasta( kmerIDs, kmers), kmerFastaFile)
	kmerBowtieResultFile <- "Temp.Kmers.BowtieResult.bam"
	
	# build a call to Bowtie2
	optT <- readOptionsTable( optionsFile)
	# force the look at options file to initialize bowtie.
	bowtie2Par.defaults( optionsFile, verbose=FALSE)
	## build the Bowtie alignment program's command line
	progName <- getOptionValue( optT, "bowtie2Program", verbose=FALSE)
	cmdline <- progName
	# next is "input options" using FASTA instead of default FASTQ
	cmdLine <- mypaste( cmdLine, "-a")
	# next is explicit alignment policy
	cmdLine <- mypaste( cmdLine, " --very-sensitive ")
	if ( quiet) cmdLine <- mypaste( cmdLine, " --quiet")
	# when not catching the unaligned, allow explicit throw-away of unaligned
	cmdLine <- mypaste( cmdLine, " --no-unal ")
	# next is threads
	nCores <- as.integer( getOptionValue( optT, "nCores", notfound="4", verbose=FALSE))
	cmdLine <- mypaste( cmdLine, " --threads", nCores)
	# the index
	index.path <- getOptionValue( optT, "bowtie2Index.path", notfound=".", verbose=FALSE)
	alignIndex <- getOptionValue( optT, "GenomicIndex", verbose=F)
	myIndexFile <- file.path( index.path, alignIndex)
	# test to see that a file is really there...
	if ( ! file.exists(  mypaste( myIndexFile, ".1.bt2", sep=""))) {
	    # new long indexex are possible
	    if ( ! file.exists(  mypaste( myIndexFile, ".1.bt2l", sep=""))) {
		stop( paste( "Bowtie2 index Path and/or File not found. \n", "Filename as given:  ", myIndexFile))
	    }
	}
	cmdLine <- mypaste( cmdLine, " -x", myIndexFile)
	# the input file
	cmdLine <- mypaste( cmdLine, " -U", kmerFastaFile)
	# lastly, the output destination
	outputTerm <- paste ( " | samtools view -bS -o", kmerBowtieResultFile, " - ")
	cmdLine <- mypaste( cmdLine, outputTerm, sep="  ")
	callBowtie2( cmdLine, verbose=verbose)
}


`kmerizeOneFastqFile` <- function( filein, kmer.size=33, buffer.size=1000000, maxReads=NULL,
				min.count=2) {

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
		seqTxt <- chunk[ seq.int( 2, length(chunk), by=4)]
		nread <- nread + length( seqTxt)
		cat( "  N_reads: ", prettyNum( as.integer(nread), big.mark=","))
		rm( chunk)

		smlTable <- kmerizeOneChunk( seqTxt, kmer.size=kmer.size)

		if ( ! length(bigTable)) {
			bigTable <- smlTable
		} else {
			cat( "  merge..")
			bigTable <- mergeTables( bigTable, smlTable)
		}
		cat( "  N_Kmer:", length( bigTable))
	}
	rm( smlTable)

	# any K-mers with too few observations are not to be kept
	tooFew <- which( bigTable < min.count)
	if ( length( tooFew)) {
		cat( "\n  Drop if count < ", min.count, " (N=", length(tooFew), ")", sep="")
		bigTable <- bigTable[ -tooFew]
		cat( "  N_Kmer:", length( bigTable))
	}
	
	# also do any RevComp resolving now too
	cat( "\nSearch for RevComp pairs to join..")
	bigTable <- kmerRevComp.1Dtable( bigTable)
	cat( "  N_Kmer:", length( bigTable))

	return( list( "nReadsIn"=nread, "Kmer.Table"=bigTable))
}


`kmerizeOneChunk` <- function( seqs, kmer.size=33) {

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


`kmerRevComp.1Dtable` <- function( kmerTable) {

	# the Kmers may be from both strands, and we only want to keep one form of each
	# so merge those where both are present
	
	# allow being given a saved file
	isFile <- FALSE
	if ( !is.table( kmerTable) && length(kmerTable) == 1 && is.character(kmerTable) && grepl( "\\.rda$", kmerTable)) {
		kmerFile <- kmerTable
		load( kmerTable)
		kmerTable <- bigKmerTable
		isFile <- TRUE
	}
	
	kmers <- names( kmerTable)
	cnts <- as.numeric( kmerTable)
	N <- length(kmerTable)
	cat( "\nN_Kmers in: ", N)

	# see what the RevComp over every Kmer is, and where it lies
	rcKmer <- findKmerRevComp( kmers)
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
	cat( "\nN_Kmers out: ", Nout)
	
	cat( "  Pct Reduction: ", round( (N-Nout) * 100 / N), "%")

	# done
	if (isFile) {
		bigKmerTable <- out
		cat( "  saving..")
		save( bigKmerTable, file=kmerFile)
		return( kmerFile)
	} else {
		return( out)
	}
}


findKmerRevComp <- function( kmers) {

	N <- length( kmers)
	out <- rep.int( "", N)
	
	# allow 'fast' lookup via a file of previously done kmers
	kmerXrefFile <- "Kmer.RevComp.Xref.rda"
	if ( file.exists( kmerXrefFile)) {
		cat( "  load RevComp Xref..")
		load( kmerXrefFile)
	} else {
		xref <- data.frame()
	}
	nXref <- nrow(xref)

	# see if any of these are known
	if ( nrow(xref)) {
		cat( "  lookup..")
		where1 <- match( kmers, xref$Kmer, nomatch=0)
		out[ where1 > 0] <- xref$RevComp[ where1]
		where2 <- match( kmers, xref$RevComp, nomatch=0)
		out[ where2 > 0] <- xref$Kmer[ where2]
		cat( "  N_Found=", sum(where1 > 0 | where2 > 0))
	}
	
	# calculate the rest
	toCalc <- which( out == "")
	if ( length( toCalc)) {
		cat( "  calc RevComps.. (N=", length(toCalc), ")", sep="")
		rcKmer <- vapply( kmers[ toCalc], FUN=myReverseComplement, FUN.VALUE="ACGT", USE.NAMES=F)
		out[ toCalc] <- rcKmer
	
		# store these for the future
		smlXref <- data.frame( "Kmer"=kmers[toCalc], "RevComp"=rcKmer, stringsAsFactors=FALSE)
		xref <- rbind( xref, smlXref)
		save( xref, file=kmerXrefFile)
	}
	return( out)
}