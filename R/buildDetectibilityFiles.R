# buildDetectabilityFile.R

# align the whole genome against itself to build a detectability profile datafile
# Note:  Only original bowtie  (bowtie1) has the behavior we need to know detectibilty

`buildDetectabilityFiles` <- function( speciesID, genomicFastaFile, 
			bowtieIndexPath=bowtiePar("IndexPath"),
			selfGenomicIndex=bowtiePar("GenomicIndex"),
			otherGenomicIndex=NULL, seqIDset=NULL, 
			outPath=".", readSize=32, chunkSize=1000000, verbose=TRUE) {
	
	outfile <- otheroutfiles <- NULL

	# load that species
	setCurrentSpecies( speciesID=speciesID)
	outfilePrefix <- getCurrentSpeciesFilePrefix()
	
	# use the first seqID as a unique term in the filenames
	outfile <- paste( outfilePrefix, "selfUniqueRegions.txt", sep=".")
	file.delete( outfile)

	nOtherGenomes <- 0
	if ( ! is.null( otherGenomicIndex)) {
		nOtherGenomes <- length( otherGenomicIndex)
		# be smarter about grabbing the prefix
		#otherNames <- substr( basename( otherGenomicIndex), 1,2)
		otherNames <- sub( "\\..+", "", basename( otherGenomicIndex))
		otheroutfiles <- vector( length=nOtherGenomes)
		for( k in 1:nOtherGenomes) {
			otheroutfiles[k] <- paste( paste( outfilePrefix, "vs", otherNames[k], sep="."), 
					"detectableRegions.txt", sep=".")
			otherIndex <- file.path( bowtieIndexPath, otherGenomicIndex[k])
			if ( ! file.exists( base::paste( otherIndex, "1.ebwt", sep="."))) {
					stop( paste( "Bowtie index for other genome not found:  ", otherIndex))
			}
		}
	}
	cat( "\nN_Other Genomes: ", nOtherGenomes, "\t", otheroutfiles)


	# local function for one SeqID
	buildDetectableOneSeqID <- function( seqID) {

		# use the seqID as a unique term in the filenames
		myoutfile <- paste( outfilePrefix, seqID, "selfUniqueRegions.txt", sep=".")
		smallFasta <- paste( outfilePrefix, seqID, "tmp.detectable.fasta", sep=".")
		smallAnswer <- paste( outfilePrefix, seqID, "tmp.detectable.bowtie", sep=".")
		file.delete( c( myoutfile, smallFasta, smallAnswer))

		# write to connections as the answer grows
		conOut1 <- file( myoutfile, open="wt")
		writeLines( paste( "SEQ_ID", "START", "END", sep="\t"), con=conOut1)

		# do the "vs Other" pass at the same time
		if ( nOtherGenomes > 0) {
			otherfiles <- otherCons <- vector( mode="list", length=nOtherGenomes)
			for( k in 1:nOtherGenomes) {
				otherfiles[[k]] <- paste( paste( outfilePrefix, "vs", otherNames[k], sep="."), 
					seqID, "detectableRegions.txt", sep=".")
				file.delete( otherfiles[[k]])
				otherCons[[k]] <- file( otherfiles[[k]], open="wt")
				writeLines( paste( "SEQ_ID", "START", "END", sep="\t"), con=otherCons[[k]])
			}
		}

		# turn this sequence into a giant vector of bases
		curSeqDNA <- strsplit( as.character( getFastaSeqFromFilePath( filePath=genomicFastaFile, seqID=seqID)),
				split="")[[1]]
		nBase <- length( curSeqDNA)
		cat( "\nChromosome: ", seqID, "\tN_bases: ", formatC(nBase,format="d", big.mark=","), "\n")

		nReadsThisSeq <- 0
		if ( nBase > 0) {
		    for( ib in seq( 1, nBase, by=chunkSize)) {

			# build a small .fasta file of N-mers
			if (verbose) cat( "making reads..")
			file.delete( smallFasta)
			ans <- syntheticReads( seqID=seqID, seqDNA=curSeqDNA, type="fasta", readSize=readSize, 
					firstBase=ib, nReads=chunkSize, outfile=smallFasta)
			nReads <- ans$nReads
			nReadsThisSeq <- nReadsThisSeq + nReads

			# make the command line by hand, and call bowtie
			if (verbose) cat( "  selfAlign..")
			# we must be unique against itself
			options <- " -v 0 -m 1 -B 1 -f -p 4 --quiet "
			# recent version sends info to stderr, trap & discard
			trapStdErr <- " 2>/dev/null"

			# make sure the given index exists
			index <- file.path( bowtieIndexPath, selfGenomicIndex)
			if ( ! file.exists( base::paste( index, "1.ebwt", sep="."))) {
				stop( paste( "Bowtie index not found:  ", index))
			}
			cmdLine <- base::paste( bowtiePar( "Program"), options, index, smallFasta, smallAnswer, trapStdErr, sep="  ")
			file.delete( smallAnswer)
			callBowtie( cmdLine, verbose=FALSE)

			# turn those unique alignments into detectable regions
			ans <- selfAlignsToDetectableRegions( infile=smallAnswer, seqID=seqID, readSize=readSize)
			 
			# every time thru, update the final output file
			write.table( ans, file=conOut1, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
			flush( conOut1)

			# also go against a second index if given one...
			if ( nOtherGenomes > 0) {
			    for( k in 1:nOtherGenomes) {
				# make the command line by hand, and call bowtie
				if (verbose) cat( " ", otherNames[k], "otherAlign..")
				# find any number of matches, but just report one
				options <- " -v 1 -k 1 -B 1 -p 4 -f --quiet "
				index <- file.path( bowtieIndexPath, otherGenomicIndex[k])
				if ( ! file.exists( base::paste( index, "1.ebwt", sep="."))) {
					stop( paste( "Bowtie index not found:  ", index))
				}
				cmdLine <- base::paste( bowtiePar( "Program"), options, index, smallFasta, 
							smallAnswer, trapStdErr, sep="  ")
				file.delete( smallAnswer)
				callBowtie( cmdLine, verbose=FALSE)
				# turn those unique alignments into detectable regions
				ans <- otherAlignsToDetectableRegions( infile=smallAnswer, readSize=readSize)
				if ( nrow(ans) > 0) {
					write.table( ans, file=otherCons[[k]], sep="\t", quote=FALSE, row.names=FALSE, 
							col.names=FALSE)
					flush( otherCons[[k]])
				}
			    }
			}

			if (verbose) cat( " ", seqID, "Bases=", formatC(nReadsThisSeq, format="d", big.mark=","), "\n")
		    }
		}
		file.delete( smallFasta)
		file.delete( smallAnswer)
		close( conOut1)
		if ( nOtherGenomes > 0) for( k in 1:nOtherGenomes) close( otherCons[[k]])

		return( nReadsThisSeq)
	}


	combinePartialDetectables <- function( outFile, otherGenomeFiles) {

		cat( "\nWriting 'selfUnique' Detectability file: ", outFile)
		allSmallFiles <- paste( outfilePrefix, seqMap$SEQ_ID, "selfUniqueRegions.txt", sep=".")
		file.combine( infiles=allSmallFiles, outfile=outFile, check.headers=TRUE)
		# now delete the temporaries
		file.delete( allSmallFiles)

		if ( nOtherGenomes > 0) {
		    for ( k in 1:nOtherGenomes) {
			cat( "\nWriting 'otherDetectable' Detectability file: ", otherGenomeFiles[k])
			allSmallFiles <- paste( paste( outfilePrefix, "vs", otherNames[k], sep="."), 
						seqMap$SEQ_ID, "detectableRegions.txt", sep=".")
			file.combine( infiles=allSmallFiles, outfile=otherGenomeFiles[k], check.headers=TRUE)
			# now delete the temporaries
			file.delete( allSmallFiles)
		    }
		}

		return()
	}
	# end of local functions...


	# allow a subset of seqIDs...
	seqMap <- getCurrentSeqMap()
	seqSet <- seqMap$SEQ_ID
	seqSet <- multicore.seqID.order( seqSet)
	if ( ! is.null( seqIDset)) seqSet <- intersect( seqSet, seqIDset)

	# do each chromosome
	ans <- multicore.lapply( seqSet, buildDetectableOneSeqID)

	# now combine all the partials
	combinePartialDetectables( outfile, otheroutfiles)

	cat( "\n")
	out <- list( "SelfUniqueFile"=outfile, "OtherDetectableFile"=otheroutfiles)
	return( out)
}


`selfAlignsToDetectableRegions` <- function( infile="bowtie.out", seqID, readSize=32) {

	# quietly catch 'no file' to allow an empty bowtie result
	if ( ! file.exists( infile)) return( data.frame())

	tmp <- readAlignmentFile( filein=infile, type="bowtie", checkLineCount=FALSE, verbose=FALSE)
	
	# double check that all are the right seqID
	tmp <- subset.data.frame( tmp, SEQ_ID == seqID)

	# catch empty...
	nRow <- nrow( tmp)
	if ( nRow < 1) return( data.frame())

	# for self uniques, the seqID is part of the result, and the alignment locataion is what we want
	seqID <- tmp$SEQ_ID[1]
	starts <- tmp$START
	return( alignsToDetectableRegions( seqID, starts, readSize))
}


`otherAlignsToDetectableRegions` <- function( infile="bowtie.out", readSize=32) {

	# quietly catch 'no file' to allow an empty bowtie result
	if ( ! file.exists( infile)) return( data.frame())

	tmp <- readAlignmentFile( filein=infile, type="bowtie", checkLineCount=FALSE, verbose=FALSE)
	nRow <- nrow( tmp)
	
	# catch empty...
	if ( nRow < 1) return( data.frame())

	# for other species detection, the seqID and the alignment locataion are in the 'readID' term that
	# was in the fasta file given to Bowtie.   They are <seqid>::<location>
	seqID <- tmp$READ_ID[1]
	seqID <- sub( "::.+$", "", seqID)

	starts <- sub( "^.+::", "", tmp$READ_ID)
	starts <- as.numeric( starts)
	return( alignsToDetectableRegions( seqID, starts, readSize))
}


`alignsToDetectableRegions` <- function( seqID, starts, readSize=32) {

	nRow <- length( starts)
	if ( nRow < 1) return( data.frame())

	froms <- tos <- vector( length=nRow)
	nreg <- 0
	seqid <- seqID

	# this traversal assumes the starts are in assending order... Force IT!
	starts <- base::sort( starts)

	for ( i in 1:nRow) {

		thisStart <- starts[i]

		# special case for first entry...
		if ( i == 1) {
			firstStart <- prevStart <- thisStart
			next
		}

		# see if this read is adjacent to the previous, if so we are still in the same region
		if ( thisStart == (prevStart + 1)) {
			prevStart <- thisStart
			next
		}

		# when here, its a new region starting...
		
		# so write the previous
		nreg <- nreg + 1
		froms[nreg] <- firstStart
		tos[nreg] <- prevStart + readSize - 1

		# now start this one
		firstStart <- prevStart <- thisStart
	}

	# when done, write the last
	nreg <- nreg + 1
	froms[nreg] <- firstStart
	tos[nreg] <- prevStart + readSize - 1

	# trim to true size, and send it back
	length( froms) <- length( tos) <- nreg
	out <- data.frame( rep( seqid, times=nreg), as.integer(froms), as.integer(tos), stringsAsFactors=FALSE)
	colnames( out) <- c( "SEQ_ID", "START", "END")
	rownames( out) <- 1:nreg
	return( out)
}

