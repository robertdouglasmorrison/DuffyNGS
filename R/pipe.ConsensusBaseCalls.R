# pipe.ConsensusBaseCalls.R -- tools for pulling sequence infrom from BAM file results


`pipe.ConsensusBaseCalls` <- function( sampleID, geneID=NULL, seqID=NULL, start=NULL, stop=NULL, 
				annotationFile="Annotation.txt", optionsFile="Options.txt", results.path=NULL,
				aaToo=TRUE, noReadCalls=c("blank","genomic"), as.cDNA=FALSE, utr.tail.length=0,
				SNP.only=FALSE, minReadCalls=NULL, minPercentSNP=NULL, verbose=TRUE) {
				
	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	genomicFastaFile <- getOptionValue( optT, "genomicFastaFile", verbose=F)

	bamfile <- file.path( results.path, "align", paste( sampleID, "genomic.bam", sep="."))

	# control what to do when no reads cover a region
	noReadCalls <- match.arg( noReadCalls)

	if ( utr.tail.length != 0) {
		if (verbose) cat( "\nAdding UTR bases..  Turning off AA and cDNA step..")
		aaToo <- FALSE
		as.cDNA <- FALSE
	}

	# we have what we need, call the lower level tool
	ans <- consensusBaseCalls( bamfile, genomicFastaFile, geneID=geneID, seqID=seqID, 
			start=start, stop=stop, aaToo=aaToo, noReadCalls=noReadCalls, utr.tail.length=utr.tail.length, 
			SNP.only=SNP.only, minReadCalls=minReadCalls, minPercentSNP=minPercentSNP, verbose=verbose)
	if ( is.null(ans)) return(NULL)

	# translate genomic info back to cDNA units if we want
	if ( ! is.null(geneID) && as.cDNA) {
		ans <- consensusBaseCallsToCDNA( ans, geneID=geneID, start=start, stop=stop, verbose=verbose)
	} else {
		ans$callsTable <- NULL
	}
	return( ans)
}


`consensusBaseCalls` <- function( bamfile, genomicFastaFile, geneID=NULL, seqID=NULL, 
				start=NULL, stop=NULL, aaToo=TRUE, noReadCalls=c("blank","genomic"),
				utr.tail.length=0, SNP.only=FALSE, minReadCalls=NULL, minPercentSNP=NULL,
				verbose=TRUE) {
				
	GENOME_BASE <- ","

	# let's be a bit more flexible with what we pass in
	# Could be a genomic range {seqID, start, stop}, a gene {geneID}, or any generic sequence...
	isRange <- isGene <- isGeneric <- FALSE
	geneMap <- getCurrentGeneMap()
	cdsMap <- getCurrentCdsMap()
	seqMap <- getCurrentSeqMap()

	# specify by either gene or range
	if ( is.null( geneID)) {
		if ( is.null(seqID)) stop( "One of 'seqID','geneID' must be not NULL")
		if ( seqID %in% seqMap$SEQ_ID) {
			isRange <- TRUE
			if ( is.null(start)) start <- 1
			if ( is.null(stop)) stop <- subset( seqMap, SEQ_ID == seqID)$LENGTH[1]
			gmap <- subset.data.frame( geneMap, SEQ_ID == seqID & stop >= POSITION & start <= END)
			cmap <- subset.data.frame( cdsMap, GENE_ID %in% gmap$GENE_ID)
			geneID <- gmap$GENE_ID[1]
			geneID <- shortGeneName( geneID, keep=1)
		} else {
			isGeneric <- TRUE
			if ( any( c( is.null(start), is.null(stop)))) stop( "Generic consensus sequence needs explicit 'start,stop'")
		}
	} else {
		isGene <- TRUE
		gmap <- subset.data.frame( geneMap, GENE_ID == geneID)
		if ( ! nrow(gmap)) {
			cat( "\nGiven GeneID not found in current species: ", geneID, "\n")
			stop( "Check GeneID and/or current species")
		}
		cmap <- subset.data.frame( cdsMap, GENE_ID == geneID)
		if ( ! nrow(cmap)) {
			if (verbose) cat( "\nGeneID not found in CDS map.. Using Exon map instead.")
			cmap <- subset.data.frame( getCurrentExonMap(), GENE_ID %in% gmap$GENE_ID)
		}
		if (is.null(start)) start <- gmap$POSITION[1]
		if (is.null(stop)) stop <- gmap$END[1]
		seqID <- gmap$SEQ_ID[1]
		geneStrand <- gmap$STRAND[1]
		if (utr.tail.length != 0) {
			start <- start - utr.tail.length
			if ( start < 1) start <- 1
			stop <- stop + utr.tail.length
			chromoLen <- subset( seqMap, SEQ_ID == seqID)$LENGTH[1]
			if( stop > chromoLen) stop <- chromoLen
		}
	}

	# load that portion of the PILEUPS
	curMPU <- BAM.mpileup( bamfile, seqID=seqID, fastaFile=genomicFastaFile, start=start, stop=stop, 
			summarize.calls=TRUE, verbose=FALSE)
	#if ( is.null(curMPU)) return(NULL)
	if ( is.null(curMPU)) curMPU <- data.frame()

	# add the reference genome
	genomicStr <- getFastaSeqFromFilePath( genomicFastaFile, seqID)
	curGenomeDNA <- strsplit( as.character(genomicStr), split="")[[1]]

	# decide how we will treat missing data
	noReadCalls <- match.arg( noReadCalls)

	# gather the reference genome bases in this range
	allBases <- start : stop
	genomeBaseText <- curGenomeDNA[ allBases]
	names(genomeBaseText) <- allBases
	if (noReadCalls == "blank") {
		genomeSNPtext <- base::rep.int( "", length(genomeBaseText))
		names(genomeSNPtext) <- allBases
	} else {
		genomeSNPtext <- genomeBaseText
		names(genomeSNPtext) <- allBases
	}

	# make a default flips matrix that is all reference genome
	oneVec <- rep.int( 1, length(allBases))
	zeroVec <- rep.int( 0, length(allBases))
	genomicNonFlips <- data.frame( "Ref"=oneVec, "A"=zeroVec, "C"=zeroVec, "G"=zeroVec, "T"=zeroVec,
					"N"=zeroVec, "Indel"=zeroVec, stringsAsFactors=F)
	colnames(genomicNonFlips)[1] <- GENOME_BASE
	rownames(genomicNonFlips) <- allBases

	hasBaseCalls <- ( nrow( curMPU) > 0)
	if ( ! hasBaseCalls) {
		genomeAminoText <- ""
		if (aaToo && isGene && noReadCalls == "genomic") {
			ans <- convertGenomicBasesToCodingAA( seqID, position=start, end=stop, 
						strand=geneStrand, dnaQuery=genomeBaseText, genomeDNA=curGenomeDNA,
						geneMap=gmap, cdsMap=cmap)
			genomeAminoText <- ans$query
			names( genomeAminoText) <- allBases
		}
		out <- list( "ref"=genomeBaseText, "callsMatrix"=genomicNonFlips, "dna.consensus"=genomeSNPtext, 
				"aa.consensus"=genomeAminoText, "indel.details"="")
		return( out)
	}

	xLocs <- curMPU$POSITION

	# turn the base calls into an explicit matrix in order: genomic,A,C,G,T,N, indel
	flips <- MPU.callStringsToMatrix( curMPU$BASE_TABLE)
	if ( nrow(flips)) rownames(flips) <- xLocs

	# get the majority base at each SNP spot, if not a SNP it will be ',' (matches the reference
	WHICH <- base::which
	WHICH.MAX <- base::which.max
	MATCH <- base::match

	snpTopBase <- apply( flips, MARGIN=1, function(x) colnames(flips)[ WHICH.MAX(x)])
	names( snpTopBase) <- xLocs

	# allow a special mode when we are doing "use genomic for missing" mode
	# if given a minimum read depth, then any spots too shallow get called as reference
	# i.e.  not enough read depth to trust the base calls
	if (noReadCalls == "genomic" && ! is.null(minReadCalls)) {
		minReadCalls <- as.numeric( minReadCalls)
		depth <- apply( flips, MARGIN=1, sum, na.rm=T)
		tooLow <- which( depth < minReadCalls)
		if ( length(tooLow)) snpTopBase[ tooLow] <- GENOME_BASE
	}

	# if given a minimum percent SNP, then any spots too 'mixed' get called reference
	# i.e.  not enough consensus to trust the base calls
	if (noReadCalls == "genomic" && ! is.null(minPercentSNP)) {
		minPercentSNP <- as.numeric( minPercentSNP)
		big <- apply( flips, MARGIN=1, max, na.rm=T)
		depth <- apply( flips, MARGIN=1, sum, na.rm=T)
		pct <- big / depth
		tooLow <- which( pct < minPercentSNP)
		if ( length(tooLow)) snpTopBase[ tooLow] <- GENOME_BASE
	}

	# the Indels are tougher, we need to get the actual bases for both the reference and the indels and figure if
	# its an insertion or deletion
	indelDetails <- MPU.getIndelDetails( flips, curMPU)
	if ( indelDetails$nIndels > 0) {
		isIndel <- indelDetails$who
		snpTopBase[isIndel] <- indelDetails$bases
	}
	indelText <- indelDetails$indelText
	names( indelText) <- xLocs

	# there may be places with no read depth at all
	# default is to use genomic... or delete those bases
	# our methods below start from the genome calls, and then undo or change them
	if ( noReadCalls == "blank") {
		genomeSNPtext <- genomeBaseText
		notSeen <- base::setdiff( allBases, xLocs)
		if ( length( notSeen)) {
			where <- MATCH( notSeen, allBases)
			genomeSNPtext[where] <- ""
		}
	}

	# see what bases are different from genomic
	whoSNP <- WHICH( snpTopBase != ",")
	if ( length( whoSNP) > 0) {
		myLocs <- as.integer( names(snpTopBase)[whoSNP])
		where <- MATCH( myLocs, allBases)
		genomeSNPtext[ where] <- snpTopBase[ whoSNP]
	}

	# we may have been asked to only allow SNPs, and disregard any indels
	if (SNP.only) {
		isIndel <- which( nchar(genomeSNPtext) != 1)
		if ( length( isIndel)) genomeSNPtext[ isIndel] <- genomeBaseText[ isIndel]
	}

	# perhaps make the protein amino acid letters too
	# only valid context when we were given a gene
	genomeAminoText <- ""
	if (aaToo && isGene) {
		ans <- convertGenomicBasesToCodingAA( seqID, position=start, end=stop, 
					strand=geneStrand, dnaQuery=genomeSNPtext, genomeDNA=curGenomeDNA,
					geneMap=gmap, cdsMap=cmap)
		genomeAminoText <- ans$query
		names( genomeAminoText) <- allBases
	}

	# we need to make the 'flips' matrix look full when we were asked to use 'genomic' for missing...
	if ( noReadCalls == "genomic") {
		# now merge this with what we were given
		if ( nrow(flips)) {
			drops <- which( rownames(genomicNonFlips) %in% rownames(flips))
			if ( length(drops)) {
				genomicNonFlips <- genomicNonFlips[ -drops, ]
			}
			flips <- rbind( flips, genomicNonFlips)
			ord <- order( rownames(flips))
			flips <- flips[ ord, ]
		} else {
			flips <- genomicNonFlips
		}
	}

	out <- list( "ref"=genomeBaseText, "callsMatrix"=flips, "dna.consensus"=genomeSNPtext, 
			"aa.consensus"=genomeAminoText, "indel.details"=indelText)
	return( out)
}


`consensusBaseCallsToCDNA` <- function( baseCallAns, geneID=NULL, start=NULL, stop=NULL, verbose=TRUE) {

	# get the gene facts we need to make cDNA
	if ( ! is.null(geneID)) {
		gmap <- subset( getCurrentGeneMap(), GENE_ID == geneID)
		cmap <- subset( getCurrentCdsMap(), GENE_ID == geneID)
		if ( any( c( nrow(gmap), nrow(cmap)) == 0)) {
			cat( "\nGeneID not found in current species..")
			cat( "\nFailed to convert to cDNA.")
			baseCallAns$callsTable <- NULL
			return( baseCallAns)
		}
		# grab the bounds and strand info
		strand <- gmap$STRAND[1]
		if ( is.null( start)) start <- gmap$POSITION[1]
		if ( is.null( stop)) stop <- gmap$END[1]
		cdsBases <- vector()
		for (j in 1:nrow(cmap)) cdsBases <- c( cdsBases, (cmap$POSITION[j] : cmap$END[j]))
		cdsBases <- sort( cdsBases)
		# trim to the range we did consensus on
		cdsBases <- cdsBases[ cdsBases %in% (start:stop)]
		names(cdsBases) <- 1:length(cdsBases)
		if (strand == "-") names(cdsBases) <- rev( 1:length(cdsBases))

		# drop the intron data, and do the RevComp if need be
		ref <- baseCallAns$ref
		keep <- which(  names(ref) %in% cdsBases)
		ref <- ref[ keep]
		if ( strand == "-") {
			dna <- paste( ref, collapse="")
			dna <- myReverseComplement(dna)
			ref <- strsplit( dna, split="")[[1]]
		}
		names(ref) <- 1:length(ref)
		baseCallAns$ref <- ref

		calls <- baseCallAns$callsMatrix
		if ( is.null( calls)) {
			cat( "\nNo base calls found..")
			cat( "\nUnable to convert to cDNA.")
			baseCallAns$callsTable <- NULL
			return( baseCallAns)
		}
		indel.details <- baseCallAns$indel.details
		keep <- which( rownames(calls) %in% cdsBases)
		calls <- calls[ keep, ]
		indel.details <- indel.details[ keep]
		if ( strand == "-" && nrow(calls)) {
			# swap A/T and C/G
			tmp <- calls[,"A"]
			calls[,"A"] <- calls[,"T"]
			calls[,"T"] <- tmp
			tmp <- calls[,"C"]
			calls[,"C"] <- calls[,"G"]
			calls[,"G"] <- tmp
			# and lastly flip the row ordering too
			calls <- calls[ rev(1:nrow(calls)), ]
			indel.details <- rev( indel.details)
		}
		# there is no gaurantee that the calls have all base locations, so 
		# rename them specifically
		if ( nrow(calls)) {
			where <- match( rownames(calls), cdsBases)
			rownames(calls) <- names(cdsBases)[where]
		}
		baseCallAns$callsMatrix <- calls

		dna.consensus <- baseCallAns$dna.consensus
		keep <- which(  names(dna.consensus) %in% cdsBases)
		dna.consensus <- dna.consensus[ keep]
		if ( strand == "-") {
			dna <- paste( dna.consensus, collapse="")
			dna <- myReverseComplement(dna)
			dna.consensus <- strsplit( dna, split="")[[1]]
		}
		names(dna.consensus) <- 1:length(dna.consensus)
		baseCallAns$dna.consensus <- dna.consensus
	} else {
		dna.consensus <- baseCallAns$dna.consensus
		indel.details <- baseCallAns$indel.details
		calls <- baseCallAns$callsMatrix
	}

	# the AA can't be converted, re-calc from first principles
	aaAns <- consensusTranslation( dna.consensus)
	aaVec <- aaAns$BestFrame
	names(aaVec) <- 1:length(aaVec)
	baseCallAns$aa.consensus <- aaVec

	# lastly, make a table that combines all the useful facts for resolving the consensus
	# trim to just entries where we have pileup data!
	colnames( calls) <- c( "Ref", colnames(calls)[2:ncol(calls)])
	callBases <- match( rownames(calls), names(dna.consensus))
	callsTable <- data.frame( "POSITION"=as.numeric(rownames(calls)), calls, "DNA"=dna.consensus[callBases], 
				"AA"=aaAns$BestFrame[callBases], aaAns[ callBases, c("Frame1","Frame2","Frame3")], 
				"IndelDetails"=indel.details, stringsAsFactors=F)

	# force this table to be fully translatable, even if we have to trim a few rows
	if (verbose) cat( "  phasing reading frames..")
	if ( nrow(callsTable)) callsTable <- forceValidTranslation( callsTable, verbose=verbose)
	baseCallAns$callsTable <- callsTable
	if (verbose) cat( "  Done.\n")
	return( baseCallAns)
}


`consensusBaseNoiseMasker` <- function( baseCallAns, min.depth=20, min.noise=0.15, dna.mask="N",
					aa.mask="X") {

	# given the base call consensus, allow masking of of high noise calls that may be unreliable,
	# due to either polyclonal material, poor genomic uniqueness, etc.
	if ( is.null(baseCallAns)) return( baseCallAns)
	cm <- baseCallAns$callsMatrix
	if ( is.null(cm)) return( baseCallAns)
	if ( ! nrow(cm)) return( baseCallAns)

	dna <- baseCallAns$dna.consensus
	aa <- baseCallAns$aa.consensus
	hasAA <- (length(dna) == length(aa))

	# see how deep and how noisy
	totalReads <- apply( cm, 1, sum, na.rm=T)
	callReads <- apply( cm, 1, max, na.rm=T)
	noiseReads <- totalReads - callReads

	# test for are we noisy enough
	pctNoise <- noiseReads / totalReads
	isNoise <- which( totalReads >= min.depth & pctNoise >= min.noise)

	out <- baseCallAns
	if ( length(isNoise)) {
		dna[ isNoise] <- dna.mask
		out$dna.consensus <- dna
		if ( hasAA) {
			aa[ isNoise] <- aa.mask
			out$aa.consensus <- aa
		}
	}
	return(out)
}


`consensusTranslation` <- function( dna.consensus) {

	Ndna <- length( dna.consensus)
	dnaNames <- names(dna.consensus)
	aa.consensus <- matrix( "", nrow=Ndna, ncol=3)
	rownames(aa.consensus) <- dnaNames
	colnames(aa.consensus) <- c( "Frame1", "Frame2", "Frame3")

	consensusDNA <- paste( dna.consensus, collapse="")
	consensusAA <- DNAtoAA( consensusDNA, clipAtStop=FALSE, readingFrames=1:3)
	nStops <- sapply( gregexpr( STOP_CODON, consensusAA, fixed=T), length)
	bestFrame <- which.min( nStops)
	aaVector <- strsplit( consensusAA, split="")
	for ( frame in 1:3) {
		theseAA <- aaVector[[frame]]
		aaLocs <- seq( frame, Ndna, by=3) + 1
		aaLocs <- aaLocs[ aaLocs <= Ndna]
		Nch <- min( length( aaLocs), length( theseAA))
		aa.consensus[ aaLocs[1:Nch], frame] <- theseAA[1:Nch]
	}

	bestAns <- aa.consensus[ , bestFrame]
	names( bestAns) <- dnaNames

	return( data.frame( "BestFrame"=bestAns, aa.consensus, stringsAsFactors=FALSE))
}


`forceValidTranslation` <- function( callsTable, verbose=TRUE) {

	# we have the output of Mpileup, with all its base details and AA translations
	# force the entire table to be a valid translation in reading frame 1

	# step 1:  if any indels (length != 1), force everybody to one DNA letter per row
	dnaVec <- callsTable$DNA
	dnaBaseLen <- nchar( dnaVec)
	zeros <- which( dnaBaseLen == 0)
	if ( length( zeros)) {
		if (verbose) cat( "  remove deletions: ", length(zeros))
		callsTable <- callsTable[ -zeros, ]
		dnaVec <- callsTable$DNA
		dnaBaseLen <- nchar( dnaVec)
	}
	oversize <- which( dnaBaseLen > 1)
	if ( nOver <- length( oversize)) {
		if (verbose) cat( "  expand insertions: ", length(nOver))
		# step along and expand each
		headDF <-  tailDF <- data.frame()
		if ( oversize[1] > 1) headDF <- callsTable[ 1:(oversize[1]-1), ]
		if ( oversize[nOver] < nrow(callsTable)) tailDF <- callsTable[ (oversize[nOver]+1):nrow(callsTable), ]
		outDF <- headDF
		for (i in 1:nOver) {
			thisRow <- oversize[i]
			myDNA <- strsplit( dnaVec[thisRow], split="")[[1]]
			nDNA <- length(myDNA)
			sml <- callsTable[ rep.int(thisRow,nDNA), ]
			sml$DNA <- myDNA
			sml$POSITION <- round( as.numeric(callsTable$POSITION[thisRow]) + (((1:nDNA) - 1) / nDNA), digits=3)
			sml$AA <- ""
			outDF <- rbind( outDF, sml)
			# if more indels, append the next non-indel segment
			if ( i != nOver) {
				otherDF <- callsTable[ (thisRow+1):(oversize[i+1]-1), ]
				outDF <- rbind( outDF, otherDF)
			}
		}
		callsTable <- rbind( outDF, tailDF)
		rm( headDF, tailDF, sml, outDF)
		dnaVec <- callsTable$DNA
		dnaBaseLen <- nchar( dnaVec)
	}

	# step 2:  find all the best longest reading frames
	dnaStr <- paste( dnaVec, collapse="")
	if ( nchar(dnaStr) != length(dnaVec)) {
		if (verbose) cat( "\nWarning: 'ConsensusBaseCalls' DNA sequence length error!")
	}
	ans <- DNAtoFrameShiftingPeptides( dnaStr, min.aa.length=2)

	# now re-concatenate just those chunks
	if (verbose) cat( "  concatenating chunks: ", nrow(ans))
	outDF <- data.frame()
	for ( i in 1:nrow(ans)) {
		from <- ans$DNA_Start[i]
		to <- ans$DNA_Stop[i]
		sml <- callsTable[ from:to, ]
		outDF <- rbind( outDF, sml)
	}

	# test it
	if (verbose) cat( "  validate..")
	testStr <- paste( outDF$DNA, collapse="")
	testAA <- DNAtoAA( testStr, clip=F, readingFrame=1)
	nStops <- sum( gregexpr( STOP_CODON_PATTERN, testAA)[[1]] > 0)
	if ( nStops > 0) {
		cat( "\nWarning: 'ConsensusBaseCalls' removing stop codons failed.. N_Stops: ", nStops)
	}

	# store it
	outDF$AA <- ""
	aaVec <- strsplit( testAA, split="")[[1]]
	outDF$AA[ seq( 2, nrow(outDF), by=3)] <- aaVec
	
	# and restore all the reading frames
	aaAns <- consensusTranslation( outDF$DNA)
	outDF[ , c("Frame1", "Frame2", "Frame3")] <- aaAns[ , c("Frame1","Frame2","Frame3")]

	return( outDF)
}
