# pipe.AssemblePseudoGenome.R

# turn raw FASTQ reads into a 'complete' genome, using SPAdes to build contigs and then
# BLAST to fing how they organize



# turn raw FASTQ into Pseudon Genome FASTA and Annotation for one sample

`pipe.AssemblePseudoGenome` <- function( sampleID=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt", 
				doCutadapt=TRUE, cutadaptProgram=Sys.which("cutadapt"), forceMatePairs=NULL, 
				doSpades=FALSE, spades.path=dirname(Sys.which("spades.py")), 
				spades.mode=c("isolate","rna","meta"), spades.kmer.sizes=NULL, spades.args="",
				doBlast=FALSE, verbose=TRUE) {

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nStarting  'AssemblePseudoGenome' on Sample:     ", sampleID,
			"\nStart Date/Time:   \t", date(), "\n")
	}
	startT <- proc.time()
	require(Biostrings)

	# FASTQ file(s) to process comes from annotation and options...
	results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	rawFastq <- getRawFastqFileNames( sampleID, annotationFile, optionsFile, verbose=FALSE)
	inputFastqFiles <- rawFastq$files
	asMatePairs <- rawFastq$asMatePairs
	if ( ! is.null(forceMatePairs) && forceMatePairs == TRUE) asMatePairs <- TRUE

	# Step 1:  perhaps use CutAdapt to clean raw reads
	# if we are doing the Cutadapt pass, check for those files, and call it if we need.
	if (doCutadapt) {
		# this called function is in the 'Kmer' tools...
		filesToDo <- checkOrCallCutadapt( inputFastqFiles, asMatePairs=asMatePairs, forceMatePairs=forceMatePairs,
						cutadaptProgram=cutadaptProgram, kmer.size=50, verbose=verbose)
	} else { 
		filesToDo <- inputFastqFiles
	}
	nFiles <- length( filesToDo)
	
	# Step 2:  run SPAdes to build contigs
	# spades needs a folder for writing all its results
	spades.program.path <- spades.path
	spades.out.path <- file.path( results.path, "SpadesContigs", sampleID, "PseudoGenome")
	if ( ! file.exists( spades.out.path)) dir.create( spades.out.path, recursive=T)	
	spadesOutFile <- file.path( spades.out.path, "contigs.fasta")
	if (doSpades || ! file.exists( spadesOutFile)) {
		spades.mode <- match.arg( spades.mode)
		doPairedEnd <- (asMatePairs && nFiles == 2)
		# make sure those FASTQ files exist
		fqExist <- file.exists( filesToDo)
		if ( ! all( fqExist)) {
			cat( "\nError:  some FASTQ files for SPAdes not found:\n")
			cat( filesToDo[ ! fqExist], "\n")
			return(NULL)
		}
		spadesAns <- makeSpadesContigs( filesToDo, outpath=spades.out.path, spades.path=spades.program.path, 
			spades.mode=spades.mode, kmerSizes=spades.kmer.sizes, 
			spades.args=spades.args, doPairedEnd=doPairedEnd, verbose=verbose)
	} else {
		if (verbose) cat( "\nUsing existings SPAdes contigs..")
	}

	# Step 3:  blast reference genes against the contigs
	genome.file <- getOptionValue( optionsFile, "genomicFastaFile", verbose=FALSE)
	blastOutfile <- file.path( spades.out.path, paste( sampleID, "BlastOut.txt", sep="."))
	if ( doBlast || !file.exists(blastOutfile)) {
		blastAns <- blastGenesToContigs( spades.out.path, sampleID=sampleID, 
						genome.file=genome.file, mode="first.pass", verbose=verbose)
	} else {
		if (verbose) cat( "\nUsing existings BLAST results..")
		blastAns <- readBlastOutput( blastOutfile, nKeep=1)
	}
	
	# Step 4:  now assemble these contigs, using reference genes as the guide
	out <- assembleContigsToGenome( spades.out.path, sampleID=sampleID, genome.file=genome.file, 
					blastAns=blastAns, verbose=verbose)
	
	cat( verboseOutputDivider)
	cat( "\n\nFinished 'AssemblePseudoGenome' on Sample:     ", sampleID, "\n")
	cat( "\nTiming Stats: \n")
	myTime <- elapsedProcTime( startT, proc.time(), N=nrow(out), what="Mapped Genes")
	print( myTime)
	
	return( out)
}


`blastGenesToContigs` <- function( out.path, sampleID, genome.file, mode=c("first.pass","final"), verbose=TRUE) {

	# given a folder of SPAdes contigs, and the reference genes, make a BLAST result 
	# that maps genes to contigs
	
	mode <- match.arg( mode)
	
	if ( mode == "first.pass") {
		# start from the final set of contigs from the Spades run
		contigFile <- file.path( out.path, "contigs.fasta")
		if ( ! file.exists( contigFile)) {
			cat( "\nSpades 'contigs.fasta' file not found:  ", contigFile)
			return(0)
		}
	} else {
		# start from the curated set of contigs from the Spades run
		contigFile <- file.path( out.path, paste( sampleID, "BlastContigs.fasta", sep="."))
		if ( ! file.exists( contigFile)) {
			cat( "\nCurated Spades 'BlastContigs.fasta' file not found:  ", contigFile)
			return(0)
		}
	}
	contigs <- loadFasta( contigFile, verbose=verbose)

	# BLAST does both strands for us, so use the Beg/End locations to deduce true strand
	contigSeqs1 <- contigs$seq
	contigDesc1 <- contigs$desc

	# get the set of reference genome genes, that we will look for in the contigs
	gmap <- subset( getCurrentGeneMap(), REAL_G == TRUE)
	genesFasta <- gene2Fasta( gmap$GENE_ID, genomicDNAfilePath=genome.file, mode="gdna")
	genesToFind <- genesFasta$desc
	
	# set up to use BLAST, to find the best hit to every gene, looking at all the contigs in both directions
	tmpFastaC <- as.Fasta( contigDesc1, contigSeqs1)
	tmpFastaCfile <- file.path( out.path, paste( sampleID, "BlastContigs.fasta", sep="."))
	if ( mode == "first.pass") writeFasta( tmpFastaC, tmpFastaCfile, line.width=100)
	genesFastaFile <- file.path( out.path, paste( sampleID, "ReferenceGenes.fasta", sep="."))
	writeFasta( genesFasta, genesFastaFile, line.width=100)
	
	# make the BLAST index, and find the one best hit for each gene
	dbname <- sub( ".fasta$", "", basename(tmpFastaCfile))
	blastOutfile <- file.path( out.path, paste( sampleID, "BlastOut.txt", sep="."))
	makeBlastDB( tmpFastaCfile, dbname=dbname, dbtype="nucl")
	callBlastn( genesFastaFile, outfile=blastOutfile, db=dbname, wordsize=11, evalue=0.1, maxhits=1, verbose=verbose)
	blastAns <- readBlastOutput( blastOutfile, nKeep=1)
	return( blastAns)
}


`assembleContigsToGenome` <- function( out.path, sampleID, genome.file, blastAns, verbose=TRUE) {

	# given the blast results, turn them into a pseudo genome
	
	# Step 0: set up 
	# get the contigs that the Blast step aligned against
	contigFile <- file.path( out.path, paste( sampleID, "BlastContigs.fasta", sep="."))
	contigFasta <- loadFasta( contigFile, verbose=F)
	NCONTIG <- length( contigFasta$desc)
	# file names for the results we will generate, and removed any previous results
	outGenome <- file.path( out.path, paste( sampleID, "PseudoGenome.fasta", sep="."))
	outGeneMap <- file.path( out.path, paste( sampleID, "PseudoGenome.GeneMap.txt", sep="."))
	file.delete( c( outGenome, outGeneMap))	
	# re-get the set of reference genome genes, that we will look for in the contigs
	gmap <- subset( getCurrentGeneMap(), REAL_G == TRUE)
	genesFasta <- gene2Fasta( gmap$GENE_ID, genomicDNAfilePath=genome.file, mode="gdna")
	genesToFind <- genesFasta$desc
	NG <- length( genesToFind)
	NSEQ <- length( unique( gmap$SEQ_ID))
	# note that a genome of only 1 chromosome is likely circular, so allow some contigs to get used twice
	isCircularGenome <- (NSEQ == 1)
	
	# Step 1:  deduce which contigs need to be rev comp'ed
	cat( "\nAssessing RevComp needs of Contigs..")
	votesFwdContig <- votesRevContig <- rep.int( 0, NCONTIG)
	votesFwdString <- votesRevString <- rep.int( "", NCONTIG)
	for( i in 1:nrow(gmap)) {
		g <- gmap$GENE_ID[i]
		gStrand <- gmap$STRAND[i]
		smlAns <- subset( blastAns, PROBE_ID == g)
		if ( ! nrow(smlAns)) next
		if ( nrow(smlAns) > 1) smlAns <- smlAns[ which.max( smallAns$SCORE), ]
		myContig <- smlAns$SEQ_ID[1]
		myPtr <- match( myContig, contigFasta$desc)
		myFrom <- smlAns$S_BEG[1]
		myTo <- smlAns$S_END[1]
		myStrand <- smlAns$STRAND[1]
		if (gStrand == myStrand) {
			votesFwdContig[myPtr] <- votesFwdContig[myPtr] + 1
			votesFwdString[myPtr] <- paste( votesFwdString[myPtr], g, sep=",")
		} else {
			votesRevContig[myPtr] <- votesRevContig[myPtr] + 1
			votesRevString[myPtr] <- paste( votesRevString[myPtr], g, sep=",")
		}
	}
	# now flip those that need it and resave
	newContigFasta <- contigFasta
	cat( "\nContigs needing RevComp: ", sum( votesRevContig > votesFwdContig))
	for ( k in 1:NCONTIG) {
		if ( votesRevContig[k] > votesFwdContig[k]) {
			newContigFasta$seq[k] <- myReverseComplement( newContigFasta$seq[k])
			newContigFasta$desc[k] <- paste( "RevComp", newContigFasta$desc[k], sep="_")
		}
	}
	# the SPAdes contig names end with coverage, as a floating point value.  Truncate some
	newContigFasta$desc <- sub( "\\.[0-9]+$", "", newContigFasta$desc)
	writeFasta( newContigFasta, contigFile, line.width=100)
	contigFasta <- newContigFasta
	
	# Step 2: redo the BLAST against the final 'some got RevComped' contigs
	cat( "\nRerunning BLAST against final Contigs..")
	blastAns <- blastGenesToContigs( out.path, sampleID=sampleID, 
					genome.file=genome.file, mode="final", verbose=verbose)
					
	# Step 2.5:  if this is a circular genome, split the contig that straddles the start into 2
	if (isCircularGenome) {
		# let's force the check that at least 2 genes at each end share one contig
		testGenes <- gmap$GENE_ID[ c( 1,2,NG-1,NG)]
		contigPtrs <- match( testGenes, blastAns$PROBE_ID)
		testIDs <- blastAns$SEQ_ID[ contigPtrs]
		if ( length( unique( testIDs)) == 1) {
			# we have proof of circle.  Break it in two
			cat( "\nSplitting Contig that circularizes genome, and re-BLAST-ing..")
			myFastaPtr <- match( testIDs[1], contigFasta$desc)
			myContigSeq <- contigFasta$seq[ myFastaPtr]
			myContigDesc <- contigFasta$desc[ myFastaPtr]
			firstContig <- substring( myContigSeq, blastAns$S_BEG[ contigPtrs[1]])
			lastContig <- substr( myContigSeq, 1, (blastAns$S_BEG[ contigPtrs[1]]-1))
			firstDesc <- paste( "UncircleHead", myContigDesc, sep="_")
			lastDesc <- paste( "UncircleTail", myContigDesc, sep="_")
			newContigDesc <- c( contigFasta$desc[1:(NCONTIG-1)], lastDesc, firstDesc)
			newContigSeq <- c( contigFasta$seq[1:(NCONTIG-1)], lastContig, firstContig)
			newContigFasta <- as.Fasta( newContigDesc, newContigSeq)
			writeFasta( newContigFasta, contigFile, line.width=100)
			contigFasta <- newContigFasta
			NCONTIG <- length( contigFasta$desc)
			blastAns <- blastGenesToContigs( out.path, sampleID=sampleID, 
						genome.file=genome.file, mode="final", verbose=verbose)
		}
	}

	# Step 3:  Build the new genome and map, one chromosome at a time
	contigHits <- rep.int( 0, NCONTIG)
	gidOut <- posOut <- endOut <- sidOut <- strandOut <- prodOut <- notesOut <- vector()
	nout <- 0
	chromoIdOut <- chromoSeqOut <- vector()
	nchromo <- 0
	SPACER <- "NNNNNNNNNN"
	sidSet <- unique( gmap$SEQ_ID)
	for ( sid in sidSet) {
		# visit each gene in this chromosome, and find it's contig, building up info about 
		# each contig's usage and location within the chromosome.  
		# Use the 'SCORE' detail to say haw many votes each gene gives to its contig
		cat( "\nAssessing Contig ordering for: ", sid, "\n")
		contigVotes <- vector( mode="list", length=NCONTIG)
		for ( k in 1:NCONTIG) contigVotes[[k]] <- numeric(0)
		smlMap <- subset( gmap, SEQ_ID == sid)
		for( i in 1:nrow(smlMap)) {
			gid <- smlMap$GENE_ID[i]
			smlAns <- subset( blastAns, PROBE_ID == gid)
			if ( ! nrow(smlAns)) next
			if ( nrow(smlAns) > 1) smlAns <- smlAns[ which.max( smallAns$SCORE), ]
			myContig <- smlAns$SEQ_ID[1]
			myPtr <- match( myContig, contigFasta$desc)
			myScore <- smlAns$SCORE[1]
			myV <- contigVotes[[myPtr]]
			myV <- c( myV, rep.int( i, round(myScore)))
			contigVotes[[myPtr]] <- myV
			contigHits[myPtr] <- contigHits[myPtr] + 1
			cat( "\rGene contig score: ", gid, myPtr, myContig, myScore, length(myV))
		}
		# summarize the relative placement of each contig
		# at this point, only use contigs that got votes
		avgVotes <- sapply( contigVotes, mean, na.rm=T)
		names(avgVotes) <- contigFasta$desc
		avgVotes[ is.nan(avgVotes)] <- 0
		avgVotes[ is.na(avgVotes)] <- 0
		drops <- which( avgVotes == 0)
		if ( length(drops)) avgVotes <- avgVotes[ -drops]
		avgVotes <- sort( avgVotes)
		# now build the chromosome, one contig at a time, and map those genes as we go
		chromoSeq <- ""
		for ( thisContigDesc in names(avgVotes)) {
			whContigFasta <- match( thisContigDesc, contigFasta$desc)
			thisContigSeq <- contigFasta$seq[ whContigFasta]
			contigHits[ whContigFasta] <- contigHits[ whContigFasta] + 1
			seqLenPrev <- nchar( chromoSeq)
			chromoSeq <- paste( chromoSeq, thisContigSeq, sep="", collapse="")
			# get the genes on this chromosome that mapped to this contig, in the current gene map's order
			genesToMap <- intersect( smlMap$GENE_ID, subset( blastAns, SEQ_ID == thisContigDesc)$PROBE_ID)
			for (gid in genesToMap) {
				# add this gene to our new map as we grow the chromosome
				whMap <- match( gid, gmap$GENE_ID)
				smlAns <- subset( blastAns, SEQ_ID == thisContigDesc & PROBE_ID == gid)
				if ( ! nrow(smlAns)) next
				if ( nrow(smlAns) > 1) smlAns <- smlAns[ which.max( smallAns$SCORE), ]
				nout <- nout + 1
				gidOut[nout] <- gid
				sidOut[nout] <- sid
				# for now, size the new gene based on the BLAST, not the gene's true length
				posOut[nout] <- seqLenPrev + smlAns$S_BEG[1]
				endOut[nout] <- seqLenPrev + smlAns$S_END[1]
				strandOut[nout] <- smlAns$STRAND[1]
				prodOut[nout] <- gmap$PRODUCT[whMap]
				notesOut[nout] <- smlAns$SEQ_ID[1]
			}
			# after each contig has been added, append a spacer
			chromoSeq <- paste( chromoSeq, SPACER, sep="", collapse="")
		}
		# each contig for this chromosome has been dealt with.  Save for final output
		nchromo <- nchromo + 1
		chromoIdOut[nchromo] <- sid
		chromoSeqOut[nchromo] <- chromoSeq
	}
	# for now, ignore any contigs that had zero genes mapped to them, as unimportant for now.
	# prehaps revisit this issue later
	
	outFasta <- as.Fasta( chromoIdOut, chromoSeqOut)
	writeFasta( outFasta, outGenome, line=100)
	outMap <- data.frame( "GENE_ID"=gidOut, "POSITION"=posOut, "END"=endOut, "SEQ_ID"=sidOut,
			"STRAND"=strandOut, "PRODUCT"=prodOut, "CONTIG_NOTES"=notesOut, stringsAsFactors=T)
	rownames(outMap) <- 1:nout
	# add a metric of relative position of each gene to the true location on the chromosome
	whMap <- match( outMap$GENE_ID, gmap$GENE_ID)
	deltaDist <- outMap$POSITION - gmap$POSITION[ whMap]
	outMap$SHIFT_DISTANCE <- deltaDist
	write.table( outMap, outGeneMap, sep="\t", quote=F, row.names=F)
	return( invisible( outMap))
}

