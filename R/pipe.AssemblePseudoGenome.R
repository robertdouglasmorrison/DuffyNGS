# pipe.AssemblePseudoGenome.R

# turn raw FASTQ reads into a 'complete' genome, using SPAdes to build contigs and then
# BLAST to decide how they align & organize


# turn raw FASTQ into a Pseudo Genome FASTA and Annotation for one sample

`pipe.AssemblePseudoGenome` <- function( sampleID=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt", 
				doCutadapt=TRUE, cutadaptProgram=Sys.which("cutadapt"), forceMatePairs=NULL, 
				doSpades=FALSE, spades.path=dirname(Sys.which("spades.py")), 
				spades.mode=c("isolate","rna","meta"), spades.kmer.sizes=NULL, spades.args="",
				doBlast=FALSE, verbose=TRUE) {

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nStarting  'AssemblePseudoGenome' on Sample:     ", sampleID,
			"\nCurrent Species: ", getCurrentSpecies(), 
			"\nStart Date/Time:   \t", date(), "\n")
	}
	startT <- proc.time()
	require(Biostrings)
	if (version$major == "4" && as.numeric( version$minor) >= 4) require( pwalign)

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
	
	# Step 5:  try to find unexpected arrangements of genes
	out <- findUnexpectedArrangements( spades.out.path, sampleID=sampleID)
	
	
	cat( verboseOutputDivider)
	cat( "\n\nFinished 'AssemblePseudoGenome' on Sample:     ", sampleID, "\n")
	cat( "\nTiming Stats: \n")
	myTime <- elapsedProcTime( startT, proc.time(), N=nrow(out), what="Mapped Genes")
	print( myTime)
	
	return( invisible( out))
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
	gidOut <- posOut <- endOut <- sidOut <- strandOut <- prodOut <- contigidOut <- vector()
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
				contigidOut[nout] <- smlAns$SEQ_ID[1]
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
			"STRAND"=strandOut, "PRODUCT"=prodOut, "CONTIG_ID"=contigidOut, stringsAsFactors=T)
	rownames(outMap) <- 1:nout
	# add a metric of relative position of each gene to the true location on the chromosome
	whMap <- match( outMap$GENE_ID, gmap$GENE_ID)
	deltaDist <- outMap$POSITION - gmap$POSITION[ whMap]
	outMap$SHIFT_DISTANCE <- deltaDist
	write.table( outMap, outGeneMap, sep="\t", quote=F, row.names=F)
	return( invisible( outMap))
}


`findUnexpectedArrangements` <- function( out.path, sampleID, min.move.dist=100, min.gap.size=100) {

	# given the assembled pseudo genome and gene map results, see if we can find any unexpected items
	
	# Step 0: set up:  get the contigs, pseudo genome, and pseudo gene map
	contigFile <- file.path( out.path, paste( sampleID, "BlastContigs.fasta", sep="."))
	contigFasta <- loadFasta( contigFile, verbose=F)
	NCONTIG <- length( contigFasta$desc)
	# file names for the results we just generated, and read them back in
	inGenome <- file.path( out.path, paste( sampleID, "PseudoGenome.fasta", sep="."))
	myGenomeFA <- loadFasta( inGenome, verbose=F)
	inGeneMap <- file.path( out.path, paste( sampleID, "PseudoGenome.GeneMap.txt", sep="."))
	myGmap <- read.delim( inGeneMap, as.is=T)
	NG <- nrow(myGmap)
	# re-get the set of reference genome genes, that we will look for in the contigs
	refGmap <- subset( getCurrentGeneMap(), REAL_G == TRUE)
	
	# Step 1:  see if any gene inside a conserved contig seems way out of place
	# 		use a set of 5 adjacent genes, inside the same raw contig, as the neighborhood.
	#		check the center one
	myGmap$Gene.Moved <- ""
	for ( g in refGmap$GENE_ID) {
		myRow <- match( g, myGmap$GENE_ID, nomatch=0)
		if ( ! myRow) next
		myNeighborRows <- (max(myRow-2, 1) : min(myRow+2, NG))
		myNeighborContigs <- myGmap$CONTIG_ID[ myNeighborRows]
		# all must be the same to trust looking at the gene
		if ( length( unique( myNeighborContigs)) > 1) next
		# see if this gene is 'more distant' than its peers
		myDist <- myGmap$SHIFT_DISTANCE[myRow]
		myNeighborDists <- myGmap$SHIFT_DISTANCE[ setdiff(myNeighborRows,myRow)]
		myDistDifferent <- setdiff( myDist, myNeighborDists)
		myDistDiff <- myDist - median( myNeighborDists)
		if ( length( myDistDifferent) && abs(myDistDiff) >= min.move.dist) {
			myGmap$Gene.Moved[ myRow] <- paste( "Moved", myDistDiff)
		}
	}
	
	# Step 2:  see if any gaps without genes inside a conserved contig seems way different size
	# 		use a set of 2 adjacent genes, inside the same raw contig, as the neighborhood.
	#		check the gap in both this map and the reference
	myGmap$Post.Gene.Gap.Feature <- ""
	myGmap$Post.Gene.Gap.Blast.Hit <- ""
	# accumulate novel DNA to blast against NT
	extraDesc <- extraSeq <- vector()
	for ( ig in 1:(nrow(refGmap)-1)) {
		g1 <- refGmap$GENE_ID[ig]
		g2 <- refGmap$GENE_ID[ig+1]
		myRow1 <- match( g1, myGmap$GENE_ID, nomatch=0)
		myRow2 <- match( g2, myGmap$GENE_ID, nomatch=0)
		if ( ! myRow1 || ! myRow2 ) next
		if ( myRow2 - myRow1 != 1) next
		if ( myGmap$CONTIG_ID[myRow1] != myGmap$CONTIG_ID[myRow2]) next
		# we have 2 adjacent genes in the same contig.  Check for large change in size.
		refGapSize <- refGmap$POSITION[ig+1] - refGmap$END[ig]
		if ( refGapSize < min.gap.size) next
		myGapSize <- myGmap$POSITION[myRow2] - myGmap$END[myRow1]
		if ( myGapSize < refGapSize/2) {
			myGmap$Post.Gene.Gap.Feature[myRow1] <- paste( "Post-Gene Gap is Smaller by", refGapSize-myGapSize)
			next
		}
		if ( myGapSize > refGapSize*1.5) {
			myGmap$Post.Gene.Gap.Feature[myRow1] <- paste( "Post-Gene Gap is Larger by", myGapSize-refGapSize)
			# since it is bigger, let's try to see what it is
			thisSeqID <- myGmap$SEQ_ID[myRow1]
			thisChromo <- myGenomeFA$seq[ match( thisSeqID, myGenomeFA$desc)]
			thisGapSeq <- substr( thisChromo, myGmap$END[myRow1] + 5, myGmap$POSITION[myRow2] - 5)
			extraDesc <- c( extraDesc, g1)
			extraSeq <- c( extraSeq, thisGapSeq)
			next
		}
	}
	
	# if we got some extra DNA, try to blast it
	if ( length( extraSeq)) {
		cat( "\n  BLASTing ", length( extraSeq), " unexpected DNA..")
		tmpFastaFile <- paste( sampleID, "Unexpected.Gaps.fasta", sep=".")
		writeFasta( as.Fasta( extraDesc, extraSeq), tmpFastaFile, line=100)
		tmpBlastFile <- paste( sampleID, "Unexpected.Gaps.Blast.Output.txt", sep=".")
		# because we build a Unix command line, embed the outfmt string in quotes
		#callBlastn( tmpFastaFile, outfile=tmpBlastFile, wordsize=11, evalue=0.1, maxhits=3, outfmt=' "6 std stitle" ')
		callBlastx( tmpFastaFile, outfile=tmpBlastFile, wordsize=4, evalue=0.1, maxhits=3, outfmt=' "6 std stitle" ')
		SAV_BLAST_ANS <<- ans <- readBlastOutput( tmpBlastFile, nKeep=1, outfmt="6 std stitle")
		# make an answer string from the result
		blastText <- paste( ans$SEQ_ID, ", Title=", ans$STITLE, "; Length=", ans$LEN_MATCH, "; Score=", ans$SCORE, sep="")
		where <- match( ans$PROBE_ID, myGmap$GENE_ID, nomatch=0)
		myGmap$Post.Gene.Gap.Blast.Hit[ where] <- blastText[ where > 0]
		#file.delete( c( tmpFastaFile, tmpBlastFile))
	}
	
	write.table( myGmap, inGeneMap, sep="\t", quote=F, row.names=F)
	return( invisible( myGmap))
}
