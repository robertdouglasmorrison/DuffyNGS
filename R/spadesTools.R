# spadesTools.R


`makeSpadesContigs` <- function( fastaFile, outpath, spades.path=dirname(Sys.which("spades.py")), 
				spades.mode=c("isolate","rna","meta"), kmerSizes=NULL, spades.args="", 
				doPairedEnd=FALSE, verbose=FALSE) {

	# in general, SPAdes wants an empty folder
	if ( file.exists( outpath)) {
		if ( file.info( outpath)$isdir) {
			system( paste( "rm -rf", outpath))
		} else { 
			file.delete( outpath)
		}
	}

	doInternal <- ! verbose

	# generate the command line args needed...
	NF <- length( fastaFile)
	if (doPairedEnd) {
		if ( NF != 2) stop( "'doPairedEnd' requires exactly 2 fastq files..")
		filesStr <- paste( "-1", fastaFile[1], "-2", fastaFile[2], sep=" ")
	} else {
		filesStr <- paste( "-s", fastaFile, sep=" ", collapse=" ")
	}

	spades.mode <- match.arg( spades.mode)
	modeStr <- paste( "--", spades.mode, sep="")

	kmerStr <- ""
	if ( ! is.null( kmerSizes)) kmerStr <- paste( "-k", paste( kmerSizes, collapse=" "))

	# put those args together
	cmdline <- paste( file.path( spades.path, "spades.py"), modeStr, filesStr, 
				kmerStr, spades.args, "-o", outpath, "  2>&1", sep=" ")

	cat( "\n\nCalling SPAdes...\n")
	if (verbose) cat( "\nCommand Line: ", cmdline, "\n")
	#SPAdes is very verbose, never show it..
	systemAns <- catch.system( cmdline, intern=TRUE)
	cat( "  Done.\n")
}


`makeSpadesPeptides` <- function( sampleID, outpath, keyword="Spades", short.desc=FALSE, verbose=TRUE) {

	require( Biostrings)

	# grab the final set of contigs from the Spades run
	contigFile <- file.path( outpath, "contigs.fasta")
	if ( ! file.exists( contigFile)) {
		cat( "\nSpades 'contigs.fasta' file not found:  ", contigFile)
		return(0)
	}
	contigs <- loadFasta( contigFile, verbose=verbose)
	cat( "\nMaking Peptides from contigs..\nN_Contigs: ", length( contigs$desc))
	if ( length( contigs$desc) < 1) return(0)
	desc <- contigs$desc
	if (short.desc) desc <- sub( "_length.+", "", desc)

	# convert to peptides
	#peps <- DNAtoBestPeptide( contigs$seq, clipAtStop=FALSE)
	mcAns <- multicore.lapply( contigs$seq, FUN=DNAtoBestPeptide, clipAtStop=FALSE, preschedule=TRUE)
	peps <- unlist( mcAns)
	rfName <- sapply( mcAns, function(x) names(x)[1])
	faOut <- as.Fasta( paste( sampleID, desc, rfName, sep="_"), peps)

	outfile <- paste( sampleID, keyword, "Peptides.fasta", sep=".")
	outfile <- file.path( outpath, outfile)
	writeFasta( faOut, file=outfile, line.width=100) 
	if (verbose) cat( "\nWrote Peptides file: ", basename(outfile))

	ncharPeps <- nchar( peps)
	meanPepLen <- mean( ncharPeps)
	if (verbose) {
		pcts <- quantile( ncharPeps, seq(0,1,0.1))
		cat( "\nPeptile Lengths, by quantile:\n")
		print( pcts)
		cat( "\nMean Peptide Length:  ", meanPepLen)
		cat( "\nTotal Amino Acids:    ", prettyNum( sum( ncharPeps), big.mark=","),"\n")
	}

	return( meanPepLen)
}


`bestSpadesProteins` <- function( sampleID, outpath, keyword="Spades", proteinFastaFile, verbose=TRUE) {

	require( Biostrings)
	data( PAM70)

	# grab the set of peptides from the Velvet run
	peptideFile <- file.path( outpath, paste( sampleID, keyword, "Peptides.fasta", sep="."))
	if ( ! file.exists( peptideFile)) {
		cat( "\nSpades 'peptides' file not found:  ", peptideFile)
		return(NULL)
	}
	if ( ! file.exists( proteinFastaFile)) {
		cat( "\nFile of proteins not found:  ", proteinFastaFile)
		return(NULL)
	}

	peptides <- loadFasta( peptideFile, verbose=verbose)
	desc <- peptides$desc
	seqs <- peptides$seq
	N <- length( desc)
	if ( N < 1) return(NULL)
	if (verbose) cat( "\nFinding best Proteins:  ", basename(peptideFile), "\nN_Peptides: ", N, "\n")

	# map peptides to proteins
	ans <- peptide2BestProtein( seqs, proteinFastaFile, substitutionMatrix=PAM70,
				details=TRUE, verbose=verbose)
	if ( is.null(ans)) return(NULL)


	# unwrap the list output
	gid <- score <- scorepaa <- starts <- stops <- editdist <- weight <- vector( length=N)

	# should be a list of lists, but if only one peptide, its not...
	if (N == 1) {
		tmp <- ans
		ans <- vector( mode="list")
		ans[[1]] <- tmp
	}

	for ( i in 1:N) {
		x <- ans[[i]]
		if ( is.null(x)) next
		gid[i] <- x$ProtName
		score[i] <- x$Score
		scorepaa[i] <- x$ScorePerAA
		starts[i] <- x$ProtStart
		stops[i] <- x$ProtStop
		editdist[i] <- x$EditDistance
		weight[i] <- x$Weight
	}

	out <- data.frame( "ContigID"=peptides$desc, "AA.Length"=nchar(peptides$seq), "Best.Match.ProteinID"=gid, "Score"=score, 
			"ScorePerAA"=round(scorepaa,digits=3), "ProteinStart"=starts, "ProteinStop"=stops, 
			"EditDistance"=editdist, "Weight"=round(weight,digits=3), stringsAsFactors=F)
	drops <- which( is.na( out$Score))
	if ( length(drops)) out <- out[ -drops, ]

	if ( nrow( out)) {
		ord <- order( out$Score, decreasing=T)
		out <- out[ ord, ]
		rownames(out) <- 1:nrow(out)
	}

	outfile <- paste( sampleID, keyword, "BestProteinHits.txt", sep=".")
	outfile <- file.path( outpath, outfile)
	write.table( out, file=outfile, sep="\t", quote=F, row.names=F)
	if (verbose) cat( "\nWrote Proteins file: ", basename(outfile))

	return( invisible( out))
}

