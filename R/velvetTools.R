# velvetTools.R


`makeVelvetContigs` <- function( fastaFile, outpath=sub( ".fa.*","", fastaFile),
		velvet.path="~/bin", buildHash=TRUE, buildContigs=TRUE,
		kmerSize=25, minLength=200, doPairedEnd=FALSE, minCoverage=NULL, 
		velveth.args="", velvetg.args="", verbose=FALSE) {

	if ( ! file.exists( outpath)) dir.create( outpath, recursive=T)

	doInternal <- ! verbose

	# step 1:  VelvelH to build the Velvet hash table
	roadmapFile <- file.path( outpath, "Roadmaps")
	if ( buildHash || !file.exists(roadmapFile)) {

		NF <- length( fastaFile)
		# deduce the args needed...
		format <- "-fasta"
		if ( regexpr( "fastq|fq", fastaFile[1]) > 0) format <- "-fastq"
		if ( regexpr( "q.gz", fastaFile[1]) > 0) format <- "-fastq.gz"

		pairedEndArgs <- ""
		if (doPairedEnd) {
			if ( NF != 2) stop( "'doPairedEnd' requires exactly 2 fastq files..")
			pairedEndArgs <- " -shortPaired -separate "
		}
		if ( NF > 1) {
			fastaFile <- paste( fastaFile, collapse="  ")
		}

		# put those args together
		cmdline <- paste( file.path( velvet.path, "velveth"), outpath, kmerSize, format, 
				pairedEndArgs, fastaFile, velveth.args, sep=" ")
		cat( "\n\nCalling VelvetH to build kmer hash table...\n")
		if (verbose) cat( "\nCommand Line: ", cmdline, "\n")
		systemAns <- catch.system( cmdline, intern=doInternal)
		if ( doInternal) {
			lines2show <- grep( "Reading FastQ", systemAns)
			lines2show <- c( lines2show, grep( "sequences found", systemAns)[1:NF])
			print( systemAns[lines2show])
		}
	}

	# step 2:  VelvetG to build the contigs
	bestCutoff <- NA
	velvetMetrics <- NULL
	graphFile <- file.path( outpath, "Graph2")

	if ( velvetg.args == "") velvetg.args <- " -exp_cov auto "

	if (buildContigs || !file.exists(graphFile)) {
		cmdline <- paste( file.path( velvet.path, "velvetg"), outpath, velvetg.args)
		cat( "\nCalling VelvetG to build contigs...\n")
		if (verbose) cat( "\nCommand Line: ", cmdline, "\n")
		systemAns <- catch.system( cmdline, intern=doInternal)
		if ( doInternal) {
			estCoverLine <- grep( "Estimated Coverage =", systemAns, value=T)[1]
			cutCoverLine <- grep( "Estimated Coverage cutoff =", systemAns, value=T)[1]
			finalGraphLine <- grep( "Final graph has", systemAns, value=T)[1]
			metricsText <- c( estCoverLine, cutCoverLine, finalGraphLine)
			metricsText[ is.na( metricsText)] <- ""
			velvetMetrics <- extractVelvetgMetrics( metricsText)
			bestCutoff <- velvetMetrics$CoverCutoff[1]
		}
	}

	# final call is to extract just the 'good' contigs
	cmdline <- paste( file.path( velvet.path, "velvetg"), outpath, velvetg.args,
				"-min_contig_lgth", minLength)

	useCutoff <- NULL
	hasCutoffLine <- TRUE
	if ( ! is.null( minCoverage)) {
		useCutoff <- minCoverage
		if ( ! is.na( bestCutoff)) useCutoff <- min( useCutoff, bestCutoff)
	}
	if ( is.null( useCutoff) || is.null( minCoverage) || useCutoff < minCoverage) {
		cmdline <- paste( cmdline," -cov_cutoff auto ")
	} else {
		cmdline <- paste( cmdline," -cov_cutoff", minCoverage)
		hasCutoffLine <- FALSE
	}
	cat( "\nCalling VelvetG to extract 'good' contigs...\n")
	if (verbose) cat( "\nCommand Line: ", cmdline, "\n")
	systemAns <- catch.system( cmdline, intern=doInternal)
	if ( doInternal) {
		estCoverLine <- grep( "Estimated Coverage =", systemAns, value=T)[1]
		cutCoverLine <- grep( "Estimated Coverage cutoff =", systemAns, value=T)[1]
		finalGraphLine <- grep( "Final graph has", systemAns, value=T)[1]
		metricsText <- c( estCoverLine, cutCoverLine, finalGraphLine)
		metricsText[ is.na( metricsText)] <- ""
		print( metricsText)
		velvetMetrics <- extractVelvetgMetrics( metricsText)
	}
	if ( !hasCutoffLine && !is.null(velvetMetrics)) velvetMetrics$CoverCutoff[1] <- bestCutoff
	cat( " Velvet done..  (Kmer=", kmerSize, ")\n", sep="")  

	return( velvetMetrics)
}


`extractVelvetgMetrics` <- function( last3VelvetgTextLines) {

	last3lines <- last3VelvetgTextLines
	estCover <- sub( "(.+Coverage = )(.+)", "\\2", last3lines[1])
	estCut <- sub( "(.+Coverage cutoff = )(.+)", "\\2", last3lines[2])
	if ( last3lines[2] == "") estCut <- NA
	nodes <- sub( "(.+graph has )([0-9]+)( nodes and n50.+)", "\\2", last3lines[3])
	n50 <- sub( "(.+nodes and n50 of )([0-9]+)(.+)", "\\2", last3lines[3])
	maxLen <- sub( "(.+, max )([0-9]+)(.+)", "\\2", last3lines[3])
	totalLen <- sub( "(.+, total )([0-9]+)(.+)", "\\2", last3lines[3])
	using <- sub( "(.+, using )([0-9]+)(/.+)", "\\2", last3lines[3])
	nreads <- sub( "(.+, using )([0-9]+)(/)([0-9]+)( reads)", "\\4", last3lines[3])
	out <- data.frame( "CoverDepth"=as.numeric(estCover), "CoverCutoff"=as.numeric(estCut), "Nodes"=as.numeric(nodes), 
			"N50"=as.numeric(n50), "Max"=as.numeric(maxLen), "Total"=as.numeric(totalLen), 
			"Usage"=(as.numeric(using)*100/as.numeric(nreads)), stringsAsFactors=FALSE)
	return( out)
}


`graphVelvetStats` <- function( path, kmerSize, n50, replot=TRUE, col=1, lwd=1, label="") {

	statsFile <- file.path( path, "stats.txt")
	if ( ! file.exists( statsFile)) {
		cat( "\nNo 'stats.txt' file for Kmer = ", kmerSize)
		return()
	}
	tbl <- read.delim( statsFile, as.is=T)
	N <- nrow(tbl)
	if ( N < 10) {
		cat( "\nNot enough contigs for Kmer = ", kmerSize)
		return()
	}
	len <- tbl$lgth + kmerSize - 1
	hgt <- tbl$short1_cov
	# drop the contigs that never grew at all!!
	keep <- which( !is.na(len) & len > 0)
	allLengths <- len[keep]
	dens <- density( allLengths, adjust=( max(allLengths) ^ 0.125))
	dens$y <- sqrt( dens$y)
	dens$y <- smooth.spline( dens$y, cv=F)$y
	#dens$y <- smooth.spline( dens$y, cv=F)$y

	# know where the n50 lands
	n50 <- as.numeric( n50)
	if ( is.na(n50)) n50 <- 100
	y50 <- dens$y[ which( dens$x >= n50)[1]]
	xlim <- c( 0, n50*1.5)
	ylim <- c( 0, y50 * 6.0)

	if (replot) {
		plot( dens, col=col, lwd=lwd, main=paste( "Velvet 'K-mer' Survey:    Contig Length Distribution:   ", 
				label), xlim=xlim, xlab="Contig Length", ylab="Density:   smoothed SQRT(y)", 
				ylim=ylim, font.axis=2, font.lab=2)
	} else {
		lines( dens, col=col, lwd=lwd)
	}

	#yuse <- min( ymax, par("usr")[4]*0.9)
	#text( xmax, yuse, paste( "Kmer = ", kmerSize), pos=3, col=col, font=2)

	lines( c(n50, n50), c(0,y50), lty=3, lwd=lwd+1, col=col)
	text( n50, y50*1.1, paste( "N50(",kmerSize,") = ", n50, sep=""), pos=3, col=col, font=2, cex=1.2)

	return( invisible( tbl))
}


`makeVelvetPeptides` <- function( sampleID, outpath, verbose=FALSE) {

	require( Biostrings)

	# grab the final set of contigs from the Velvet run
	contigFile <- file.path( outpath, "contigs.fa")
	if ( ! file.exists( contigFile)) {
		cat( "\nVelvet 'contigs.fa' file not found:  ", contigFile)
		return(0)
	}
	contigs <- loadFasta( contigFile, verbose=verbose)
	cat( "\nMaking Peptides from contigs..\nN_Contigs: ", length( contigs$desc))
	if ( length( contigs$desc) < 1) return(0)

	# convert to peptides
	#peps <- DNAtoBestPeptide( contigs$seq, clipAtStop=FALSE)
	peps <- unlist( multicore.lapply( contigs$seq, FUN=DNAtoBestPeptide, clipAtStop=FALSE,
			preschedule=TRUE))
	faOut <- as.Fasta( paste( sampleID, contigs$desc, sep="_"), peps)

	outfile <- paste( sampleID, "velvetPeptides.fasta", sep=".")
	outfile <- file.path( outpath, outfile)
	writeFasta( faOut, file=outfile, line.width=80) 
	if (verbose) cat( "\nWrote Peptides file: ", basename(outfile))

	ncharPeps <- nchar( peps)
	meanPepLen <- mean( ncharPeps)
	pcts <- quantile( ncharPeps, seq(0,1,0.1))
	cat( "\nPeptile Lengths, by quantile:\n")
	print( pcts)
	cat( "\nMean Peptide Length:  ", meanPepLen)
	cat( "\nTotal Amino Acids:    ", prettyNum( sum( ncharPeps), big.mark=","),"\n")

	return( meanPepLen)
}


`bestVelvetProteins` <- function( sampleID, outpath, proteinFastaFile, verbose=TRUE) {

	require( Biostrings)
	data( PAM70)

	# grab the set of peptides from the Velvet run
	peptideFile <- file.path( outpath, paste( sampleID, "velvetPeptides.fasta", sep="."))
	if ( ! file.exists( peptideFile)) {
		cat( "\nVelvet 'peptides' file not found:  ", peptideFile)
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
	cat( "\nFinding best Proteins:  ", basename(peptideFile),
		"\nN_Peptides: ", N, "\n")

	# map peptides to proteins
	ans <- peptide2BestProtein( seqs, proteinFastaFile, substitutionMatrix=PAM70,
				details=FALSE, verbose=verbose)
	if ( is.null(ans)) return(NULL)


	# unwrap the list output
	gid <- score <- scorepaa <- vector( length=N)

	# should be a list of lists, but if only one peptide, its not...
	if (N == 1) {
		tmp <- ans
		ans <- vector( mode="list")
		ans[[1]] <- tmp
	}

	for ( i in 1:N) {
		x <- ans[[i]]
		if ( is.null(x)) next
		gid[i] <- x$Protein
		score[i] <- x$Score
		scorepaa[i] <- x$ScorePerAA
	}

	out <- data.frame( "ContigID"=peptides$desc, "ProteinID"=gid, "Score"=score, 
			"ScorePerAA"=scorepaa, "Sequence"=peptides$seq, stringsAsFactors=F)
	drops <- which( is.na( out$Score))
	if ( length(drops)) out <- out[ -drops, ]

	if ( nrow( out)) {
		ord <- order( out$Score, decreasing=T)
		out <- out[ ord, ]
		rownames(out) <- 1:nrow(out)
	}

	outfile <- paste( sampleID, "velvetProteins.txt", sep=".")
	outfile <- file.path( outpath, outfile)
	write.table( out, file=outfile, sep="\t", quote=F, row.names=F)
	if (verbose) cat( "\nWrote Proteins file: ", basename(outfile))

	return( invisible( out))
}


`bestVelvetVar2csaDomains` <- function( sampleID, outpath, minScorePerAA=1, verbose=TRUE) {

	# grab the set of peptides from the Velvet run
	peptideFile <- file.path( outpath, paste( sampleID, "velvetPeptides.fasta", sep="."))
	if ( ! file.exists( peptideFile)) {
		cat( "\nVelvet 'peptides' file not found:  ", peptideFile)
		return(NULL)
	}

	peptides <- loadFasta( peptideFile, verbose=verbose)
	desc <- peptides$desc
	seqs <- peptides$seq
	N <- length( desc)
	if ( N < 1) return(NULL)

	cat( "\nFinding best Var2csa Domains:  ", "\nN_Peptides: ", N, "\n")

	# search peptides for VAR domains
	ans <- lapply( seqs, findVar2csaDomains, minScorePerAA=minScorePerAA)

	# unwrap the list output
	out <- data.frame()
	nPepHits <- 0
	for ( i in 1:N) {
		x <- ans[[i]]
		if ( is.null(x)) next
		if ( ! nrow(x)) next
		nPepHits <- nPepHits + 1
		sml <- data.frame( "CONTIG"=desc[i], x, stringsAsFactors=F)
		out <- rbind( out, sml)
	}

	if ( nrow(out)) {
		ord <- order( out$SCORE_PER_AA, decreasing=T)
		out <- out[ ord, ]
		rownames(out) <- 1:nrow(out)
	}

	outfile <- paste( sampleID, "velvetVar2csaDomains.txt", sep=".")
	outfile <- file.path( outpath, outfile)
	write.table( out, file=outfile, sep="\t", quote=F, row.names=F)
	cat( "\nN_Var2csa Nodes:   ", nPepHits, "\nN_Var2csa Domains: ", nrow(out), "\n")
	return( out)
}
