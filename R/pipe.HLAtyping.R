# pipe.HLAtyping.R

# turn RNA-seq data into a profile of HLA calls

`pipe.HLAtyping` <- function( sampleID=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt",
				results.path=NULL, doHLAtyping=TRUE, fastqMode=c("RawReads","HLAlocus"),
				python=Sys.which("python"), seq2HLA.path="~/seq2HLA", 
				verbose=TRUE) {

	# path for all results
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	}
	HLAresults.path <- file.path( results.path, "HLA.typing", sampleID)

	# list of Hs_grc HLA genes to harvest
	HLAgenes <- c( "HLA-F:GI3134:06:29723340", "HLA-F-AS1:GI285830:06:29726601", "HLA-V:GI352962:06:29791753", 
			"HLA-P:GI352963:06:29800044", "HLA-G:GI3135:06:29826979", "HLA-H:GI3136:06:29887760", 
			"HLA-T:GI352964:06:29896443", "HLA-K:GI3138:06:29926659", "HLA-U:GI352965:06:29933764", 
			"HLA-A:GI3105:06:29942470", "HLA-W:GI352966:06:29955834", "HLA-J:GI3137:06:30005971", 
			"HLA-L:GI3139:06:30259562", "HLA-N:GI267014:06:30351074", "HLA-E:GI3133:06:30489406", 
			"HLA-C:GI3107:06:31268749", "HLA-B:GI3106:06:31353868", "HLA-S:GI267015:06:31381569", 
			"HLA-X:GI267016:06:31461846", "HLA-DRA:GI3122:06:32439842", "HLA-DRB9:GI3132:06:32459820", 
			"HLA-DRB5:GI3127:06:32517374", "HLA-DRB6:GI3128:06:32552713", "HLA-DRB1:GI3123:06:32578769", 
			"HLA-DQA1:GI3117:06:32637396", "HLA-DQB1:GI3119:06:32659464", "HLA-DQA2:GI3118:06:32741386", 
			"HLA-DQB2:GI3120:06:32756094", "HLA-DOB:GI3112:06:32812763", "HLA-Z:GI267017:06:32896402", 
			"HLA-DMB:GI3109:06:32934629", "HLA-DMA:GI3108:06:32948614", "HLA-DOA:GI3111:06:33004182", 
			"HLA-DPA1:GI3113:06:33064569", "HLA-DPB1:GI3115:06:33075926", "HLA-DPA2:GI646702:06:33091482", 
			"HLA-DPB2:GI3116:06:33112516", "HLA-DPA3:GI267013:06:33131197")

	# use either just the HLA locus reads, or the raw FASTQ reads
	fastqMode <- match.arg( fastqMode)
	didLocusGather <- FALSE

	if ( fastqMode == "HLAlocus") {
		HLAlocusFile <- paste( sampleID, "HLAlocus", "fastq.gz", sep=".")
		HLAlocusFile <- file.path( results.path, "fastq", HLAlocusFile)
		if ( ! file.exists(HLAlocusFile)) {
			# this is only valid for human for now...
			if (getCurrentSpecies() != "Hs_grc") stop( "HLA locus gather only valid for 'Hs_grc' species")
			# step 1:   gather all aligned reads that land in any HLA gene locus
			cat( "\nGathering all HLA aligned reads..\n")
			pipe.GatherGeneAlignments( sampleID, HLAgenes, asFASTQ=TRUE, fastq.keyword="HLAlocus")
			didLocusGather <- TRUE
		}
	} else {
		fastq.path <- getOptionValue( optionsFile, "fastqData.path", notfound=".", verbose=F)
		isPaired <- getAnnotationTrue( annotationFile, sampleID, "PairedEnd", notfound=FALSE)
		if ( ! isPaired) {
			cat( "\nHLA typing from raw reads needs paired end fastq files.")
			cat( "\nTry 'HLAlocus' mode which treats reads as unpaired.")
			cat( "\nNo HLA typing performed..")
			return(NULL)
		}
		fastqFiles <- getAnnotationValue( annotationFile, sampleID, "Filename")
		fastqFiles <- strsplit( fastqFiles, split=", *")[[1]]
		if ( length(fastqFiles) != 2) stop( "Expected two comma separated Fastq filenames for this sampleID")
		fastqFiles <- file.path( fastq.path, fastqFiles)
	}

	# step 2:  Throw those reads at seq2HLA
	if (didLocusGather || doHLAtyping) {

		# pre-test the Python installation and Bio packages
		cat( "\n\nTesting python install & 'Bio' package..")
		testCmd <- paste( python, ' 2>&1  -c "from Bio import SeqIO"')
		ans <- system( testCmd, intern=TRUE)
		if ( length(ans)) {
			cat( "\nPython Test Failed..   Traceback:\n")
			cat( ans, sep="\n")
			stop( "Unable to use Python to call 'seq2HLA'")
		}

		cat( "\n\nCalling 'seq2HLA'..")
		seq2HLA.program <- if (fastqMode == "HLAlocus") "seqUnpaired2HLA.py" else "seq2HLA.py"
		hlaProgram <- file.path( seq2HLA.path, seq2HLA.program)
		outFilePrefix <- file.path( HLAresults.path, sampleID)

		if (fastqMode == "HLAlocus") {
			cmdline <- paste( python, " ", hlaProgram, "  -u", HLAlocusFile, " -r", outFilePrefix, sep=" ")
		} else {
			cmdline <- paste( python, " ", hlaProgram, "  -1", fastqFiles[1], " -2", fastqFiles[2], 
					" -r", outFilePrefix, sep=" ")
		}
		if (verbose) cat( "\nHLA typing command line:\n", cmdline)

		system( cmdline)
	}

	# there are several files of text results...  grab them and combine...
	cat( "\n\nSummarizing 'seq2HLA' result files..")
	fPrefix <- paste( sampleID, "-Class", sep="")

	# there are 2 files of expression info
	type1 <- readLines( file.path( HLAresults.path, paste( fPrefix, "I.expression.txt", sep="")))
	terms <- strsplit( type1, split=" ", fixed=T)
	type1type <- sapply( terms, "[[", 1)
	type1rpkm <- sapply( terms, "[[", 2)
	type2 <- readLines( file.path( HLAresults.path, paste( fPrefix, "II.expression.txt", sep="")))
	terms <- strsplit( type2, split=" ", fixed=T)
	type2type <- sapply( terms, "[[", 1)
	type2rpkm <- sapply( terms, "[[", 2)

	locus <- c( type1type, type2type) 
	rpkm <- round( as.numeric( c( type1rpkm, type2rpkm)), digits=2) 
	out <- data.frame( "HLA_Locus"=locus, "RPKM"=rpkm, stringsAsFactors=F)
	ord <- order( out$HLA_Locus)
	out <- out[ ord, ]
	outfile <- file.path( HLAresults.path, paste( sampleID, "HLA.Expression.txt", sep="."))
	write.table( out, outfile, sep="\t", quote=F, row.names=F)
	cat( "\nWrote file:  ", outfile)

	# there are 4 files of genotype calls
	type1.2 <- readLines( file.path( HLAresults.path, paste( fPrefix, "I.HLAgenotype2digits.txt", sep="")))
	type1.2 <- type1.2[ -grep( "^#", type1.2)]
	terms <- strsplit( type1.2, split="\t", fixed=T)
	type1.2type <- paste( sapply( terms, "[[", 1), "2digit", sep=":")
	type1.2allele1 <- sapply( terms, "[[", 2)
	type1.2pval1 <- sapply( terms, "[[", 3)
	type1.2allele2 <- sapply( terms, "[[", 4)
	type1.2pval2 <- sapply( terms, "[[", 5)

	type1.4 <- readLines( file.path( HLAresults.path, paste( fPrefix, "I.HLAgenotype4digits.txt", sep="")))
	type1.4 <- type1.4[ -grep( "^#", type1.4)]
	terms <- strsplit( type1.4, split="\t", fixed=T)
	type1.4type <- paste( sapply( terms, "[[", 1), "4digit", sep=":")
	type1.4allele1 <- sapply( terms, "[[", 2)
	type1.4pval1 <- sapply( terms, "[[", 3)
	type1.4allele2 <- sapply( terms, "[[", 4)
	type1.4pval2 <- sapply( terms, "[[", 5)

	type2.2 <- readLines( file.path( HLAresults.path, paste( fPrefix, "II.HLAgenotype2digits.txt", sep="")))
	type2.2 <- type2.2[ -grep( "^#", type2.2)]
	terms <- strsplit( type2.2, split="\t", fixed=T)
	type2.2type <- paste( sapply( terms, "[[", 1), "2digit", sep=":")
	type2.2allele1 <- sapply( terms, "[[", 2)
	type2.2pval1 <- sapply( terms, "[[", 3)
	type2.2allele2 <- sapply( terms, "[[", 4)
	type2.2pval2 <- sapply( terms, "[[", 5)

	type2.4 <- readLines( file.path( HLAresults.path, paste( fPrefix, "II.HLAgenotype4digits.txt", sep="")))
	type2.4 <- type2.4[ -grep( "^#", type2.4)]
	terms <- strsplit( type2.4, split="\t", fixed=T)
	type2.4type <- paste( sapply( terms, "[[", 1), "4digit", sep=":")
	type2.4allele1 <- sapply( terms, "[[", 2)
	type2.4pval1 <- sapply( terms, "[[", 3)
	type2.4allele2 <- sapply( terms, "[[", 4)
	type2.4pval2 <- sapply( terms, "[[", 5)

	locus <- c( type1.2type, type1.4type, type2.2type, type2.4type) 
	allele1 <- c( type1.2allele1, type1.4allele1, type2.2allele1, type2.4allele1) 
	pval1 <- round( as.numeric( c( type1.2pval1, type1.4pval1, type2.2pval1, type2.4pval1)), digits=6)
	allele2 <- c( type1.2allele2, type1.4allele2, type2.2allele2, type2.4allele2) 
	pval2 <- round( as.numeric( c( type1.2pval2, type1.4pval2, type2.2pval2, type2.4pval2)), digits=6)

	out <- data.frame( "HLA_Locus"=locus, "Allele_1"=allele1, "Pvalue_1"=pval1,
			"Allele_2"=allele2, "Pvalue_2"=pval2, stringsAsFactors=F)
	ord <- order( out$HLA_Locus)
	out <- out[ ord, ]
	outfile <- file.path( HLAresults.path, paste( sampleID, "HLA.Genotype.txt", sep="."))
	write.table( out, outfile, sep="\t", quote=F, row.names=F)
	cat( "\nWrote file:  ", outfile)
	return()
}


`pipe.HLAsummary` <- function( sampleIDs, annotationFile="Annotation.txt", optionsFile="Options.txt",
				results.path=NULL) {

	# path for all results
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	}

	# gather up the HLA data for these samples
	outID <- vector()
	outHLA <- matrix( "", nrow=length(sampleIDs), ncol=6)
	colnames(outHLA) <- c( "A", "B", "C", "DQA", "DQB", "DRB")

	for ( s in sampleIDs) {
		HLAresults.path <- file.path( results.path, "HLA.typing", s)
		f <- file.path( HLAresults.path, paste( s, "HLA.Genotype.txt", sep="."))
		if ( ! file.exists( f)) {
			cat( "\nHLA results file not found: ", f)
			next
		}

		tbl <- read.delim( f, as.is=T)
		is2digit <- grep( "2digit", tbl$HLA_Locus)
		is4digit <- grep( "4digit", tbl$HLA_Locus)
		NHLA <- length( is2digit)
		str <- rep.int( "", NHLA)
		for ( j in 1:NHLA) {
			all1.2 <- tbl$Allele_1[is2digit[j]]
			all2.2 <- tbl$Allele_2[is2digit[j]]
			all1.4 <- tbl$Allele_1[is4digit[j]]
			all2.4 <- tbl$Allele_2[is4digit[j]]
			if ( all1.4 != "no") all1.2 <- all1.4
			if ( all2.4 != "no") all2.2 <- all2.4
			str[j] <- paste( all1.2, all2.2, sep="/")
		}
		outID <- c( outID, s)
		outHLA[ length(outID), ] <- str
	}

	out <- data.frame( "SampleID"=outID, outHLA[ 1:length(outID), ], stringsAsFactors=F)
	rownames(out) <- 1:nrow(out)
	out
}



`pipe.HLA.byGroup` <- function( sampleIDs, annotationFile="Annotation.txt", optionsFile="Options.txt",
				results.path=NULL, groupColumn="Group", colorColumn="Color",
				mainText="", N_SIMUL=1000) {

	# grab the table of HLA summary calls
	hlaIn <- pipe.HLAsummary( sampleIDs, annotationFile=annotationFile, optionsFile=optionsFile,
				results.path=results.path)

	# gather up the HLA data for these samples
	outID <- vector()
	outHLA <- matrix( "", nrow=length(sampleIDs), ncol=6)

	# the annotation holds the group calls
	annT <- readAnnotationTable( annotationFile)
	allGroups <- annT[[ groupColumn]]
	myGroup <- checkGroupNames( allGroups[ match( sampleIDs, annT$SampleID)])
	groupNames <- sort( unique( myGroup))
	nGroups <- length( groupNames)
	groupColors <- annT[[ colorColumn]][ match( groupNames, allGroups)]

	# clean the HLA data
	# make a 2x sized table, to keep both of the 2 calls per subject
	hlaOut <- hlaIn
	hlaTerms <- rbind( hlaIn, hlaIn)
	N_In <- nrow( hlaIn)
	HLA_Groups <- c( "A", "B", "C", "DQA", "DQB", "DRB")
	for ( HLAgrp in HLA_Groups) {
		tmp <- hlaIn[, HLAgrp]
		# some have "hom(", X, ")" motifs
		tmp <- sub( "(hoz\\()(.+)(\\)$)", "\\2", tmp)
		# put the 2 calls into sorted order
		terms <- strsplit( tmp, split="/", fixed=T)
		terms <- lapply( terms, function(x) sort( x))
		term1 <- sapply( terms, `[`, 1)
		term2 <- sapply( terms, `[`, 2)
		tmp <- sapply( terms, function(x) paste( x, collapse="/"))
		hlaOut[, HLAgrp] <- tmp
		hlaTerms[ 1:N_In, HLAgrp] <- term1
		hlaTerms[ (N_In+1):(N_In*2), HLAgrp] <- term2
	}

	# OK, set the group for each row
	where <- match( hlaTerms$SampleID, sampleIDs)
	hlaTerms$Group <- myGroup[ where]
	
	# set up to do a simulated FDR, of the Infect/Protect calls, and perhaps the HLA calls themselves too
	SIMUL_P <- TRUE
	SIMUL_HLA <- TRUE
	
	checkX11()
	padFactor <- 1.6
	par( mai=c(1,1,0.8,0.4))

	for ( HLAgrp in HLA_Groups) {
		allTbl <- table( hlaTerms[ , HLAgrp])
		labels <- names( allTbl)
		# get the percentages to many digits for later simulations
		allPcts <- round( allTbl * 100 / sum(allTbl), digits=5)
		# make a matrix to hold the breakdown by group
		m <- matrix( 0, nrow=length(labels), ncol=nGroups)
		colnames(m) <- groupNames
		rownames(m) <- labels
		for ( igrp in 1:nGroups) {
			who <- which( hlaTerms$Group == groupNames[igrp])
			thisTbl <- table( factor( hlaTerms[ who, HLAgrp], levels=labels))
			thisPcts <- round( thisTbl * 100 / sum(thisTbl), digits=4)
			m[ , igrp] <- thisPcts
		}
	
		barAns <- barplot( t(m), beside=T, col=groupColors, ylim=c(0,max(m)*1.24), xlim=c(0.4,(nrow(m)*ncol(m)*padFactor)),
				main=paste( mainText, ":     HLA Group = ", HLAgrp, sep=""),
				ylab="Percentage of all HLA type Calls", xlab=NA, las=3,
				font.axis=2, font.lab=2, cex.axis=1.1, cex.lab=1.1)
		legend( 'topright', colnames(m), fill=groupColors, bg='white', cex=1.1)
	
		if (SIMUL_P && N_SIMUL > 0) {
			simPcts <- array( 0, dim=c( N_SIMUL, length(labels), nGroups))
			# build random simulated percentages
			cat( "\nSimulating ", HLAgrp, "\n")
			for( k in 1:N_SIMUL) {
				# randomize the group calls for each subject
				randRows <- sample( nrow(hlaTerms))
				randTerms <- hlaTerms[ , HLAgrp]
				if (SIMUL_HLA) randTerms <- sample( labels, size=length(randTerms), replace=TRUE, prob=allPcts)
				for (igrp in 1:nGroups) {
					who <- which( hlaTerms$Group == groupNames[ igrp])
					randWho <- randRows[ who]
					randTbl <- table( factor( randTerms[ randWho], levels=labels))
					randPcts <- round( randTbl * 100 / sum(randTbl), digits=4)
					simPcts[ k, , igrp] <- randPcts
				}
				if( k %% 100 == 0) cat( "\r", k)
			}
			cat( "  Done.")
			# now see how often the random was 'better' than the actual
			for (lvl in 1:nrow(m)) {
				thisValueSet <- m[ lvl, ]
				thisDiff <- diff( range( thisValueSet))
				randValueM <- simPcts[ , lvl, ]
				randDiffs <- apply( randValueM, 1, function(x) diff(range(x)))
				nRandBetter <- sum( randDiffs > thisDiff)
				fdr <- round( nRandBetter / N_SIMUL, digits=3)
				if ( fdr <= 0.2) {
					ptxt <- paste( "P=", fdr, sep="")
					ytxt <- max( m[ lvl, ])
					xtxt <- mean( barAns[ ,lvl])
					# show the P-value either horizontal, or vertical if too many HLA types
					if ( nrow(m) < 9) {
						text( xtxt, ytxt, ptxt, cex=0.85, pos=3)
					} else {
						xAdjust <- (nGroups * 0.4) - 0.1
						text( xtxt-xAdjust, ytxt+0.25, ptxt, cex=0.85, srt=90, pos=4, offset=1)
					}
				}
			}
		}
		dev.flush();  Sys.sleep(1)
		dev.print( png, paste( "HLA.TypeByGroup", HLAgrp, "png", sep="."), width=1000, height=650)
	}
}


cleanupHLAtypingFiles <- function( path="results/HLA.typing") {

	# delete any large unneeded files left behind by Velvet...
	filePatterns <- c( ".aligned$")

	totalBytes <- 0
	nfiles <- 0

	for (patt in filePatterns) {
		files <- dir( path, pattern=patt, recursive=T, full.name=T)
		if ( length( files) < 1) next
		for ( f in files) {
			bytes <- file.info( f)$size
			nfiles <- nfiles + 1
			cat( "\r", nfiles, "  Size: ", bytes, "  File: ", f)
			totalBytes <- totalBytes + bytes
			file.delete( f)
		}
		cat( "\n")
	}
	cat( "\nDeleted_Files: ", nfiles, "\tDeleted Bytes: ", prettyNum( totalBytes, big.mark=","), "\n")
}
