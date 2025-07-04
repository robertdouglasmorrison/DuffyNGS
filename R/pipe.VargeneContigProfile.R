# pipe.VargeneContigProfile.R

# turn RNA-seq data into a profile of Vargene constructs & expression

# giving 2 tools to do the task:  Velvet or SPAdes


`pipe.VargeneContigProfile.Spades` <- function( sampleID, vargeneFastaFile, annotationFile="Annotation.txt", optionsFile="Options.txt",
				results.path=NULL, spades.path=dirname(Sys.which("spades.py")), doSpades=FALSE,
				spades.mode=c("isolate","rna","meta"), kmerSizes=NULL, spades.args="", 
				doCutadapt=TRUE, cutadaptProgram=Sys.which("cutadapt"), keyword="PfEMP1", 
				min.aa.length=100, min.score.per.aa=2, verbose=TRUE) {

	if ( is.null( results.path)) results.path <- getOptionValue( optionsFile, "results.path", 
				notfound=".", verbose=F)

	# make sure the reference file of var genes is readable
	if ( ! file.exists( vargeneFastaFile)) {
		cat( "\nError: 'vargeneFastaFile' of reference proteins not found: ", vargeneFastaFile)
		return(NULL)
	}

	# path for all results
	spades.output.path <- file.path( results.path, "SpadesContigs", sampleID, "PfEMP1")
	contigsFastaFile <- file.path( spades.output.path, "contigs.fasta")
	peptidesFastaFile <- file.path( spades.output.path, paste( sampleID, keyword, "Peptides.fasta", sep="."))
	domainsFastaFile <- file.path( spades.output.path, paste( sampleID, keyword, "Domains.fasta", sep="."))
	proteinsFastaFile <- sub( "Domains.fasta$", "Proteins.fasta", domainsFastaFile)
	domainsTextFile <- file.path( spades.output.path, paste( sampleID, keyword, "DomainDetails.txt", sep="."))
	proteinsTextFile <- sub( "DomainDetails.txt$", "BestProteinHits.txt", domainsTextFile)
	
	# remove any files we know we will remake
	if (doSpades) file.delete( c( contigsFastaFile, peptidesFastaFile, domainsTextFile, domainsFastaFile,
					proteinsFastaFile, proteinsTextFile))

	# get the universe of PfEMP1 var genes
	vargeneDomainMap <- getVargeneDomainMap()
	vargenes <- sort( unique( vargeneDomainMap$GENE_NAME))
	alignedReadsFile <- file.path( results.path, "fastq", paste( sampleID, "PfEMP1.Hits.fastq.gz", sep="."))
	nohitReadsFile <- file.path( results.path, "fastq", paste( sampleID, "noHits.fastq.gz", sep="."))

	if ( doSpades || !file.exists(contigsFastaFile)) {

		# verify we see the denovo tool executable location
		if ( spades.path == "") {
			cat( "\nError:  path to SPAdes executable program not found..")
			return(NULL)
		}

		# step 1:   gather all aligned reads that land in/near the vargene loci
		if ( ! file.exists( alignedReadsFile)) {
			cat( "\nGathering Var Gene aligned reads..\n")
			pipe.GatherGeneAlignments( sampleID, genes=vargenes, asFASTQ=TRUE, mode="best.one",
						fastq.keyword="PfEMP1.Hits")
		}

		# in most cases call cutadapt to trim off primer ends
		# if we are doing the Cutadapt pass, check for those files, and call it if we need.
		fastqFiles <- c( "PfEMP1.Hits", "noHits")
		if (doCutadapt) {
			inputFastqFiles <- c( alignedReadsFile, nohitReadsFile)
			filesToDo <- checkOrCallCutadapt( inputFastqFiles, asMatePairs=FALSE, forceMatePairs=FALSE,
						cutadaptProgram=cutadaptProgram, kmer.size=45, verbose=verbose)
			fastqFiles <- paste( fastqFiles, "trimmed", sep=".")
		}
	
		# step 2:  create contigs of anything that may be var gene reads
		spades.mode <- match.arg( spades.mode)
		pipe.SpadesContigs( sampleID, fastqSource=fastqFiles, folderName="PfEMP1", 
				spades.mode=spades.mode, kmerSizes=kmerSizes, spades.args=spades.args, 
				spades.path=spades.path, keyword=keyword, makePep=T)

		# cleanup SPAdes temp files
		cleanupSpadesFiles( path=spades.output.path, verbose=F)
	} else {
		cat( "\n\nUsing previously calculated Spades contigs..\n")
	}

	# because of reading frame and stop codon trimming, and that SPAdes has no length cutoff,
	# we may have peptides that are now shorter than the minimun we want.
	pepFA <- loadFasta( peptidesFastaFile, short=T, verbose=F)
	len <- nchar( pepFA$seq)
	drops <- which( len < min.aa.length)
	if ( length( drops)) {
		pepFA <- as.Fasta( pepFA$desc[-drops], pepFA$seq[-drops])
		writeFasta( pepFA, peptidesFastaFile, line=100)
	}
	NPEP <- length( pepFA$desc)

	# step 3:  Throw those contigs as peptides against the given set of full length var gene proteins
	cat( "\n\nSearching Contigs for PfEMP1 constructs..")
	proteinAns <- bestSpadesProteins( sampleID, outpath=spades.output.path, proteinFastaFile=vargeneFastaFile,
				keyword=keyword, verbose=verbose)
	# visit each peptide's best match explicitly, and either drop, keep, and maybe trim
	protDesc <- protSeq <- vector()
	NPROT <- 0
	for ( i in 1:nrow(proteinAns)) {
		thisScore <- proteinAns$ScorePerAA[i]
		if ( thisScore < min.score.per.aa) next
		thisID <- proteinAns$ContigID[i]
		where <- match( thisID, pepFA$desc, nomatch=0)
		if ( ! where) next
		NPROT <- NPROT + 1
		protDesc[NPROT] <- thisID
		protSeq[NPROT] <- pepFA$seq[where]
		# see if we should trim any untranslated UTR from the peptide
		thisProtStart <- proteinAns$ProteinStart[i]
		if ( thisProtStart < 1) protSeq[NPROT] <- substring( protSeq[NPROT], abs(thisProtStart)+2)
	}
	if (NPROT) {
		protFA <- as.Fasta( protDesc, protSeq)
		writeFasta( protFA, proteinsFastaFile, line=100)
		cat( "\nWrote", NPROT, "PfEMP1 proteins as FASTA to: ", basename(proteinsFastaFile))
	}

	# at this point, we have a small chance that none look like var genes
	if ( ! NPROT) {
		cat( "\nNo peptides look like PfEMP1 proteins..")
		proteinAns$Domain.Architecture <- ""
		proteinAns$Cassette.Architecture <- ""
		write.table( proteinAns, proteinsTextFile, sep="\t", quote=F, row.names=F)
		return( invisible( proteinAns))
	}

	# step 4:  find the var gene domains for each
	cat( "\n\nSearching for PfEMP1 Domains..\n")
	domStrs <- cassStrs <- rep.int( "", NPROT)
	domainAns <- multicore.lapply( 1:NPROT, function(x) {
				mySeq <- protFA$seq[x]
				myDesc <- protFA$desc[x]
				domAns <- findVsaDomains( mySeq, keepVSApattern="3D7|DD2|HB3|IGH|IT4")
				if ( ! nrow(domAns)) return(domAns)
				# drop the domain columns we do not need
				domAns <- domAns[ , -grep("^REF_",colnames(domAns))]
				# make the stirngs about the architecture
				myDomStr <- paste( domAns$DOMAIN_ID, collapse="-")
				myCassStr <- paste( domAns$CASSETTE, collapse="-")
				if (verbose) cat( "\r", x, myDesc, myDomStr, "  ")
				# attach the contig name to the domain ans
				domAns <- data.frame( "CONTIG_ID"=myDesc, domAns, 
						"Domain.Architecture"=myDomStr, "Cassette.Architecture"=myCassStr, 
						stringsAsFactors=F)
				return( domAns)
			})

	# organize the results.  Merge all the domain details, and append the archtecture details to the proteins file
	out2 <- data.frame()
	for ( j in 1:length(domainAns)) {
		sml <- domainAns[[j]]
		if ( ! nrow(sml)) next
		# extract the architecture strings and remove them
		domStrs[j] <- sml$Domain.Architecture[1]
		cassStrs[j] <- sml$Cassette.Architecture[1]
		sml <- sml[ , -grep( "Architecture",colnames(sml))]
		out2 <- rbind( out2, sml)
	}
	proteinAns$Domain.Architecture <- ""
	proteinAns$Cassette.Architecture <- ""
	whereProtein <- match( proteinAns$ContigID, protFA$desc, nomatch=0)
	proteinAns$Domain.Architecture[ whereProtein > 0] <- domStrs[ whereProtein]
	proteinAns$Cassette.Architecture[ whereProtein > 0] <- cassStrs[ whereProtein]

	# rewrite the proteins file with the new extra architecture info
	write.table( proteinAns, proteinsTextFile, sep="\t", quote=F, row.names=F)
	# write out the file of domain details, and also write the domains as fasta too
	if ( nrow(out2)) {
		write.table( out2, domainsTextFile, sep="\t", quote=F, row.names=F)
		domDesc <- paste( out2$CONTIG_ID, out2$DOMAIN_ID, sep="_")
		domSeq <- out2$QUERY_SEQ
		domFA <- as.Fasta( domDesc, domSeq)
		writeFasta( domFA, domainsFastaFile, line=100)
	}
	return( invisible( proteinAns))
}


`pipe.VargeneContigSummary.Spades` <- function( sampleIDset, optionsFile="Options.txt", 
						results.path=NULL, keyword="PfEMP1") {

	if ( is.null( results.path)) results.path <- getOptionValue( optionsFile, "results.path", 
				notfound=".", verbose=F)

	# gather up all the vargene contig files for the named samples
	spades.path <- file.path( results.path, "SpadesContigs")
	N <- length( sampleIDset)
	cat( "\nChecking", N, "folders of Spades results..")

	domainTextFiles <- domainFastaFiles <- proteinTextFiles <- proteinFastaFiles <- vector()
	for ( i in 1:N) {
		thisSID <- sampleIDset[i]
		this.path <- file.path( spades.path, thisSID, keyword)
		if ( ! file.exists(this.path)) {
			cat( "\nWarn: no Spades vargene results folder found for:", thisSID)
			next
		}
		domainTextFiles[i] <- file.path( this.path, paste( thisSID, keyword, "DomainDetails.txt", sep="."))
		proteinTextFiles[i] <- file.path( this.path, paste( thisSID, keyword, "BestProteinHits.txt", sep="."))
		domainFastaFiles[i] <- file.path( this.path, paste( thisSID, keyword, "Domains.fasta", sep="."))
		proteinFastaFiles[i] <- file.path( this.path, paste( thisSID, keyword, "Proteins.fasta", sep="."))
	}

	cat( "\nGathering Best Protein Hits..\n")
	protDF <- data.frame()
	proteinFile <- "VargeneContigs.PfEMP1.BestProteinHits.txt"
	for ( f in proteinTextFiles) {
		if ( is.na(f)) next
		if ( ! file.exists(f)) next
		sml <- read.delim( f, as.is=T)
		cat( "  ", basename(f), "=", nrow(sml),sep="")
		if ( ! nrow(sml)) next
		protDF <- rbind( protDF, sml)
	}
	write.table( protDF, proteinFile, sep="\t", quote=F, row.names=F)
	cat( "\nWrote: ", proteinFile, " \tN_Protein_Hits: ", nrow(protDF))

	cat( "\nGathering Domain Details..\n")
	domDF <- data.frame()
	domainFile <- "VargeneContigs.PfEMP1.DomainDetails.txt"
	for ( f in domainTextFiles) {
		if ( is.na(f)) next
		if ( ! file.exists(f)) next
		sml <- read.delim( f, as.is=T)
		cat( "  ", basename(f), "=", nrow(sml),sep="")
		if ( ! nrow(sml)) next
		domDF <- rbind( domDF, sml)
	}
	write.table( domDF, domainFile, sep="\t", quote=F, row.names=F)
	cat( "\nWrote: ", domainFile, " \tN_Domain_Hits: ", nrow(domDF))

	cat( "\nGathering Protein FASTA..\n")
	desc <- seq <- vector()
	proteinFile <- "VargeneContigs.PfEMP1.Proteins.fasta"
	for ( f in proteinFastaFiles) {
		if ( is.na(f)) next
		if ( ! file.exists(f)) next
		sml <- loadFasta( f, short=F, verbose=F)
		cat( "  ", basename(f), "=", length(sml$desc),sep="")
		if ( ! length(sml$desc)) next
		desc <- c( desc, sml$desc)
		seq <- c( seq, sml$seq)
	}
	writeFasta( as.Fasta( desc, seq), proteinFile, line=100)
	cat( "\nWrote: ", proteinFile, " \tN_Proteins: ", length(desc))

	cat( "\nGathering Domain FASTA..\n")
	desc <- seq <- vector()
	domainFile <- "VargeneContigs.PfEMP1.Domains.fasta"
	for ( f in domainFastaFiles) {
		if ( is.na(f)) next
		if ( ! file.exists(f)) next
		sml <- loadFasta( f, short=F, verbose=F)
		cat( "  ", basename(f), "=", length(sml$desc),sep="")
		if ( ! length(sml$desc)) next
		desc <- c( desc, sml$desc)
		seq <- c( seq, sml$seq)
	}
	writeFasta( as.Fasta( desc, seq), domainFile, line=100)
	cat( "\nWrote: ", domainFile, " \tN_Domains: ", length(desc))

	return( N)
}


`pipe.VargeneContigProfile.Velvet` <- function( sampleID, vargeneFastaFile, annotationFile="Annotation.txt", optionsFile="Options.txt",
				results.path=NULL, kmerSize=55, doVelvet=FALSE, velvet.path=dirname(Sys.which("velveth")), 
				doCutadapt=TRUE, cutadaptProgram=Sys.which("cutadapt"), keyword="PfEMP1", 
				min.aa.length=200, min.score.per.aa=2, minCoverage=3, verbose=TRUE) {

	if ( is.null( results.path)) results.path <- getOptionValue( optionsFile, "results.path", 
				notfound=".", verbose=F)

	# verify we see the denovo tool executable location
	if ( velvet.path == "") {
		cat( "\nError:  path to Velvet executable program not found..")
		return(NULL)
	}

	# make sure the reference file of var genes is readable
	if ( ! file.exists( vargeneFastaFile)) {
		cat( "\nError: 'vargeneFastaFile' of reference proteins not found: ", vargeneFastaFile)
		return(NULL)
	}

	# path for all results
	velvet.output.path <- file.path( results.path, "VelvetContigs", sampleID, "PfEMP1")
	contigsFastaFile <- file.path( velvet.output.path, "contigs.fa")
	peptidesFastaFile <- file.path( velvet.output.path, paste( sampleID, keyword, "Peptides.fasta", sep="."))
	domainsTextFile <- file.path( velvet.output.path, paste( sampleID, keyword, "DomainDetails.txt", sep="."))
	
	# remove any files we know we will remake
	if (doVelvet) file.delete( c( contigsFastaFile, peptidesFastaFile, domainsTextFile))

	# get the universe of PfEMP1 var genes
	vargeneDomainMap <- getVargeneDomainMap()
	vargenes <- sort( unique( vargeneDomainMap$GENE_NAME))
	alignedReadsFile <- file.path( results.path, "fastq", paste( sampleID, "PfEMP1.Hits.fastq.gz", sep="."))
	nohitReadsFile <- file.path( results.path, "fastq", paste( sampleID, "noHits.fastq.gz", sep="."))

	# note the the kmerSize can be a triplet or a single value
	kmerSizeCutadapt <- max( kmerSize, na.rm=T)

	min.dna.length <- (min.aa.length - 1) * 3
	if ( doVelvet || !file.exists(contigsFastaFile)) {

		# step 1:   gather all aligned reads that land in/near the vargene loci
		if ( ! file.exists( alignedReadsFile)) {
			cat( "\nGathering Var Gene aligned reads..\n")
			pipe.GatherGeneAlignments( sampleID, genes=vargenes, asFASTQ=TRUE, mode="best.one",
						fastq.keyword="PfEMP1.Hits")
		}

		# in most cases call cutadapt to trim off primer ends
		# if we are doing the Cutadapt pass, check for those files, and call it if we need.
		fastqFiles <- c( "PfEMP1.Hits", "noHits")
		if (doCutadapt) {
			inputFastqFiles <- c( alignedReadsFile, nohitReadsFile)
			filesToDo <- checkOrCallCutadapt( inputFastqFiles, asMatePairs=FALSE, forceMatePairs=FALSE,
						cutadaptProgram=cutadaptProgram, kmer.size=kmerSizeCutadapt, verbose=verbose)
			fastqFiles <- paste( fastqFiles, "trimmed", sep=".")
		}
	
		# step 2:  create contigs of anything that may be var gene reads
		pipe.VelvetContigs( sampleID, fastqSource=fastqFiles, kmerSize=kmerSize, folderName="PfEMP1", keyword=keyword,
					velvet.path=velvet.path, minCoverage=minCoverage, minLength=min.dna.length, makePep=T)

		# cleanup Velvet temp files
		tempFiles <- file.path( velvet.output.path, c( "Roadmaps", "Sequences", "Graph", "Graph2", "PreGraph", "LastGraph"))
		file.delete( tempFiles)
	}

	# because of reading frame and stop codon trimming, we may have peptides that are now shorter than the minimun we want.
	pepFA <- loadFasta( peptidesFastaFile, short=F)
	len <- nchar( pepFA$seq)
	drops <- which( len < min.aa.length)
	if ( length( drops)) {
		pepFA <- as.Fasta( pepFA$desc[-drops], pepFA$seq[-drops])
		writeFasta( pepFA, peptidesFastaFile, line=100)
	}
	NPEP <- length( pepFA$desc)

	# step 3:  Throw those contigs as peptides against the given set of proteins
	cat( "\n\nSearching Contigs for PfEMP1 constructs..")
	ans <- bestVelvetProteins( sampleID, outpath=velvet.output.path, proteinFastaFile=vargeneFastaFile,
				keyword=keyword, verbose=verbose)
	drops <- which( ans$ScorePerAA < min.score.per.aa)
	if ( length( drops)) ans <- ans[ -drops, ]

	# step 4:  find the var gene domains for each
	cat( "\n\nFinding PfEMP1 Domains..\n")
	domainString <- unlist( multicore.lapply( 1:NPEP, function(x) {
				i <- x
				mySeq <- pepFA$seq[i]
				myDesc <- pepFA$desc[i]
				domAns <- findVsaDomains( mySeq, keepVSApattern="3D7|DD2|HB3|IGH|IT4")
				myDomStr <- if ( nrow(domAns)) paste( domAns$DOMAIN_ID, collapse="-") else "none"
				cat( "\r", i, myDesc, myDomStr, "  ")
				return( myDomStr)
			}))

	# place these where they belong in the final output, and reformat a bit
	# note the the output may be less rows, due to low ScorePerAA, so match accordingly
	where <- match( ans$ContigID, pepFA$desc)
	ansDomStr <- domainString[ where]
	out <- cbind( ans[,1:2], "Domain.Architecture"=ansDomStr, ans[,3:ncol(ans)], stringsAsFactors=F)
	write.table( out, domainsTextFile, sep="\t", quote=F, row.names=F)
	return( invisible( out))
}


`genePlotWithContigs` <- function( sampleID, geneID, tail=2500, pileup.col=2, text.col=3, cex=1.65, pos=NULL, 
				results.path="results", ...) {

	pipe.PlotGene( sampleID, geneID, tail=tail, col=pileup.col, label=paste( "  Sample:  ",sampleID))
	overlayContigs( sampleID, seqID=subset( getCurrentGeneMap(), GENE_ID==geneID)$SEQ_ID, 
			results.path=results.path, col=text.col, cex=cex, pos=pos)
}


`overlayContigs` <- function( sampleID, seqID=getCurrentSeqMap()$SEQ_ID[1], results.path="results",
		contigSubfolder="Vargenes", cex=1.0, col=3, pos=NULL, min.length=200) {

	velvetPath <- file.path( results.path, "VelvetContigs", sampleID, contigSubfolder)
	contigsFastaFile <- file.path( velvetPath, "contigs.fa")
	bestHitFile <- file.path( velvetPath, "Pf.BestBlastHit.csv")
	hits <- read.csv( bestHitFile, as.is=T)

	# don't show any that are too short to be useful
	# there is no one value to know the contig length, so deduce as best you can
	len1 <- hits$LEN_MATCH
	len2 <- (hits$P_LAST - hits$P_FIRST + 1)
	len3 <- (hits$S_BEG - hits$S_END + 1)
	len <- pmax( len1, len2, len3)
	drops <- which( len < min.length)
	if ( length(drops) > 0) hits <- hits[ -drops, ]

	plotLimits <- par("usr")
	loBase <- as.integer( plotLimits[1])
	hiBase <- as.integer( plotLimits[2])
	loY <- plotLimits[3]
	hiY <- plotLimits[4]
	NX <- hiBase - loBase + 1

	# see which contigs touch this plot window
	who <- which( hits$SEQ_ID == seqID & hits$S_BEG < hiBase & hits$S_END > loBase)
	if ( length( who) < 1) return()

	# visit them in  chromosomal order...
	hits <- hits[ who, ]
	ord <- order( hits$S_BEG)
	hits <- hits[ ord, ]
	who <- 1:nrow(hits)

	# got at least one, so load the contigs
	fa <- loadFasta( contigsFastaFile)

	# keep track of how 'deep' we have contigs on the plot
	topdepth <- rep.int( 0, NX)
	botdepth <- rep.int( hiY, NX)
	names(topdepth) <- names(botdepth)  <- loBase:hiBase

	# visit each, and draw it on the plot
	for (k in who) {

		# get details from the Blast Hit
		id <- hits$PROBE_ID[k]
		nodeID <- sub( "_length.+","",id)
		coverage <- as.numeric( sub( ".+cov_", "", id))
		pStart <- hits$P_FIRST[k]
		pStop <- hits$P_LAST[k]
		sStart <- hits$S_BEG[k]
		sStop <- hits$S_END[k]
		strand <- hits$STRAND[k]

		# get details from the contig
		ptr <- match( id, fa$desc, nomatch=0)
		if ( ptr == 0) stop( paste( "Can't find Contig: ", id))
		seq <- fa$seq[ptr]
		len <- nchar( seq)
		if ( strand == "-") {
			contigEnd <- sStop + (pStart - 1)
			contigBeg <- sStart - (len - pStop)
			seq <- myReverseComplement(seq)
			nodeID <- paste( "RevC(", nodeID,")",sep="")
		} else {
			contigBeg <- sStart - pStart + 1
			contigEnd <- sStop + len - pStop
		}

		# trim to the current window
		contigBeg <- as.integer( max( contigBeg, loBase))
		contigEnd <- as.integer( min( contigEnd, hiBase))

		# what is the tallest/ lowest contig yet shown in this region, and update for this contig
		depthPtrs <- match(contigBeg,names(topdepth)) : match(contigEnd,names(topdepth)) 
		deepestTop <- max( topdepth[ depthPtrs])
		lowestBot <- min( botdepth[ depthPtrs])
		myDepth <- max( coverage, (hiY-loY)*0.1) * 0.9
		if ( deepestTop > 0) {
			# there is something else here...
			if ( lowestBot > (myDepth * 2)) {
				deepest <- 0
			} else {
				deepest <- deepestTop + (hiY-loY)*0.01
			}
		} else {
			deepest <- 0
		}
		contigYlo <- deepest
		contigYhi <- contigYlo + myDepth
		topdepth[ depthPtrs] <- contigYhi
		botdepth[ depthPtrs] <- contigYlo

		# draw the entire contig as an outline
		rect( contigBeg, contigYlo, contigEnd, contigYhi, border=1, lwd=4)
		# draw the Blast hit part as filled
		rect( sStart, contigYlo, sStop, contigYhi, border=1, density=15, lwd=4)
		textX <- (sStart+sStop)/2
		if (textX < loBase) textX <- loBase + (hiBase-loBase)*0.05
		if (textX > hiBase) textX <- hiBase - (hiBase-loBase)*0.05
		textY <- (contigYlo+contigYhi)/2
		if ( !is.null(pos) && pos == 1) textY <- contigYlo
		if ( !is.null(pos) && pos == 3) textY <- contigYhi
		text( textX, textY, nodeID, cex=cex, col=col, pos=pos, font=2)
		#boxed.labels( textX, textY, nodeID, bg=1, cex=cex, col=3, font=2, ypad=1.9)
	}
}
 



`showCGpct` <- function( seqID=getCurrentSeqMap()$SEQ_ID[1], window=100, col=3, lwd=1, bases=T, axis=TRUE) {

	plotLimits <- par("usr")
	loBase <- as.integer( plotLimits[1])
	hiBase <- as.integer( plotLimits[2])
	loY <- plotLimits[3]
	hiY <- plotLimits[4]
	NX <- hiBase - loBase + 1

	chrDNA <- getFastaSeqFromFilePath( genomicFastaFile, seqID=seqID)
	myBases <- strsplit( substr( chrDNA, loBase, hiBase), split="")[[1]]
	names( myBases) <- loBase:hiBase
	half <- round( window/2)
	halfm1 <- half - 1
	step <- ceiling( window/10)
	myXpts <- seq( half, length(myBases)-half, by=step)
	myCGpct <- sapply( myXpts, function(x) {
				tbl <- table( myBases[ (x-halfm1) : (x+halfm1)])
				myAT <- sum( tbl[ match( c("A","T"), names(tbl), nomatch=0)])
				myCG <- sum( tbl[ match( c("C","G"), names(tbl), nomatch=0)])
				return( myCG * 100 / (myCG + myAT))
			})
	# scale this 0..100 to fit the plot
	showCG <- myCGpct * (hiY/100)
	lines( as.numeric( names( myBases[myXpts])), showCG, col=col, lwd=lwd)

	if (bases) {
		baseYs <- loY * c( 0.15,0.10, 0.10,0.15)
		colorYs <- c( 2,4,"goldenrod",3)
		basePtrs <- tapply( 1:length(myBases), factor(myBases), FUN=NULL)
		points( as.numeric(names(myBases)), baseYs[basePtrs], pch=".", col=colorYs[basePtrs], cex=4)
	}

	if (axis) {
		yShow <- seq( 0, 100, by=20)
		at <- seq( loY, hiY, length.out=length(yShow))
		axis( side=4, at=at, label=yShow, col=col, col.tick=col, col.axis=col)
		mtext( "CG Percentage", side=4, line=-1, col=col, font=2)
	}
}


`interestingContigs` <- function( sampleID, results.path=file.path( resultsPath, "VelvetContigs"), min.bases=5) {

	velvetPath <- file.path( results.path, sampleID)
	contigsFastaFile <- file.path( velvetPath, "contigs.fa")
	bestHitFile <- file.path( velvetPath, "Pf.BestBlastHit.csv")
	hits <- read.csv( bestHitFile, as.is=T)

	# load the contigs
	fa <- loadFasta( contigsFastaFile)
	contigLens <- nchar( fa$seq)

	# visit each, and see if there is some fragment of the contig that Blast couldn't align well
	gmap <- subset( getCurrentGeneMap(), REAL_G == TRUE)
	out <- data.frame()
	extraBases <- vector()

	cat( "\n")
	for (k in 1:nrow(hits)) {

		# get details from the Blast Hit
		id <- hits$PROBE_ID[k]
		nodeID <- sub( "_length.+","",id)
		coverage <- as.numeric( sub( ".+cov_", "", id))
		pStart <- hits$P_FIRST[k]
		pStop <- hits$P_LAST[k]
		sStart <- hits$S_BEG[k]
		sStop <- hits$S_END[k]
		strand <- hits$STRAND[k]

		# get details from the contig
		ptr <- match( id, fa$desc, nomatch=0)
		if ( ptr == 0) stop( paste( "Can't find Contig: ", id))
		len <- contigLens[ptr]
		if ( strand == "-") {
			contigBeg <- sStart - len + pStop   #sStart - pStart + 1
			contigEnd <- sStop + pStart - 1     #sStop + len - pStop
			nodeID <- paste( "RevC(", nodeID,")",sep="")
		} else {
			contigBeg <- sStart - pStart + 1
			contigEnd <- sStop + len - pStop
		}

		# if the contig is bigger than the Blast aligned part, that is interesting
		nExtra <- 0
		if ( (extra <- (sStart - min.bases - contigBeg)) > 0) nExtra <- extra
		if ( (extra <- (contigEnd - sStop - min.bases)) > 0) nExtra <- nExtra + extra
		if ( nExtra > 0) {
			out <- rbind( out, hits[ k, ])
			extraBases <- c( extraBases, nExtra)
			cat( "\r", k, nExtra)
		}
	}

	# give the proximal GeneID for each interesting contig
	out$GENE_ID <- NA
	out$EXTRA_BASES <- extraBases
	midpts <- (out$S_BEG + out$S_END) / 2
	tapply( 1:nrow(out), factor( out$SEQ_ID), function(x) {
		gmp <- subset( gmap, SEQ_ID == out$SEQ_ID[ x[ 1]])
		ptrs <- findInterval( midpts[x], gmp$POSITION, all.inside=T)
		out$GENE_ID[x] <<- gmp$GENE_ID[ptrs]
	})

	# weighted rank by score and unaligned bases...
	rank1 <- rank2 <- 1:nrow(out)
	ord <- order( out$EXTRA_BASES, decreasing=T)
	rank2[ ord] <- 1:nrow(out)
	ord <- order( out$SCORE, decreasing=T)
	rank1[ ord] <- 1:nrow(out)
	rank <- (rank1 + rank2)/2
	ord <- order( rank)
	out <- out[ ord, ]
	rownames(out) <- 1:nrow(out)

	return( out)
}
