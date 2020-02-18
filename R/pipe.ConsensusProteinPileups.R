# pipe.ConsensusProteinPileups -- use raw reads as peptides to determine true protein sequence in a field isolate

`pipe.ConsensusProteinPileups` <- function( sampleID, geneID, geneName=geneID, optionsFile="Options.txt",
						results.path=NULL, exon=NULL, maxNoHits.pileup=1000000,
						max.depth=60, txt.cex=0.25, forceSetup=FALSE, maxNoHits.setup=1000000,
						mode=c("normal","realigned"), plotOnly=FALSE, max.drawnPerSite=3,
						trim5.aligns=0, trim3.aligns=0, trim5.nohits=0, trim3.nohits=0,
						draw.box=FALSE, chunkSize.pileup=50000, useCutadapt=FALSE) {

	require(Biostrings)

	# make sure we were given valid parameters..
	if ( length(sampleID) > 1) {
		sampleID <- sampleID[1]
		cat( "\nWarning:  only one sampleID allowed")
	}
	gmap <- getCurrentGeneMap()
	if ( ! (geneID %in% gmap$GENE_ID)) {
		cat( "\nNot a known geneID: ", geneID)
		return(NULL)
	}

	if ( is.null(results.path)) results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	peptide.path <- file.path( results.path, "ConsensusProteins", sampleID)
	if ( ! file.exists( peptide.path)) dir.create( peptide.path, recursive=T)

	# step 1: make sure all needed peptide and pileup files are in place
	nAA <- pipe.ConsensusProteinSetup( sampleID, geneID=geneID, geneName=geneName, results.path=results.path, 
				exon=exon, maxNoHits=maxNoHits.setup, forceSetup=forceSetup,
				trim5.aligns=trim5.aligns, trim3.aligns=trim3.aligns, 
				trim5.nohits=trim5.nohits, trim3.nohits=trim3.nohits,
				useCutadapt=useCutadapt)
	if ( nAA < 1) {
		cat( "\nSetting up for CPP gave an empty protein sequence..")
		return(NULL)
	}

	consensusAAfile <- file.path( peptide.path, paste( sampleID, geneName, "ConsensusAA.fasta", sep="."))
	
	mode <- match.arg( mode)
	if ( mode == "realigned") {
		consensusAAfile <- file.path( peptide.path, paste( sampleID, geneName, "RealignedConsensusAA.fasta", sep="."))
	}

	# step 2:  run the tool that aligns peptides to a construct

	# better luck drawing if we close the current graphics window first
	if ( dev.cur() > 1) dev.off()
	x11Width <- max( round( nAA/550 * 25), 10)
	x11Height <- max( round( max.depth/10, digits=1), 6)
	pdfWidth <- max( round( nAA/550 * 20), 10)
	pdfHeight <- max( round( max.depth/12, digits=1), 5)
	checkX11( bg='white', width=x11Width, height=x11Height)

	consensusAns <<- proteinConstructPeptidePileups( sampleID, geneName=geneName, constructFile=consensusAAfile, 
				peptide.path=peptide.path, txt.cex=txt.cex, maxNoHits=maxNoHits.pileup, 
				max.depth=max.depth, max.drawnPerSite=max.drawnPerSite, mode=mode, draw.box=draw.box,
				chunkSize=chunkSize.pileup)

	# record metrics to the audit file
	writeAuditRecord( peptide.path, sampleID, geneName, mode="Pileup", info=consensusAns)

	pdfFile1 <- file.path( peptide.path, paste( sampleID, geneName, "PeptidePileups.pdf", sep="."))
	pdfFile2 <- file.path( peptide.path, paste( sampleID, geneName, "FinalConsensus.pdf", sep="."))
	if ( mode == "realigned") {
		pdfFile1 <- file.path( peptide.path, paste( sampleID, geneName, "RealignedPeptidePileups.pdf", sep="."))
		pdfFile2 <- file.path( peptide.path, paste( sampleID, geneName, "RealignedFinalConsensus.pdf", sep="."))
	}
	dev.print( pdf, pdfFile1, width=pdfWidth, height=pdfHeight)
	if ( plotOnly) return()

	consensusSaveFile <- file.path( peptide.path, paste( sampleID, geneName, "ConsensusAnswer.rda", sep="."))
	save( consensusAns, file=consensusSaveFile)

	# step 3:  summarize how the consensus differs from what it used as its construct
	differences <- proteinConstructPileupSummary( consensusSaveFile, sampleID=sampleID, geneName=geneName, 
						txt.cex=txt.cex)
	dev.print( pdf, pdfFile2, width=pdfWidth, height=pdfHeight)
	return( invisible( differences))
}


`pipe.ConsensusProteinSetup` <- function( sampleID, geneID, geneName=geneID, maxNoHits=1000000, 
					optionsFile="Options.txt", results.path=NULL, exon=NULL, 
					trim5.aligns=0, trim3.aligns=0, trim5.nohits=0, trim3.nohits=0,
					forceSetup=FALSE, useCutadapt=FALSE, ...) {

	if ( is.null(results.path)) results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	peptide.path <- file.path( results.path, "ConsensusProteins", sampleID)
	if ( ! file.exists( peptide.path)) dir.create( peptide.path, recursive=T)

	madeAnyFiles <- FALSE

	# Step 1:  make sure that all the needed raw read and peptide files are ready
	#		allow use of trimmed reads for the noHits by default
	nohitReadsFile <- file.path( results.path, "fastq", paste( sampleID, "noHits.trimmed.fastq.gz", sep="."))
	nohitReadsFile2 <- file.path( results.path, "fastq", paste( sampleID, "noHits.fastq.gz", sep="."))
	nohitPeptidesFile <- file.path( peptide.path, paste( sampleID, "NoHits", "RawReadPeptides.txt", sep="."))
	geneReadsFile <- file.path( results.path, "fastq", paste( sampleID, geneName, "fastq.gz", sep="."))
	geneReadsTrimmedFile <- file.path( results.path, "fastq", paste( sampleID, geneName, "trimmed.fastq.gz", sep="."))
	genePeptidesFile <- file.path( peptide.path, paste( sampleID, geneName, "RawReadPeptides.txt", sep="."))
	auditFile <- file.path( peptide.path, paste( sampleID, geneName, "AuditTrail.txt", sep="."))

	if ( useCutadapt) {
		if ( ! file.exists( nohitReadsFile)) {
			cat( "\nCutting Adapters off 'NoHit' reads..")
			cutadapt( file1=basename(nohitReadsFile2), path=dirname(nohitReadsFile))
		}
	}
	if ( ! file.exists( nohitReadsFile)) nohitReadsFile <- nohitReadsFile2
	if ( ! file.exists( nohitReadsFile)) {
		stop( "Alignment pipeline 'NoHits' fastq file not found.  Run main pipeline first.")
	}

	if ( forceSetup || ! file.exists( geneReadsFile)) {
		cat( "\nExtracting Gene alignments..")
		nAligns <- pipe.GatherGeneAlignments( sampleID, geneID, tail=500, asFASTQ=T, fastq.keyword=geneName)
		if (nAligns < 1) {
			cat( "\nZero reads aligned to gene..")
			return( 0)
		}
		# new file, so make sure we remake its decendants
		file.delete( genePeptidesFile)
		madeAnyFiles <- TRUE
	}
	if ( forceSetup || ! file.exists( genePeptidesFile)) {
		cat( "\nConverting Gene alignments to peptides..")
		if (useCutadapt) {
			cat( "\nCutting Adapters off 'Gene' reads..")
			cutadapt( file1=basename(geneReadsFile), path=dirname(geneReadsFile))
			geneReadsFile <- geneReadsTrimmedFile
		}
		nPeptides <- fastqToPeptides( geneReadsFile, genePeptidesFile, chunk=100000, lowComplexityFilter=FALSE,
				trim5=trim5.aligns, trim3=trim3.aligns)
		if (nPeptides < 1) {
			cat( "\nGene alignments gave zero peptides..")
			return( 0)
		}
		madeAnyFiles <- TRUE
	}
	if ( forceSetup || ! file.exists( nohitPeptidesFile)) {
		cat( "\nConverting 'NoHit' reads to peptides..  Takes quite a while..")
		cat( "\nBase Trimming:     5' =", trim5.nohits, "    3' =", trim3.nohits)
		multicore.setup(10)
		fastqToPeptides( nohitReadsFile, nohitPeptidesFile, chunk=100000, maxPeptides=maxNoHits,
				lowComplexityFilter=TRUE, trim5=trim5.nohits, trim3=trim3.nohits, ...)
		madeAnyFiles <- TRUE
	}
	if ( forceSetup || ! file.exists( auditFile)) {
		createAuditFile( peptide.path, sampleID, geneName)
	}

	# Step 2:  extract the consensus sequence from the Bowtie alignment results
	exonSuffix <- "Full"
	gmap <- subset( getCurrentGeneMap(), GENE_ID == geneID)
	# using CDS instead of Exons
	emap <- subset( getCurrentCdsMap(), GENE_ID == geneID)
	exonStart <- min( emap$POSITION)
	exonStop <- max( emap$END)
	if (! is.null(exon)) {
		Nexons <- nrow( emap)
		thisExon <- as.integer(exon)
		if( ! (thisExon %in% 1:Nexons)) {
			cat( "\nError:  'exon' must be an integer in 1..",Nexons)
			stop()
		}
		exonSuffix <- paste( "Exon", thisExon, sep="")
		if ( gmap$STRAND == "-") {
			ord <- order( emap$POSITION, decreasing=T)
			myExon <- ord[ thisExon]
			exonStart <- emap$POSITION[myExon]
			exonStop <- emap$END[myExon]
		} else {
			ord <- order( emap$POSITION, decreasing=F)
			myExon <- ord[ thisExon]
			exonStart <- emap$POSITION[myExon]
			exonStop <- emap$END[myExon]
		}
	}
	consensusAAfile <- file.path( peptide.path, paste( sampleID, geneName, "ConsensusAA.fasta", sep="."))
	consensusDNAfile <- file.path( peptide.path, paste( sampleID, geneName, "ConsensusDNA.fasta", sep="."))
	consensusBASEfile <- file.path( peptide.path, paste( sampleID, geneName, "ConsensusBaseCalls.txt", sep="."))
	if ( forceSetup || ! all( file.exists( c( consensusBASEfile, consensusAAfile, consensusDNAfile)))) {
		cat( "\nTurning alignment pileups into consensus base calls..")
		ans <- pipe.ConsensusBaseCalls( sampleID, geneID, start=exonStart, stop=exonStop, as.cDNA=TRUE, 
				noReadCalls="blank")
		callsTable <- ans$callsTable
		myAAvec <- callsTable$AA
		myAAvec[ is.na(myAAvec)] <- ""
		myDNAvec <- callsTable$DNA
		myDNAvec[ is.na(myDNAvec)] <- ""
		myAA <- paste( myAAvec, collapse="")
		myDNA <- paste( myDNAvec, collapse="")
		myDesc <- paste( sampleID, geneName, exonSuffix, sep="_")
		writeFasta( as.Fasta( myDesc, myAA), consensusAAfile, line.width=100)
		writeFasta( as.Fasta( myDesc, myDNA), consensusDNAfile, line.width=100)
		write.table( callsTable, consensusBASEfile, sep="\t", quote=F, row.names=FALSE)
		madeAnyFiles <- TRUE
		nAA <- nchar( myAA)
	} else {
		fa <- loadFasta( consensusAAfile, verbose=F)
		nAA <- nchar( fa$seq)
	}
	if (madeAnyFiles) cat( "\nConsensus Prep done for sample: ", sampleID, "\tGene: ", geneID, "\n")
	return( nAA)
}


`inspectConsensus` <- function( sampleID, geneName="Var2csa", context=NULL, extra.rows=3,
				readingFrame=c("BestFrame","Frame1","Frame2","Frame3"), 
				optionsFile="Options.txt", results.path=NULL) {

	# get the current consensus base calls file
	if ( is.null(results.path)) results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	peptide.path <- file.path( results.path, "ConsensusProteins", sampleID)
	consensusBASEfile <- file.path( peptide.path, paste( sampleID, geneName, "ConsensusBaseCalls.txt", sep="."))
	calls <- read.delim( consensusBASEfile, as.is=T)
	aaAns <- consensusTranslation( calls$DNA)
	readingFrame <- match.arg( readingFrame)
	aaSeq <- paste( aaAns[[readingFrame]], collapse="")
	dnaSeq <- paste( calls$DNA, collapse="")
	calls$AA <- aaAns[[readingFrame]]
	rownames(calls) <- 1:nrow(calls)

	# there is a tiny chance of no indels at all, which makes that column of "" look like NA
	# when we read in the call table
	if ( all( is.na( calls$IndelDetails))) calls$IndelDetails <- ""

	# find that AA context
	context <- toupper( context)
	contextStart <- gregexpr( context, aaSeq, fixed=T)[[1]]
	if ( contextStart[1] < 1) {
		cat( "\nFailed to find 'context' in AA sequence: ", context)
		return(NULL)
	}
	if ( length(contextStart) > 1) {
		cat( "\nNon-unique 'context' in AA sequence: ", context, "\tN_times: ", length(contextStart))
		cat( "\nSites at: ", contextStart)
		return(NULL)
	}
	contextStop <- contextStart + nchar( context) - 1
	contextStr <- substr( aaSeq, contextStart, contextStop)
	cat( "\nContext location: ", contextStart, contextStop, contextStr)

	dnaContextStart <- (contextStart-1) * 3 + 1 - extra.rows
	dnaContextStop <- dnaContextStart + nchar(context) * 3 + extra.rows*2
	cat( "\nBase Call Profile:\n")
	print( calls[ rownames(calls) %in% (dnaContextStart:dnaContextStop), ])
}


`modifyConsensus` <- function( sampleID, geneName="Var2csa", 
				command=c("delete", "insert", "replace", "frameshift","append"), 
				location=NULL, seq, readingFrame=c("BestFrame","Frame1","Frame2","Frame3"), 
				optionsFile="Options.txt", results.path=NULL) {

	# get the current consensus base calls file
	if ( is.null(results.path)) results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	peptide.path <- file.path( results.path, "ConsensusProteins", sampleID)
	consensusBASEfile <- file.path( peptide.path, paste( sampleID, geneName, "ConsensusBaseCalls.txt", sep="."))
	calls <- read.delim( consensusBASEfile, as.is=T)
	aaAns <- consensusTranslation( calls$DNA)
	readingFrame <- match.arg( readingFrame)
	aaSeq <- paste( aaAns[[readingFrame]], collapse="")
	dnaSeq <- paste( calls$DNA, collapse="")
	calls$AA <- aaAns[[readingFrame]]
	rownames(calls) <- 1:nrow(calls)	

	# there is a tiny chance of no indels at all, which makes that column of "" look like NA
	# when we read in the call table
	if ( all( is.na( calls$IndelDetails))) calls$IndelDetails <- ""

	# what command do we want
	command <- match.arg( command)
	
	# what rows do we want
	if ( is.null(location)) stop( "Explicit row number location (range) required")
	if( is.character(location)) {
		rowNumbers <- as.integer( strsplit( gsub(" ","",location), split=":", fixed=T)[[1]])
		if ( length( rowNumbers) == 2) rowNumbers <- rowNumbers[1] : rowNumbers[2]
	} else {
		rowNumbers <- location
	}
	rowNumbers <- intersect( rowNumbers, 1:nrow(calls))

	# force complete codon edits!
	if ( (!(command %in% c("frameshift","append"))) && rowNumbers[1] %% 3 != 1) {
		cat( "\nLocation 1 for edits must be exact codon start boundary!!")
		stop()
	}
	if ( command == "append" && rowNumbers[1] %% 3 != 0) {
		cat( "\nLocation for append must be exact codon end boundary!!")
		stop()
	}
	if ( command %in% c( "delete","replace") && rowNumbers[length(rowNumbers)] %% 3 != 0) {
		cat( "\nLocation 2 for edits must be exact codon end boundary!!")
		stop()
	}
	cat( "\nCommand: ", command, "\tRows: ", range(rowNumbers))

	newCalls <- NULL
	oldSeq <- ""
	newSeq <- toupper(seq)
	protLen <- sum( calls$AA != "", na.rm=T)

	if ( command %in% c( "delete", "frameshift")) {
		newCalls <- calls[ -rowNumbers, ]
		cat( "\nDeleted ", length(rowNumbers), " rows.")
		oldSeq <- calls$AA[ rowNumbers]
	} else {
		# we need the sequence
		newAAseq <- toupper( seq)
		newDNAseq <- AAtoCodonOptimizedDNA( newAAseq, dnaSeq)
		newDNAvector <- strsplit( newDNAseq, split="")[[1]]
		nNewBases <- length(newDNAvector)
		newAAvector <- rep.int( "", nNewBases)
		newAAvector[ seq( 2, nNewBases, by=3)] <- strsplit( newAAseq, split="")[[1]]
	
		if ( command == "replace") {
			if ( nNewBases != length( rowNumbers)) {
				cat( "\nReplacement sequence is wrong size for location range!")
				stop()
			}
			# OK to replace
			newCalls <- calls
			newCalls$Ref[ rowNumbers] <- 1
			newCalls[ rowNumbers, c("A","C","G","T","N","Indel")] <- 0
			newCalls$DNA[ rowNumbers] <- newDNAvector
			newCalls$AA[ rowNumbers] <- newAAvector
			newCalls$IndelDetails[ rowNumbers] <- ""
			cat( "\nReplaced ", nNewBases, " rows with new sequence: ", newAAseq)
			oldSeq <- calls$AA[ rowNumbers]
		}
		if ( command == "insert") {	# insert 'before' X
			part1 <- part2 <- data.frame()
			if (rowNumbers[1] > 1) part1 <- calls[ 1:(rowNumbers[1]-1), ]
			if (rowNumbers[1] <= nrow(calls)) part2 <- calls[ rowNumbers[1]:nrow(calls), ]
			copyFrom <- min( rowNumbers[1], nrow(calls))
			newPart <- calls[ rep.int( copyFrom, nNewBases), ]
			newPart$Ref <- 1
			newPart[ , c("A","C","G","T","N","Indel")] <- 0
			newPart$DNA <- newDNAvector
			newPart$AA <- newAAvector
			newCalls$IndelDetails <- ""
			newCalls <- rbind( part1, newPart, part2)
		}
		if ( command == "append") {	# insert 'after' X
			part1 <- part2 <- data.frame()
			part1 <- calls[ 1:rowNumbers[1], ]
			if (rowNumbers[1] < nrow(calls)) part2 <- calls[ (rowNumbers[1]+1):nrow(calls), ]
			copyFrom <- min( rowNumbers[1], nrow(calls))
			newPart <- calls[ rep.int( copyFrom, nNewBases), ]
			newPart$Ref <- 1
			newPart[ , c("A","C","G","T","N","Indel")] <- 0
			newPart$DNA <- newDNAvector
			newPart$AA <- newAAvector
			newCalls$IndelDetails <- ""
			newCalls <- rbind( part1, newPart, part2)
		}
	}
	if ( ! is.null( newCalls)) {
		cat( "\nUpdating sequence after edits...")
		aaAns <- consensusTranslation( newCalls$DNA)
		dnaSeq <- paste( newCalls$DNA, collapse="")

		# if we did an explicit frame shift, force keeping RF 1, don't let the consensus pick the best!
		if ( command == "frameshift") {
			newCalls$AA <- aaAns$Frame1
			aaSeq <- paste( aaAns$Frame1, collapse="")
		} else {
			newCalls$AA <- aaAns[[readingFrame]]
			aaSeq <- paste( aaAns[[readingFrame]], collapse="")
		}
		newCalls[ ,c('Frame1','Frame2','Frame3')] <- aaAns[,c("Frame1","Frame2","Frame3")]
		rownames(newCalls) <- 1:nrow(newCalls)

		# save the current prior to overwriting..
		cat( "\nSaving previous .FASTA and BaseCalls..")
		consensusAAfile <- file.path( peptide.path, paste( sampleID, geneName, "ConsensusAA.fasta", sep="."))
		consensusDNAfile <- file.path( peptide.path, paste( sampleID, geneName, "ConsensusDNA.fasta", sep="."))
		saveAAfile <- file.path( peptide.path, paste( sampleID, geneName, "SaveConsensusAA.fasta", sep="."))
		saveDNAfile <- file.path( peptide.path, paste( sampleID, geneName, "SaveConsensusDNA.fasta", sep="."))
		saveBASEfile <- file.path( peptide.path, paste( sampleID, geneName, "SaveConsensusBaseCalls.txt", sep="."))
		file.copy( consensusAAfile, saveAAfile, overwrite=T)
		file.copy( consensusDNAfile, saveDNAfile, overwrite=T)
		file.copy( consensusBASEfile, saveBASEfile, overwrite=T)
		oldFA <- loadFasta( consensusAAfile, verbose=F)
		myDesc <- oldFA$desc
		writeFasta( as.Fasta( myDesc, aaSeq), consensusAAfile, line.width=100)
		writeFasta( as.Fasta( myDesc, dnaSeq), consensusDNAfile, line.width=100)
		write.table( newCalls, consensusBASEfile, sep="\t", quote=F, row.names=TRUE)
		cat( "\nWrote new edited .FASTA and BaseCalls for sampleID: ", sampleID, "\n")
		protLen <- sum( newCalls$AA != "", na.rm=T)
	} else {
		cat( "\nNo modifications made...\n")
	}

	# record what we did
	writeAuditRecord( peptide.path, sampleID, geneName, mode="Modify", 
			info=list( "Command"=command, "Location"=rowNumbers, "OldSequence"=oldSeq,
			"NewSequence"=newSeq, "ReadingFrame"=readingFrame, 
			"Length_AA"=protLen))
}


`realignConsensus` <- function( sampleID, geneID="PF3D7_1200600", geneName="Var2csa", 
				readingFrame=c("BestFrame","Frame1","Frame2","Frame3"), 
				optionsFile="Options.txt", annotationFile="Annotation.txt", results.path=NULL,
				exon=NULL, extra.fastq.keyword=NULL, useCutadapt=FALSE) {

	require(Biostrings)

	# given a manually improved consensus sequence, redo the Bowtie alignment step to this custome protein target
	if ( is.null(results.path)) results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	peptide.path <- file.path( results.path, "ConsensusProteins", sampleID)
	consensusAAfile <- file.path( peptide.path, paste( sampleID, geneName, "ConsensusAA.fasta", sep="."))
	consensusDNAfile <- file.path( peptide.path, paste( sampleID, geneName, "ConsensusDNA.fasta", sep="."))
	consensusBASEfile <- file.path( peptide.path, paste( sampleID, geneName, "ConsensusBaseCalls.txt", sep="."))

	# step 1:  make sure the consensus data match the FASTA files of the same data
	cat( "\nValidating current consensus constructs..")
	calls <- read.delim( consensusBASEfile, as.is=T)
	aaAns <- consensusTranslation( calls$DNA)
	readingFrame <- match.arg( readingFrame)
	aaSeq <- paste( aaAns[[readingFrame]], collapse="")
	dnaSeq <- paste( calls$DNA, collapse="")
	aaFA <- loadFasta( consensusAAfile, verbose=F)
	dnaFA <- loadFasta( consensusDNAfile, verbose=F)
	myDesc <- dnaFA$desc
	aaMatch <- (aaSeq == aaFA$seq)
	# for the DNA, ignore last bases after last AA
	ncDNAcheck <- nchar( aaSeq) * 3
	dnaMatch <- ( substr( dnaSeq, 1, ncDNAcheck) ==  substr( dnaFA$seq, 1, ncDNAcheck))
	if ( ! aaMatch) cat( "\nAmino acid FASTA sequence does not match consensus base calls AA column!")
	if ( ! dnaMatch) cat( "\nNucleotide FASTA sequence does not match consensus base calls DNA column!")
	if ( ! all( c( aaMatch, dnaMatch))) stop( "Fix sequence disparities before realignment can run")

	# step 2:  make a Bowtie index of this tiny sequence
	# but we want to add a bit of flanking data to make the alignments at the edges be more robust...

	# step 2.1: find the flanks of reference DNA around this construct
	gmap <- subset( getCurrentGeneMap(), GENE_ID == geneID)
	if ( nrow(gmap) < 1) stop( "'geneID' problem:  no gene detected..")
	myStart <- gmap$POSITION[1]
	myStop <- gmap$END[1]
	cat( "\nExtracting genomic flanking regions..")
	flankLen <- 50
	refAns <- pipe.ConsensusBaseCalls( sampleID, geneID=geneID, start=myStart, stop=myStop, 
			aaToo=FALSE, as.cDNA=FALSE, utr.tail.length=flankLen, noReadCalls="genomic")
	myRefV <- refAns$ref
	myRefStr <- paste( myRefV, collapse="")
	if ( gmap$STRAND[1] == "-") myRefStr <- myReverseComplement(myRefStr)
	# we need the 'exon' info to grab/find the correct flanking data
	emap <- subset( getCurrentCdsMap(), GENE_ID == geneID)
	Nexons <- nrow( emap)
	if (! is.null(exon)) {
		exon <- as.integer(exon)
		if( ! (exon %in% 1:Nexons)) stop( paste( "\nError:  'exon' must be an integer in 1..",Nexons))
	}
	# if we are doing a single exon, then our construct will be found directly
	if ( is.integer( exon) || Nexons == 1) {
		cat( "\nGetting single exon flanking regions..")
		paAns <- pairwiseAlignment( dnaSeq, myRefStr, type="global-local")
		if ( score( paAns) < nchar(dnaSeq)/2) cat( "  low 'flank search alignment score' warning: ", score(paAns))
		myRefStart <- start( subject( paAns))
		myRefStop <- width( subject( paAns)) + myRefStart - 1
		leftFlank <- substr( myRefStr, max( 1, myRefStart-flankLen), myRefStart-1)
		rightFlank <- substr( myRefStr, myRefStop+1, min( myRefStop+flankLen, nchar(myRefStr)))
	} else {
	# but if doing a multi-exon construct, we have to go find the correct edges separately
		cat( "\nGetting multi exon flanking regions..")
		if (gmap$STRAND[1] == "-") emap <- emap[ rev( 1:Nexons), ]
		firstExonStart <- 1
		firstExonStop <- emap$END[1] - emap$POSITION[1] + 1
		firstExonSeq <- substr( dnaSeq, firstExonStart, firstExonStop)
		paAns <- pairwiseAlignment( firstExonSeq, myRefStr, type="global-local")
		if ( score( paAns) < nchar(firstExonSeq)/2) cat( "  low 'flank search alignment score' warning: ", score(paAns))
		myRefStart <- start( subject( paAns))
		myRefStop <- width( subject( paAns)) + myRefStart - 1
		leftFlank <- substr( myRefStr, max( 1, myRefStart-flankLen), myRefStart-1)
		lastExonStop <- nchar( dnaSeq)
		lastExonStart <- lastExonStop - (emap$END[Nexons] - emap$POSITION[Nexons])
		lastExonSeq <- substr( dnaSeq, lastExonStart, lastExonStop)
		paAns <- pairwiseAlignment( lastExonSeq, myRefStr, type="global-local")
		if ( score( paAns) < nchar(lastExonSeq)/2) cat( "  low 'flank search alignment score' warning: ", score(paAns))
		myRefStart <- start( subject( paAns))
		myRefStop <- width( subject( paAns)) + myRefStart - 1
		rightFlank <- substr( myRefStr, myRefStop+1, min( myRefStop+flankLen, nchar(myRefStr)))
	}
	leftFlankLen <- nchar( leftFlank)
	rightFlankLen <- nchar( rightFlank)

	# step 2.2:  make a temp sequence that has these flanks
	tmpdna <- paste( leftFlank, dnaSeq, rightFlank, sep="", collapse="")
	#tmpdnaV <- strsplit( tmpdna, split="")[[1]]
	consensusDNAflankFile <- sub( "ConsensusDNA", "ConsensusDNAwithFlanks", consensusDNAfile)
	writeFasta( as.Fasta( myDesc, tmpdna), consensusDNAflankFile, line.width=100)
	indexFile <- file.path( peptide.path, "ConsensusProteinIndex")
	cat( "\nMaking Bowtie index from consensus sequence with flanks..")
	callBowtie2Build( buildBowtie2BuildCommandLine( inputFastaFile=consensusDNAflankFile, outputIndexFile=indexFile, 
			optionsFile=optionsFile, verbose=F))

	# step 3:  realign the raw reads against this new target
	cat( "\nCalling Bowtie against consensus sequence..")
	nohitReadsFile <- file.path( results.path, "fastq", paste( sampleID, "noHits.fastq.gz", sep="."))
	geneReadsFile <- file.path( results.path, "fastq", paste( sampleID, geneName, "fastq.gz", sep="."))
	if (useCutadapt) {
		nohitReadsFile <- file.path( results.path, "fastq", paste( sampleID, "noHits.trimmed.fastq.gz", sep="."))
		geneReadsTrimFile <- file.path( results.path, "fastq", paste( sampleID, geneName, "trimmed.fastq.gz", sep="."))
	}
	allReadsFiles <- c(geneReadsFile, nohitReadsFile)
	if ( ! is.null( extra.fastq.keyword)) {
		cat( "  extra FASTQ:  ")
		for ( keyword in extra.fastq.keyword) {
			extraFile <- file.path( results.path, "fastq", paste( sampleID, keyword, "fastq.gz", sep="."))
			allReadsFiles <- c( allReadsFiles, extraFile)
			cat( "  ", keyword)
		}
	}
	bamFile <- file.path( peptide.path, "ConsensusProtein.bam")
	nhFile <- file.path( peptide.path, "ConsensusProtein.noHits.fastq.gz")
	bowtieAns <- fastqToBAM( inputFastqFile=allReadsFiles, outputFile=bamFile, sampleID=sampleID, 
			optionsFile=optionsFile, annotationFile=annotationFile, 
			alignIndex=indexFile, index.path=".", noHitsFile=nhFile, verbose=F)
	nReadsAligned <- bowtieAns$UniqueReads + bowtieAns$MultiReads
	
	# step 4:  get the new consensus from this BAM file
	cat( "\nExtract new consensus from re-aligning reads to previous consensus..")
	ans <- consensusBaseCalls( bamfile=bamFile, genomicFastaFile=consensusDNAflankFile, seqID=myDesc, 
				geneID=NULL, start=leftFlankLen+1, stop=leftFlankLen+nchar(dnaSeq), aaToo=TRUE, 
				noReadCalls="genomic")
	#ans <- consensusBaseCallsToCDNA( ans, geneID=NULL, start=leftFlankLen+1, stop=leftFlankLen+nchar(dnaSeq))
	ans <- consensusBaseCallsToCDNA( ans, geneID=NULL)

	# step 5:  write these out
	calls <- ans$callsTable
	dnaSeq <- paste( calls$DNA, collapse="")
	aaSeq <- paste( calls$AA, collapse="")
	consensusAAfile <- file.path( peptide.path, paste( sampleID, geneName, "RealignedConsensusAA.fasta", sep="."))
	consensusDNAfile <- file.path( peptide.path, paste( sampleID, geneName, "RealignedConsensusDNA.fasta", sep="."))
	consensusBASEfile <- file.path( peptide.path, paste( sampleID, geneName, "RealignedConsensusBaseCalls.txt", sep="."))
	writeFasta( as.Fasta( myDesc, aaSeq), consensusAAfile, line.width=100)
	writeFasta( as.Fasta( myDesc, dnaSeq), consensusDNAfile, line.width=100)
	write.table( calls, consensusBASEfile, sep="\t", quote=F, row.names=TRUE)
	cat( "\nWrote new realigned .FASTA and BaseCalls for sampleID: ", sampleID, "\n")

	# step 6:  turn the alignments back to fastq and then to peptides
	bamFile <- file.path( peptide.path, "ConsensusProtein.sorted.bam")
	fastqFile <- file.path( peptide.path, "ConsensusProtein.AlignedReads.fq.gz")
	tempPeptidesFile <- file.path( peptide.path, "ConsensusProtein.Peptides.txt")
	genePeptidesFile <- file.path( peptide.path, paste( sampleID, geneName, "RawReadPeptides.txt", sep="."))
	bam2fastq( bamfile=bamFile, outfile=fastqFile, verbose=F)
	# since the raw DNA reads wt got aligned are already cutAdapt'ed, no need to do it again...
	fastqToPeptides( filein=fastqFile, fileout=tempPeptidesFile)
	# merge these new peptides into the existing set, instead of overwriting...
	mergePeptideFiles( tempPeptidesFile, genePeptidesFile, outfile=genePeptidesFile, 
			mergeCountsMode="File2")

	# record this step to audit trail
	writeAuditRecord( peptide.path, sampleID, geneName, mode="Realign", 
			info=list( "Length_AA"=nchar(aaSeq), "N_Peptides"=nReadsAligned))

	# lastly, clean up
	system( paste( "rm ", file.path( peptide.path, "ConsensusProtein*bam*")))
	system( paste( "rm ", file.path( peptide.path, "ConsensusProteinIndex*")))
	system( paste( "rm ", file.path( peptide.path, "ConsensusProtein.noHits.fastq.gz")))
	system( paste( "rm ", file.path( peptide.path, "ConsensusProtein.AlignedReads.fq.gz")))
	system( paste( "rm ", consensusDNAflankFile))
	system( paste( "rm ", file.path( peptide.path, paste( sampleID, geneName, "ConsensusDNAwithFlanks.fasta.fai", sep="."))))
	system( paste( "rm ", tempPeptidesFile))
}


`acceptRealignedConsensus` <- function( sampleID, geneID="PF3D7_1200600", geneName="Var2csa", optionsFile="Options.txt", 
				results.path=NULL) {

	if ( is.null(results.path)) results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	peptide.path <- file.path( results.path, "ConsensusProteins", sampleID)

	consensusAAfile <- file.path( peptide.path, paste( sampleID, geneName, "ConsensusAA.fasta", sep="."))
	consensusDNAfile <- file.path( peptide.path, paste( sampleID, geneName, "ConsensusDNA.fasta", sep="."))
	consensusBASEfile <- file.path( peptide.path, paste( sampleID, geneName, "ConsensusBaseCalls.txt", sep="."))

	realignedAAfile <- file.path( peptide.path, paste( sampleID, geneName, "RealignedConsensusAA.fasta", sep="."))
	realignedDNAfile <- file.path( peptide.path, paste( sampleID, geneName, "RealignedConsensusDNA.fasta", sep="."))
	realignedBASEfile <- file.path( peptide.path, paste( sampleID, geneName, "RealignedConsensusBaseCalls.txt", sep="."))
	if ( ! all( file.exists( c( realignedAAfile, realignedDNAfile, realignedBASEfile)))) {
		cat( "\nRe-aligned files not found!  Tried: ", basename( realignedAAfile))
		stop()
	}

	saveAAfile <- file.path( peptide.path, paste( sampleID, geneName, "SaveConsensusAA.fasta", sep="."))
	saveDNAfile <- file.path( peptide.path, paste( sampleID, geneName, "SaveConsensusDNA.fasta", sep="."))
	saveBASEfile <- file.path( peptide.path, paste( sampleID, geneName, "SaveConsensusBaseCalls.txt", sep="."))
	file.copy( consensusAAfile, saveAAfile, overwrite=T)
	file.copy( consensusDNAfile, saveDNAfile, overwrite=T)
	file.copy( consensusBASEfile, saveBASEfile, overwrite=T)

	file.copy( realignedAAfile, consensusAAfile, overwrite=T)
	file.copy( realignedDNAfile, consensusDNAfile, overwrite=T)
	file.copy( realignedBASEfile, consensusBASEfile, overwrite=T)
	file.delete( c( realignedAAfile, realignedDNAfile, realignedBASEfile))

	consensusPILEUPfile <- file.path( peptide.path, paste( sampleID, geneName, "PeptidePileups.pdf", sep="."))
	consensusFINALfile <- file.path( peptide.path, paste( sampleID, geneName, "FinalConsensus.pdf", sep="."))
	realignedPILEUPfile <- file.path( peptide.path, paste( sampleID, geneName, "RealignedPeptidePileups.pdf", sep="."))
	realignedFINALfile <- file.path( peptide.path, paste( sampleID, geneName, "RealignedFinalConsensus.pdf", sep="."))
	if ( file.exists( realignedPILEUPfile)) {
		file.copy( realignedPILEUPfile, consensusPILEUPfile, overwrite=T)
		file.delete( realignedPILEUPfile)
	}
	if ( file.exists( realignedFINALfile)) {
		file.copy( realignedFINALfile, consensusFINALfile, overwrite=T)
		file.delete( realignedFINALfile)
	}
}


mergePeptideFiles <- function( infile1, infile2, outfile, mergeCountsMode=c("File1","File2","Add")) {

	df1 <- read.delim( infile1, as.is=T)
	df2 <- read.delim( infile2, as.is=T)

	outPep <- union( df1$Peptide, df2$Peptide)
	outCounts1 <- outCounts2 <- rep.int( 0, length(outPep))

	wh1 <- match( outPep, df1$Peptide, nomatch=0)
	wh2 <- match( outPep, df2$Peptide, nomatch=0)
	outCounts1[ wh1 > 0] <- df1$Count[wh1]
	outCounts2[ wh2 > 0] <- df2$Count[wh2]

	mergeCountsMode <- match.arg( mergeCountsMode)
	if ( mergeCountsMode == "File1") {
		missing <- which( outCounts1 == 0)
		outCounts1[ missing] <- outCounts2[ missing]
		outCounts <- outCounts1
	} else if ( mergeCountsMode == "File2") {
		missing <- which( outCounts2 == 0)
		outCounts2[ missing] <- outCounts1[ missing]
		outCounts <- outCounts2
	} else {
		outCounts <- outCounts1 + outCounts2
	}

	out <- data.frame( "Peptide"=outPep, "Count"=outCounts, stringsAsFactors=FALSE)
	write.table( out, outfile, sep="\t", quote=F, row.names=F)
}


`searchForMisalignedReads` <- function( sampleID, geneName="Var2csa", context=NULL, flank.size=50,
				optionsFile="Options.txt", results.path=NULL) {

	# get the current consensus base calls file
	if ( is.null(results.path)) results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	peptide.path <- file.path( results.path, "ConsensusProteins", sampleID)
	consensusBASEfile <- file.path( peptide.path, paste( sampleID, geneName, "ConsensusBaseCalls.txt", sep="."))
	calls <- read.delim( consensusBASEfile, as.is=T)
	aaAns <- consensusTranslation( calls$DNA)
	readingFrame <- "Frame1"
	aaSeq <- paste( aaAns[[readingFrame]], collapse="")
	dnaSeq <- paste( calls$DNA, collapse="")
	calls$AA <- aaAns[[readingFrame]]
	rownames(calls) <- 1:nrow(calls)

	# find that AA context
	context <- toupper( context)
	contextStart <- gregexpr( context, aaSeq, fixed=T)[[1]]
	if ( contextStart[1] < 1) {
		cat( "\nFailed to find 'context' in AA sequence: ", context)
		return(NULL)
	}
	if ( length(contextStart) > 1) {
		cat( "\nNon-unique 'context' in AA sequence: ", context, "\tN_times: ", length(contextStart))
		cat( "\nSites at: ", contextStart)
		return(NULL)
	}
	contextStop <- contextStart + nchar( context) - 1
	contextStr <- substr( aaSeq, contextStart, contextStop)
	cat( "\nContext location: ", contextStart, contextStop, contextStr)

	# this gives us the DNA in the area around that AA context, plus flanks
	dnaContextStart <- (contextStart-1) * 3 + 1 - flank.size
	dnaContextStop <- dnaContextStart + nchar(context) * 3 + flank.size * 2
	dnaContext <- substr( dnaSeq, dnaContextStart, dnaContextStop)
	dnaContextRC <- myReverseComplement( dnaContext)
	cat( "\nDNA Context to look for: ", dnaContext)

	# get the DNA for every gene
	geneDNAfile <- file.path( peptide.path, "AllGenes.GenomeDNA.fasta")
	if ( ! file.exists( geneDNAfile)) {
		cat( "\nBuilding FASTA of DNA for all genes..")
		genomicFastaFile <- getOptionValue( optionsFile, "genomicFastaFile", verbose=F)
		gmap <- subset( getCurrentGeneMap(), REAL_G == TRUE)
		allGenes <- gmap$GENE_ID
		fa <- gene2Fasta( allGenes, genomicFastaFile, mode="gdna", verbose=T)
		writeFasta( fa, geneDNAfile, line=100)
	} else {
		fa <- loadFasta( geneDNAfile, verbose=F)
	}

	# find the best hits of this DNA context to any gene in the genome
	cat( "\nSearching DNA of all genes..")
	require( Biostrings)
	submat <- nucleotideSubstitutionMatrix()
	scores1 <- pairwiseAlignment( fa$seq, dnaContext, type="local", substitutionMatrix=submat, scoreOnly=T)
	names(scores1) <- fa$desc
	scores2 <- pairwiseAlignment( fa$seq, dnaContextRC, type="local", substitutionMatrix=submat, scoreOnly=T)
	names(scores2) <- fa$desc

	allScores <- sort( c( scores1, scores2), decreasing=T)
	ans <- allScores[ 1:10]

	out <- data.frame( "GENE_ID"=names(ans), "SCORE"= as.numeric(ans), "PRODUCT"=gene2Product( names(ans)), 
				stringsAsFactors=F)
	return(out)
}


`createAuditFile` <- function( peptide.path, sampleID, geneName) {

	auditFile <- auditFileName( peptide.path, sampleID, geneName)
	con <- file( auditFile, open="wt")
	headerTerms <- c( "Command", "DateTime", "Length_AA", "N_Peptides", "SubCommand", "Location_DNA", 
			"OldSeqFrag", "NewSeqFrag")
	headerText <- paste( headerTerms, collapse="\t")
	writeLines( headerText, con=con)
	close( con)
}


`auditFileName` <- function( peptide.path, sampleID, geneName) {

	f <- paste( sampleID, geneName, "AuditTrail.txt", sep=".")
	return( file.path( peptide.path, f))
}


`writeAuditRecord` <- function( peptide.path, sampleID, geneName, mode, info) {

	auditFile <- auditFileName( peptide.path, sampleID, geneName)
	con <- file( auditFile, open="at")

	dat <- date()
	cmnd <- mode
	cmd2 <- protLen <- nPep <- locStr <- oldSeqFrag <- newSeqFrag <- ""
	if ( mode == "Pileup") {
		# info is a list of details from the protein pileup step
		nPep <- info$N_Peptides
		# the protein is a vector of single characters
		protLen <- length( info$Construct)
	}
	if ( mode == "Modify") {
		# info is a list of details from the protein modification step
		cmd2 <- info$Command
		loc <- as.numeric( info$Location)
		locStr <- paste( min(loc,na.rm=T), max(loc,na.rm=T), sep=":")
		oldSeqFrag <- gsub( " ", "", paste( info$OldSequence, collapse=""))
		newSeqFrag <- gsub( " ", "", paste( info$NewSequence, collapse=""))
		protLen <- info$Length_AA
	}
	if ( mode == "Realign") {
		# info is a list of details from the bowtie re-alignment step
		nPep <- info$N_Peptides
		cmd2 <- "Bowtie2"
		protLen <- info$Length_AA
	}
	auditText <- paste( c( cmnd, dat, protLen, nPep, cmd2, locStr, oldSeqFrag, 
				newSeqFrag), collapse="\t")
	writeLines( auditText, con=con)
	close( con)
}


CPP.AuditSummary <- function( sampleID, geneName="Varcsa", results.path=getOptionValue( "Options.txt", "results.path", verbose=F)) {

	path <- file.path( results.path, "ConsensusProteins", sampleID)
	f <- auditFileName( path, sampleID, geneName=geneName)
	if ( ! file.exists( f)) {
		cat( "\nCPP audit file not found: ", f)
		return( NULL)
	}

	tbl <- read.delim( f, as.is=T)
	cat( "\n\nAudit Summary: ", sampleID, "\n")
	cat( "\nCounts by CPP Command Type:")
	print( ans1 <- table( tbl$Command))
	cat( "\nCounts by CPP 'Modify' Sub-Command Type:")
	print( ans2 <- table( tbl$SubCommand[ tbl$Command == "Modify"]))

	# if the gene is VAR2CSA, do more detailed analysis
	if ( toupper(geneName) == "VAR2CSA") {
		dmap <- getVar2csaDomainMap( strain="3D7")
		# make a "findInterval" construct, using the midpoint between domains
		domStarts <- dmap$REF_START
		prevDomStops <- c( 1, dmap$REF_STOP[1:(nrow(dmap)-1)])
		domStarts <- round( (prevDomStops + domStarts) / 2)
		names( domStarts) <- dmap$DOMAIN_ID
		# let's look at all the modifcation lines
		isMOD <- which( tbl$Command == "Modify")
		# turn the "Location_DNA" info into  a AA location
		dnaLocTerms <- strsplit( tbl$Location_DNA[ isMOD], split=":")
		aaCenter <- sapply( dnaLocTerms, function( x) return( round( mean( as.numeric(x), na.rm=T) / 3)))
		# bacause the CPP construct tends to be short early, and grow to full size, just using AA locs
		# will have a bias.  Convert to percentages to do the find, show we are being fair
		domStarts <- domStarts / max( domStarts)
		# more precisely, the length to use is the one from "before" this edit, as it was the location this edit occured on
		aaLens <- as.numeric( tbl$Length_AA[ isMOD-1])
		# we except a typical Var2csa to be ~2650 long.  And we know that the worst, most likely missing region is DBL6
		# thus, lets assume the length at a minimum to be more fair with our fractional guesses.
		aaLens <- pmax( aaLens, rep.int(2600,length(aaLens)))
		aaCenter <- aaCenter / aaLens
		hits <- findInterval( aaCenter, domStarts, all.inside=T)
		domHits <- names( domStarts)[hits]
		domHitTbl <- table( factor( domHits, levels=dmap$DOMAIN_ID))
		domHitPct <- round( domHitTbl * 100 / sum(domHitTbl), digits=1)
		cat( "\nLocations of Modify Operations By Var2csa Domain: \n")
		ans3 <- data.frame( "DomainID"=dmap$DOMAIN_ID, "N_Modify_Ops"=as.numeric(domHitTbl), 
				"Pct_Modify_Ops"=as.numeric(domHitPct), stringsAsFactors=F)
		rownames(ans3) <- 1:nrow(ans3)
		print( ans3)

		return( list( "CPP Commands"=ans1, "Modify_Subcommands"=ans2, "Var2csa Domains"=ans3))
	} else {
		return( list( "CPP Commands"=ans1, "Modify_Subcommands"=ans2))
	}
}

