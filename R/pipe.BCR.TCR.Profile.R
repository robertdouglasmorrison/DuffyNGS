# pipe.BCR.TCR.Profile.R

# turn RNA-seq data into BCR & TCR sequences, to search for V-D-J and CDR3 motifs


`pipe.BCR.TCR.Profile` <- function( sampleID, bcrtcrFastaFile=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt",
				results.path=NULL, spades.path=dirname(Sys.which("spades.py")), doSpades=FALSE,
				spades.mode=c("isolate","rna","meta"), kmerSizes=NULL, spades.args="", 
				doCutadapt=TRUE, cutadaptProgram=Sys.which("cutadapt"), keyword="BCR.TCR", 
				min.aa.length=90, max.aa.length=250, min.v.score=100, verbose=TRUE) {

	if ( is.null( results.path)) results.path <- getOptionValue( optionsFile, "results.path", 
				notfound=".", verbose=F)

	# make sure the reference file of var genes is readable
	if ( ! is.null( bcrtcrFastaFile)) {
		if ( ! file.exists( bcrtcrFastaFile)) {
			cat( "\nError: 'bcrtcrFastaFile' of example BCR.TCR proteins not found: ", bcrtcrFastaFile)
			return(NULL)
		}
	}

	# path for all results
	spades.output.path <- file.path( results.path, "SpadesContigs", sampleID, "BCR.TCR")
	contigsFastaFile <- file.path( spades.output.path, "contigs.fasta")
	peptidesFastaFile <- file.path( spades.output.path, paste( sampleID, keyword, "Peptides.fasta", sep="."))
	proteinsFastaFile <- sub( "Peptides.fasta$", "Proteins.fasta", peptidesFastaFile)
	proteinsTextFile <- sub( "Proteins.fasta$", "BestProteinHits.txt", proteinsFastaFile)
	
	# remove any files we know we will remake
	if (doSpades) file.delete( c( contigsFastaFile, peptidesFastaFile, proteinsFastaFile, proteinsTextFile))

	# get the universe of BCR & TCR genes from the annotation
	setCurrentSpecies( "Hs_grc")
	gmap <- getCurrentGeneMap()
	tcrMap <- gmap[ grep( "^TR[ABDG][VDJC]", gmap$GENE_ID), ]
	bcrMap <- gmap[ grep( "^IG[HKL][VDJCG]", gmap$GENE_ID), ]
	tcrgenes <- sort( unique( c( bcrMap$GENE_ID, tcrMap$GENE_ID)))
	alignedReadsFile <- file.path( results.path, "fastq", paste( sampleID, "BCR.TCR.Hits.fastq.gz", sep="."))
	nohitReadsFile <- file.path( results.path, "fastq", paste( sampleID, "noHits.fastq.gz", sep="."))

	if ( doSpades || !file.exists(contigsFastaFile)) {

		# verify we see the denovo tool executable location
		if ( spades.path == "") {
			cat( "\nError:  path to SPAdes executable program not found..")
			return(NULL)
		}

		# step 1:   gather all aligned reads that land in/near the TCR genes loci
		if ( ! file.exists( alignedReadsFile)) {
			cat( "\nGathering BCR/TCR Gene aligned reads..\n")
			pipe.GatherGeneAlignments( sampleID, genes=tcrgenes, asFASTQ=TRUE, mode="best.one",
						fastq.keyword="BCR.TCR.Hits")
		}

		# in most cases call cutadapt to trim off primer ends
		# if we are doing the Cutadapt pass, check for those files, and call it if we need.
		fastqFiles <- c( "BCR.TCR.Hits", "noHits")
		if (doCutadapt) {
			inputFastqFiles <- c( alignedReadsFile, nohitReadsFile)
			filesToDo <- checkOrCallCutadapt( inputFastqFiles, asMatePairs=FALSE, forceMatePairs=FALSE,
						cutadaptProgram=cutadaptProgram, kmer.size=45, verbose=verbose)
			fastqFiles <- paste( fastqFiles, "trimmed", sep=".")
		}
	
		# step 2:  create contigs of anything that may be TCR reads
		spades.mode <- match.arg( spades.mode)
		pipe.SpadesContigs( sampleID, fastqSource=fastqFiles, folderName="BCR.TCR", 
				spades.mode=spades.mode, kmerSizes=kmerSizes, spades.args=spades.args, 
				spades.path=spades.path, keyword=keyword, makePep=T)

		# cleanup SPAdes temp files
		cleanupSpadesFiles( path=spades.output.path, verbose=F)
	} else {
		cat( "\n\nUsing previously calculated Spades contigs..\n")
	}

	# because of reading frame and stop codon trimming, and that SPAdes has no length cutoff,
	# we may have peptides that are now shorter than the minimun we want, or too long to be a VDJ construct.
	pepFA <- loadFasta( peptidesFastaFile, short=T, verbose=F)
	len <- nchar( pepFA$seq)
	drops1 <- which( len < min.aa.length)
	drops2 <- which( len > max.aa.length)
	drops <- sort( union( drops1, drops2))
	if ( length( drops)) {
		pepFA <- as.Fasta( pepFA$desc[-drops], pepFA$seq[-drops])
		writeFasta( pepFA, peptidesFastaFile, line=100)
	}
	NPEP <- length( pepFA$desc)

	# step 3:  Throw those contigs as peptides against the given set of full length TCR gene proteins
	if ( ! is.null( bcrtcrFastaFile)) {
		cat( "\n\nSearching Contigs for BCR/TCR constructs..")
		proteinAns <- bestSpadesProteins( sampleID, outpath=spades.output.path, proteinFastaFile=bcrtcrFastaFile,
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
			cat( "\nWrote", NPROT, "BCR/TCR proteins as FASTA to: ", basename(proteinsFastaFile))
		}
	
		# at this point, we have a small chance that none look like BCR/TCR genes
		if ( ! NPROT) {
			cat( "\nNo peptides look like BCR/TCR proteins..")
			write.table( proteinAns, proteinsTextFile, sep="\t", quote=F, row.names=F)
			return( invisible( proteinAns))
		}
	}

	# step 4:  call IgBlast
	# at this point, we know which contigs are at least the correct length range, and perhap which look like example BCR/TCR,
	# so we can narrow down the raw contigs to those most likely to be true BCR/TCR constructs
	igblastOutputFile1 <- file.path( spades.output.path, "IgBlast.BCR.Results.txt")
	igblastOutputFile2 <- file.path( spades.output.path, "IgBlast.TCR.Results.txt")
	resultFile1 <- file.path( spades.output.path, paste( sampleID, "Final.BCR.Results.txt", sep="."))
	resultFile2 <- file.path( spades.output.path, paste( sampleID, "Final.TCR.Results.txt", sep="."))

	usableContigIDs <- pepFA$desc
	if ( exists( proteinAns) && NPROT) usableContigIDs <- protFA$desc
	usableContigIDs <- sub( "_FR[1:6]$", "", usableContigs)
	contigsFA <- loadFasta( contigFastaFile)
	keep <- which( contigFA$desc %in% usableContigIDs)
	if ( ! length(keep)) {
		cat( "\n  Warning: no contigs look like BCR/TCR sequences. Not calling IgBlast.")
		file.delete( c( resultFile1, resultFile2))
		return(NULL)
	}
	igblastFA <- as.Fasta( contigFA$desc[keep], contigFA$seq[keep])
	igblastContigsFile <- file.path( spades.output.path, "IgBlast.Contigs.fasta")
	writeFasta( igblastFA, igblastContigsFile, line=100)

	# look for both TCR and BCR by calling the IgBlast tool
	callIgBlast( igblastContigsFile, outfile=igblastOutputFile1, db="IG")
	ansIG <- readIgBlastOutput( igblaseOutputFile1)
	ansIG <- subset( ansIG, V_SCORE >= min.v.score)
	cat( "\n  N_IG Constructs:  ", nrow(ansIG))
	callIgBlast( igblastContigsFile, outfile=igblastOutputFile2, db="TCR")
	ansTCR <- readIgBlastOutput( igblaseOutputFile2)
	ansTCR <- subset( ansTCR, V_SCORE >= min.v.score)
	cat( "\n  N_TCR Constructs: ", nrow(ansTCR))

	# write out final results
	write.table( ansIG, resultFile1, sep="\t", quote=F, row.names=F)
	write.table( ansTCR, resultFile2, sep="\t", quote=F, row.names=F)

	out <- c( "N_BCR"=nrow(ansIG), "N_TCR"=nrow(ansTCR))
	return(out)
}

