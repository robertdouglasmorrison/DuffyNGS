# pipe.HLA.ContigProfile.R

# turn RNA-seq data into a profile of HLA. constructs & expression
# use the SPADES denovo assembly, on anything that looks like an HLA gene


`pipe.HLA.ContigProfile` <- function( sampleID, hlaFastaFile, annotationFile="Annotation.txt", optionsFile="Options.txt",
				results.path=NULL, spades.path=dirname(Sys.which("spades.py")), doSpades=FALSE,
				spades.mode=c("isolate","rna","meta"), kmerSizes=NULL, spades.args="", 
				doCutadapt=TRUE, cutadaptProgram=Sys.which("cutadapt"), keyword="HLA", 
				min.aa.length=100, min.score.per.aa=2, verbose=TRUE) {

	if ( is.null( results.path)) results.path <- getOptionValue( optionsFile, "results.path", 
				notfound=".", verbose=F)

	# make sure the reference file of var genes is readable
	if ( ! file.exists( hlaFastaFile)) {
		cat( "\nError: 'hlaFastaFile' of reference proteins not found: ", hlaFastaFile)
		return(NULL)
	}

	# path for all results
	spades.output.path <- file.path( results.path, "SpadesContigs", sampleID, "HLA")
	contigsFastaFile <- file.path( spades.output.path, "contigs.fasta")
	peptidesFastaFile <- file.path( spades.output.path, paste( sampleID, keyword, "Peptides.fasta", sep="."))
	proteinsFastaFile <- sub( "Peptides.fasta$", "Proteins.fasta", peptidesFastaFile)
	proteinsTextFile <- sub( "Peptides.fasta$", "BestProteinHits.txt", peptidesFastaFile)
	
	# remove any files we know we will remake
	if (doSpades) file.delete( c( contigsFastaFile, peptidesFastaFile, proteinsFastaFile, proteinsTextFile))

	# get the universe of HLA genes
	if ( getCurrentSpecies() != "Hs_grc") stop( "HLA.ContigProfiling is only defind for human samples")
	gmap <- getCurrentGeneMap()
	hlaMap <- gmap[ grep( "^HLA\\-", gmap$GENE_ID), ]
	hlagenes <- sort( unique( hlaMap$GENE_ID))
	alignedReadsFile <- file.path( results.path, "fastq", paste( sampleID, "HLA.Hits.fastq.gz", sep="."))
	nohitReadsFile <- file.path( results.path, "fastq", paste( sampleID, "noHits.fastq.gz", sep="."))

	if ( doSpades || !file.exists(contigsFastaFile)) {

		# verify we see the denovo tool executable location
		if ( spades.path == "") {
			cat( "\nError:  path to SPAdes executable program not found..")
			return(NULL)
		}
		spadesExecFile <= file.path( spades.path, "spades.py")
		if ( ! file.exists( spadesExecFile)) {
			cat( "\nError:  Failed to find SPAdes executable.  Tried: ", spadesExecFile)
			return(NULL)
		}

		# step 1:   gather all aligned reads that land in/near the HLA loci
		if ( ! file.exists( alignedReadsFile)) {
			cat( "\nGathering HLA Gene aligned reads..\n")
			pipe.GatherGeneAlignments( sampleID, genes=hlagenes, asFASTQ=TRUE, mode="best.one",
						fastq.keyword="HLA.Hits")
		}

		# in most cases call cutadapt to trim off primer ends
		# if we are doing the Cutadapt pass, check for those files, and call it if we need.
		fastqFiles <- c( "HLA.Hits", "noHits")
		if (doCutadapt) {
			inputFastqFiles <- c( alignedReadsFile, nohitReadsFile)
			filesToDo <- checkOrCallCutadapt( inputFastqFiles, asMatePairs=FALSE, forceMatePairs=FALSE,
						cutadaptProgram=cutadaptProgram, kmer.size=45, verbose=verbose)
			fastqFiles <- paste( fastqFiles, "trimmed", sep=".")
		}
	
		# step 2:  create contigs of anything that may be HLA gene reads
		spades.mode <- match.arg( spades.mode)
		pipe.SpadesContigs( sampleID, fastqSource=fastqFiles, folderName="HLA", 
				spades.mode=spades.mode, kmerSizes=kmerSizes, spades.args=spades.args, 
				spades.path=spades.path, keyword=keyword, makePep=T)

		# cleanup SPAdes temp files
		cleanupSpadesFiles( path=spades.output.path, verbose=F)
	} else {
		cat( "\n\nUsing previously calculated Spades HLA contigs..\n")
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

	# step 3:  Throw those contigs as peptides against the given set of full length HLA gene proteins
	cat( "\n\nSearching Contigs for HLA constructs..")
	proteinAns <- bestSpadesProteins( sampleID, outpath=spades.output.path, proteinFastaFile=hlaFastaFile,
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
		cat( "\nWrote", NPROT, "HLA proteins as FASTA to: ", basename(proteinsFastaFile))
	}

	# at this point, we have a small chance that none look like HLA genes
	if ( ! NPROT) {
		cat( "\nNo peptides look like HLA proteins..")
		write.table( proteinAns, proteinsTextFile, sep="\t", quote=F, row.names=F)
		return( invisible( proteinAns))
	}
	return( invisible( proteinAns))
}


`pipe.HLA.ContigSummary` <- function( sampleIDset, optionsFile="Options.txt", 
						results.path=NULL, keyword="HLA") {

	if ( is.null( results.path)) results.path <- getOptionValue( optionsFile, "results.path", 
				notfound=".", verbose=F)

	# gather up all the HLA contig files for the named samples
	spades.path <- file.path( results.path, "SpadesContigs")
	N <- length( sampleIDset)
	cat( "\nChecking", N, "folders of Spades results..")

	proteinTextFiles <- proteinFastaFiles <- vector()
	for ( i in 1:N) {
		thisSID <- sampleIDset[i]
		this.path <- file.path( spades.path, thisSID, keyword)
		if ( ! file.exists(this.path)) next
		proteinTextFiles[i] <- file.path( this.path, paste( thisSID, keyword, "BestProteinHits.txt", sep="."))
		proteinFastaFiles[i] <- file.path( this.path, paste( thisSID, keyword, "Proteins.fasta", sep="."))
	}

	cat( "\nGathering Best Protein Hits..\n")
	protDF <- data.frame()
	proteinFile <- "HLA.Contigs.BestProteinHits.txt"
	for ( f in proteinTextFiles) {
		if ( is.na(f)) next
		if ( ! file.exists(f)) next
		sml <- read.delim( f, as.is=T)
		cat( "\r", basename(f), " ", nrow(sml))
		if ( ! nrow(sml)) next
		protDF <- rbind( protDF, sml)
	}
	write.table( protDF, proteinFile, sep="\t", quote=F, row.names=F)
	cat( "\nWrote: ", proteinFile, " \tN_Protein_Hits: ", nrow(protDF))

	cat( "\nGathering Protein FASTA..\n")
	desc <- seq <- vector()
	proteinFile <- "HLA.Contigs.Proteins.fasta"
	for ( f in proteinFastaFiles) {
		if ( is.na(f)) next
		if ( ! file.exists(f)) next
		sml <- loadFasta( f, short=F, verbose=F)
		cat( "\r", basename(f), " ", length(sml$desc))
		if ( ! length(sml$desc)) next
		desc <- c( desc, sml$desc)
		seq <- c( seq, sml$seq)
	}
	writeFasta( as.Fasta( desc, seq), proteinFile, line=100)
	cat( "\nWrote: ", proteinFile, " \tN_Proteins: ", length(desc))
	return( N)
}

