# pipe.BAMgenome.R -- get the consensus genome of a sample directly from the BAM file results


# convert a BAM file into a FASTA of the complete genome

`pipe.BAMgenome` <- function( sampleID, seqIDset=NULL, annotationFile="Annotation.txt",
				optionsFile="Options.txt", speciesID=getCurrentSpecies(), results.path=NULL,
				noReadCalls=NULL, chunk.size=100000, verbose=TRUE) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	if ( is.null( noReadCalls)) {
		dataType <- getAnnotationValue( annotationFile, key=sampleID, columnArg="DataType", notfound="DNA-seq", verbose=verbose)
		if (dataType == "DNA-seq") noReadCalls <- "blank"
		if (dataType == "RNA-seq") {
			cat( "\nWarning: extracting genomid DNA consensus from RNA-seq data is unwise.")
			cat( "\n         Default is to use reference genome for regions with no RNA-seq reads.")
			noReadCalls <- "genomic"
		}
	}
	if ( is.null(noReadCalls) || !(noReadCalls %in% c("blank","genomic"))) {
		cat( "\nArgument 'noReadCalls' must be one of 'blank' or 'genomic'")
		stop( "See command 'pipe.BAMgenome()'")
	}

	# make sure we have the BAM file already sorted
	bamfile <- paste( sampleID, "genomic.bam", sep=".")
	bamfile <- file.path( results.path, "align", bamfile)
	sortedbamfile <- BAM.verifySorted( bamfile, index=TRUE)
	if ( is.null( sortedbamfile)) return(NULL)

	if ( speciesID != getCurrentSpecies()) setCurrentSpecies( speciesID)
	seqMap <- getCurrentSeqMap()
	allSeqs <- sort( unique( seqMap$SEQ_ID))
	if ( is.null( seqIDset)) {
		seqIDset <- allSeqs
	} else {
		seqIDset <- intersect( seqIDset, allSeqs)
	}
	if ( length( seqIDset) < 1) {
		cat( "\nNo Chromosomes selected")
		return(NULL)
	}

	require( Biostrings)
	if (version$major == "4" && as.numeric( version$minor) >= 4) require( pwalign)

	# get the genome into vectors of bases
	genomicFastaFile <- getOptionValue( optionsFile, "genomicFastaFile", verbose=verbose)
	fa <- loadFasta( genomicFastaFile, verbose=verbose)
	baseVectors <- strsplit( fa$seq, split="")
	names(baseVectors) <- fa$desc

	# set up where the result FASTA will be written to
	fa.path <- file.path( results.path, "ConsensusGenomes", sampleID)
	if ( ! file.exists( fa.path)) dir.create( fa.path, recursive=T)
	fa.file <- file.path( fa.path, paste( "BAM.Genome", sampleID, "fasta", sep="."))

	if (verbose) cat( "\nExtracting genome from BAM consensus:  N_Chromosomes =", length(seqIDset), "\n")
	
	# do and write each chromosome one at a time, to save memory, etc.
	PASTE <- base::paste
	
	nBaseOut <- 0
	for ( sid in seqIDset) {
		where <- MATCH( sid, seqMap$SEQ_ID)
		seqStop <- seqMap$LENGTH[where]		
		myBaseVecPt <- match( sid, names(baseVectors))

		# get the pileups for the entire chromosome, but do it in chunks to save memory, etc.
		cat( "\nGetting pileups for: ", sid, "  N_Bases =", seqStop)
		
		bigDNA <- ""
		curStop <- 0
		repeat {
			curStart <- curStop + 1
			if ( curStart > seqStop) break
			curStop <- curStop + chunk.size
			if ( curStop > seqStop) curStop <- seqStop
			
			# grab this chunk
			ans <- pipe.ConsensusBaseCalls( sampleID, geneID=NULL, seqID=sid, start=curStart, stop=curStop, 
						annotationFile=annotationFile, genomicFastaFile=genomicFastaFile,
						genomicVector=baseVectors[[myBaseVecPt]],
						optionsFile=optionsFile, results.path=results.path, noReadCalls=noReadCalls,
						aaToo=FALSE, as.cDNA=FALSE, best.frame=FALSE, SNP.only=FALSE, 
						minReadCalls=NULL, minPercentSNP=NULL, verbose=verbose)
					
			# turn it from vector back to one giant string
			curDNA <- PASTE( ans$dna.consensus, collapse="")
			# and join it to what we have already
			bigDNA <- PASTE( bigDNA, curDNA, sep="")
			cat(".")
			rm( ans)
		}

		# write it out
		out <- as.Fasta( "desc"=paste( sid, sampleID, sep="_"), "seq"=bigDNA)
		doAPPEND <- (nBaseOut > 0)
		writeFasta( out, fa.file, line.width=100, append=doAPPEND)
		nBaseOut <- nBaseOut + nchar( bigDNA)
	}
	if (verbose) cat( "\nWrote Fasta Genome file: ", fa.file, "  N_Bases =", nBaseOut, "\n")
}

