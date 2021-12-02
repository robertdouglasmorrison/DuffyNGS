# pipe.RnaEditingSearch.R -- compare a DNA and RNA sample pair to find differences

`pipe.RnaEditingSearch` <- function( dnaSampleID, rnaSampleID, geneIDset=NULL, min.depth=5, exon.edge.trim=0,
				results.path=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt", 
				speciesID=getCurrentSpecies(), verbose=TRUE) {

	# get options that we need
	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	bam.path <- file.path( results.path, "align")

	# make sure our 2 samples are what we expect
	dataType1 <- getAnnotationValue( annotationFile, key=dnaSampleID, columnArg="DataType", notfound="DNA-seq", verbose=verbose)
	dataType2 <- getAnnotationValue( annotationFile, key=rnaSampleID, columnArg="DataType", notfound="RNA-seq", verbose=verbose)
	if ( dataType1 != "DNA-seq") stop( "First SampleID required to be genomic DNA sample.")
	if ( dataType2 != "RNA-seq") stop( "Second SampleID required to be messenger RNA sample.")

	# make sure we have the BAM files already sorted
	bamfile1 <- paste( dnaSampleID, "genomic.bam", sep=".")
	bamfile1 <- file.path( results.path, "align", bamfile1)
	sortedbamfile1 <- BAM.verifySorted( bamfile1, index=TRUE)
	if ( is.null( sortedbamfile1)) return(NULL)
	bamfile2 <- paste( rnaSampleID, "genomic.bam", sep=".")
	bamfile2 <- file.path( results.path, "align", bamfile2)
	sortedbamfile2 <- BAM.verifySorted( bamfile2, index=TRUE)
	if ( is.null( sortedbamfile2)) return(NULL)

	# get the set of genes to investigate
	if ( speciesID != getCurrentSpecies()) setCurrentSpecies( speciesID)
	geneMap <- getCurrentGeneMap()
	cdsMap <- getCurrentCdsMap()
	allGenes <- sort( unique( cdsMap$GENE_ID))
	if ( is.null( geneIDset)) {
		geneIDset <- allGenes
	} else {
		geneIDset <- intersect( geneIDset, allGenes)
	}
	if ( length( geneIDset) < 1) {
		cat( "\nNo Genes selected")
		return(NULL)
	}

	# get the genome into vectors of bases
	genomicFastaFile <- getOptionValue( optT, "genomicFastaFile", verbose=FALSE)
	fa <- loadFasta( genomicFastaFile, verbose=FALSE)
	baseVectors <- strsplit( fa$seq, split="")
	names(baseVectors) <- fa$desc

	MATCH <- base::match
	NCHAR <- base::nchar
	PASTE <- base::paste
	WHICH <- base::which


	`getSequenceOneGene` <- function( gid, sid) {

		where <- MATCH( gid, geneMap$GENE_ID)
		myStart <- geneMap$POSITION[where]
		myStop <- geneMap$END[where]
		myStrand <- geneMap$STRAND[where]
		mySID <- geneMap$SEQ_ID[where]
		myBaseVecPt <- match( mySID, names(baseVectors))
		
		ans <- pipe.ConsensusBaseCalls( sid, geneID=gid, seqID=mySID, start=myStart, stop=myStop, 
					annotationFile=annotationFile, genomicFastaFile=genomicFastaFile,
					genomicVector=baseVectors[[myBaseVecPt]],
					optionsFile=optionsFile, results.path=results.path, noReadCalls="genomic",
					aaToo=TRUE, as.cDNA=TRUE, best.frame=FALSE, SNP.only=FALSE, 
					minReadCalls=NULL, minPercentSNP=NULL, verbose=FALSE)
		return( ans)
	}
	
	require( Biostrings)
	data( BLOSUM62)

	if (verbose) cat( "\nExtracting cDNA from BAM consensus pileups:  N_Genes =", length(geneIDset), "\n")

	out <- data.frame()
	
	for ( i in 1:length(geneIDset)) {
		g <- geneIDset[i]
		cat( "\r", i, g)
		ans1 <<- getSequenceOneGene( g, dnaSampleID)
		ans2 <<- getSequenceOneGene( g, rnaSampleID)

		# decide which parts of the gene are comparable, by read coverage, etc.
		cdnaDepth1 <- apply( ans1$callsMatrix, 1, sum, na.rm=T)
		cdnaDepth2 <- apply( ans2$callsMatrix, 1, sum, na.rm=T)
		aaDepth1 <- sapply( seq( 1, length(cdnaDepth1), by=3), function(x) { y <- min( x+2, length(cdnaDepth1)); min( cdnaDepth1[x:y])})
		aaDepth2 <- sapply( seq( 1, length(cdnaDepth2), by=3), function(x) { y <- min( x+2, length(cdnaDepth2)); min( cdnaDepth2[x:y])})

		# get the protein call
		aa1 <- paste( ans1$callsTable$AA, collapse="")
		aa2 <- paste( ans2$callsTable$AA, collapse="")

		# do we see any differences? Use pairwise alignment
		pa <- pairwiseAlignment( aa1, aa2, type="global", substitutionMatrix=BLOSUM62)
		aln1 <- as.character( alignedPattern( pa))
		aln2 <- as.character( alignedSubject( pa))
		alnV1 <- strsplit( aln1, split="")[[1]]
		alnV2 <- strsplit( aln2, split="")[[1]]
		diffs <- which( alnV1 != alnV2)
		if ( ! length(diffs)) next

		sml <- data.frame( "AA_POSITION"=diffs, "DNA_CALL"=alnV1[diffs], "RNA_CALL"=alnV2[diffs], 
					"DNA_DEPTH"=aaDepth1[diffs], "RNA_DEPTH"=aaDepth2[diffs], 
					stringsAsFactors=F)
		# don't let very low coverage sites get kept
		drops <- which( sml$DNA_DEPTH < min.depth | sml$RNA_DEPTH < min.depth)
		if ( length( drops)) sml <- sml[ -drops, ]
		if ( ! nrow(sml)) next

		cat( "\nFound Differences: ", i, g, "  N=", length(diffs), "\n")
		print( head( sml))
		cat( "\n")

		out <- rbind( out, cbind( "GENE_ID"=g, sml, stringsAsFactors=F))
		break
	}

	if ( nrow(out)) rownames(out) <- 1:nrow(out)
	return( out)
}

