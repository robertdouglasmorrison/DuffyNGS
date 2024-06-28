# pipe.HLA.Expression.R

# turn RNA-seq data into HLA allele specific abundance data.

# Uses the HLA.ProteinCalls results about alleles, to build a custom Bowtie2 target 'genome'
# and then aligns any reads that hit HLA genes against this target, to get allele-specific
# expression


`pipe.HLA.Expression` <- function( sampleIDset, hlaFastaFile, annotationFile="Annotation.txt", optionsFile="Options.txt",
				results.path=NULL, verbose=TRUE) {

	if ( is.null( results.path)) results.path <- getOptionValue( optionsFile, "results.path", 
				notfound=".", verbose=F)

	# make sure the reference file of HLA genes is readable
	if ( ! file.exists( hlaFastaFile)) {
		cat( "\nError: 'hlaFastaFile' of reference DNA not found: ", hlaFastaFile)
		return(NULL)
	}
	hlaFA <- loadFasta( hlaFastaFile, verbose=verbose)

	# path for all inputs and results
	HLAresults.path <- file.path( results.path, "HLA.ProteinCalls")
	if ( ! file.exists( HLAresults.path)) {
		cat( "\nError: Unable to find needed HLA allele calls directory: ", HLAresults.path)
		return(NULL)
	}
	bowtieFastaInputFile <- "Bowtie2.HLA.Allele.Inputs.fasta"
	bowtieTargetFile <- "Bowtie2.HLA.Allele.Target_idx"

	# get the universe of HLA genes
	if ( getCurrentSpecies() != "Hs_grc") stop( "HLA.Expression is only defind for human samples")
	gmap <- getCurrentGeneMap()
	hlaMap <- gmap[ grep( "^HLA\\-", gmap$GENE_ID), ]
	hlagenes <- sort( unique( hlaMap$GENE_ID))

	# get the union of all alleles seen by all these samples
	nSamples <- length( sampleIDset)
	if ( nSamples < 2) cat( "\nWarning:  expected 2 or more sample IDs, to create a matrix of allele expression.")
	allelesSeen <- vector()
	for (sid in sampleIDset) {
		allelesFile <- file.path( HLAresults.path, sid, paste( sid, "Merged.HLA.Calls.csv", sep="."))
		if ( ! file.exists( allelesFile)) {
			cat( "\nWarning: file of allele calls not found: ", basename(allelesFile), "  Skipping..")
			next
		}
		tmp <- read.csv( allelesFile, as.is=T)
		allelesSeen <- c( allelesSeen, tmp$Allele)
	}
	allelesSeen <- sort( unique( allelesSeen))
	nAlleles <- length(allelesSeen)
	if ( verbose) cat( "\nTotal number of HLA alleles: ", nAlleles)
	
	# make a Bowtie2 index from exactly those
	if (verbose) cat( "\nBuilding custom HLA Bowtie2 target index..")
	where <- match( allelesSeen, hlaFA$desc, nomatch=0)
	alleleFasta <- as.Fasta( hlaFA$desc[where], hlaFA$seq[where])
	writeFasta( alleleFasta, bowtieFastaInputFile, line=100)
	buildCmd <- buildBowtie2BuildCommandLine( bowtieFastaInputFile, bowtieTargetFile, optionsFile="Options.txt", verbose=verbose)
	catch.system( command=buildCmd, wait=T)
	file.delete( bowtieFastaInputFile)
	
	# before we visit each sample, set up the storage and expression scaling info we will need
	exprM <- matrix( NA, nrow=nAlleles, ncol=nSamples)
	colnames(exprM) <- sampleIDset
	rownames(exprM) <- allelesSeen
	nBasesPerAllele <- nchar( alleleFasta$seq)
	
	for ( j in 1:nSamples) {
		sid <- sampleIDset[j]
		alignedReadsFile <- file.path( results.path, "fastq", paste( sid, "HLA.Hits.fastq.gz", sep="."))

		# step 1:   gather all aligned reads that land in/near the HLA loci
		if ( ! file.exists( alignedReadsFile)) {
			cat( "\nGathering HLA Gene aligned reads for sample: ", sid, "\n")
			pipe.GatherGeneAlignments( sid, genes=hlagenes, asFASTQ=TRUE, mode="best.one", fastq.keyword="HLA.Hits")
		}
	
		# step 2:  call Bowtie2 with these raw reads
		if (verbose) cat( "\nCalling Bowtie for sample: ", sid)
		# put the sample name in these temp files, to prevent collisions
		bamFile <- paste( sid, "HLA.AlleleHits.bam", sep=".")
		bowtieAns <- fastqToBAM( inputFastqFile=alignedReadsFile, outputFile=bamFile, sampleID=sid, 
					optionsFile=optionsFile, annotationFile=annotationFile, 
					alignIndex=bowtieTargetFile, index.path=".", noHitsFile=NULL, verbose=F)
		nReadsAligned <- bowtieAns$UniqueReads + bowtieAns$MultiReads
		
		# step 3:  extract the read counts per allele.   Since we treated each allele as its own chromosome,
		# this is just count by SEQ_ID.
		countsAns <- countReadsBySequence( bamFile, readBufferSize=100000, verbose=F)
		
		# step 4:  normalize and store
		tpmUnits <- tpm( as.numeric( countsAns), geneLen=nBasesPerAllele)
		wh <- match( allelesSeen, names(countsAns), nomatch=0)
		exprM[ wh > 0, j] <- tmpUnits[ wh]
		
		# step 5: clean up after each sample
		file.delete( bamFile)
	}
	# step 6: clean up after all samples
	file.delete( c( bowtieFastaInputFile, bowtieTargetFile))

	out <- data.frame( "Allele"=allelesSeen, round(exprM, digits=3), stringsAsFactors=F)
	rownames(out) <- 1:nrow(out)
	return(out)
}
