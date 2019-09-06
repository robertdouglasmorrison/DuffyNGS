# pipe.BAMtoProteins.R -- get the consensus proteins directly from the BAM file results


`pipe.BAMtoProteins` <- function( sampleID, geneIDset=NULL, annotationFile="Annotation.txt",
				optionsFile="Options.txt", speciesID=getCurrentSpecies(), results.path=NULL,
				SNP.only=TRUE, minReadCalls=NULL, minPercentSNP=NULL, verbose=TRUE) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}

	# make sure we have the BAM file already sorted
	bamfile <- paste( sampleID, "genomic.bam", sep=".")
	bamfile <- file.path( results.path, "align", bamfile)
	sortedbamfile <- BAM.verifySorted( bamfile, index=TRUE)
	if ( is.null( sortedbamfile)) return(NULL)

	if ( speciesID != getCurrentSpecies()) setCurrentSpecies( speciesID)
	geneMap <- getCurrentGeneMap()
	cdsMap <- getCurrentCdsMap()
	allGenes <- sort( unique( cdsMap$GENE_ID))
	if ( is.null( geneIDset)) {
		geneIDset <- allGenes
	} else {
		geneIDset <- intersect( geneIDset, allGenes)
	}

	require( Biostrings)


	`getProteinOneGene` <- function( gid) {

		where <- match( gid, geneMap$GENE_ID)
		myStart <- geneMap$POSITION[where]
		myStop <- geneMap$END[where]
		myStrand <- geneMap$STRAND[where]
		mySID <- geneMap$SEQ_ID[where]

		ans <- pipe.ConsensusBaseCalls( sampleID, geneID=gid, seqID=mySID, start=myStart, stop=myStop, 
					annotationFile=annotationFile,
					optionsFile=optionsFile, results.path=results.path,
					aaToo=TRUE, noReadCalls="genomic", as.cDNA=TRUE, SNP.only=SNP.only, 
					minReadCalls=minReadCalls, minPercentSNP=minPercentSNP, verbose=FALSE)

		# by default, indels can hose the translation into proteins
		# do we want to prevent that?
		if ( SNP.only) {
			# indels will be detectable as elements that are not a single character long.
			# to prevent them from causing trouble, let's find them and replace with the reference
			dna <- ans$dna.consensus
			isIndel <- which( nchar(dna) != 1)
			if ( length( isIndel)) {
				refDNA <- ans$ref[ isIndel]
				refDNA[ nchar(refDNA) != 1] <- "N"
				dna[ isIndel] <- refDNA
			}
			dna <- paste( dna, collapse="")
			myProt <- DNAtoAA( dna, clipAtStop=F, readingFrame=1)
		} else {
			myAA <- ans$aa.consensus
			myProt <- paste( myAA, collapse="")
		}

		cat( "\r", gid, nchar(myProt))
		return( myProt)
	}
	

	cat( "\nExtracting proteins from BAM consensus: ", length(geneIDset), "\n")
	ans <- multicore.lapply( geneIDset, FUN=getProteinOneGene)
	
	allProts <- unlist( ans)
	names( allProts) <- geneIDset
	return( allProts)
}

