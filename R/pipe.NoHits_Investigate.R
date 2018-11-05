# pipe.NoHits_Investigate.R


`pipe.NoHits_Investigate` <- function( sampleID, annotationFile="Annotation.txt",
		optionsFile="Options.txt", nUSRkeep=1000000, reload=FALSE, verbose=!interactive()) {

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\n\nStarting 'No Hits Investigate' for Sample:     ", sampleID, "\n\n")
	}

	optT <- readOptionsTable( optionsFile)
	setCurrentTarget( getOptionValue( optT, "targetID", notfound="HsPf"))
	resultsPath <- getOptionValue( optT, "results.path", notfound=".", verbose=verbose)
	fastqPath <- getOptionValue( optT, "fastqData.path", verbose=verbose)
	trim5 <- as.numeric( getOptionValue( optT, "trim5", notfound=0)) 
	trim3 <- as.numeric( getOptionValue( optT, "trim3", notfound=0))
	inputFastqFile <- getAnnotationValue( annotationFile, sampleID, "Filename")

	# allow comma sep files for paired reads
	inputFastqFileSet <- strsplit( inputFastqFile, split=", *")[[1]]

	# get the sampleID and true filename for this pair
	inputFastqFile <- inputFastqFileSet[ 1]
	inputFastqFile <- file.path( fastqPath, inputFastqFile)
	inputFastqFile <- allowCompressedFileName( inputFastqFile)

	startT <- proc.time()
	
	# get the file of noHit reads
	infile <- paste( sampleID, "noHits", "fastq", sep=".") 
	infile <- file.path( resultsPath, "fastq", infile)
	infile <- allowCompressedFileName( infile)
	if ( ! file.exists( infile)) {
		notgenomicfile <- paste( sampleID, "not.genomic", "fastq", sep=".") 
		notgenomicfile <- file.path( resultsPath, "fastq", notgenomicfile)
		cat( "\nFile not found:  ", infile, "\nTrying: ", notgenomicfile)
		infile <- notgenomicfile
	}
	nNoHits <- getFileLineCount( infile, sampleID) / 4

	# also get the number of original raw reads...
	nRawReads <- getFileLineCount( inputFastqFile, sampleID) / 4

	outTxt <- vector()
	outTxt[1] <- paste( "'NoHit' reads: ", prettyNum( nNoHits, big.mark=","), as.percent(1), 
			as.percent( nNoHits, big.value=nRawReads), sep="\t")
	outTxt[2] <- " "
	nOut <- length( outTxt)
	
	#ans <- noHit.contaminants( sample, annotationFile, optionsFile, nRawReads=nRawReads) 
	#outTxt <- base::append( outTxt, ans$text)
	#nOut <- length( outTxt)
	#nHitContaminants <- ans$nHits
	nHitContaminants <- 0

	# set up for adapter tail and poly A search...
	ans <- USR_setup( infile, sampleID, resultsPath, Nkeep=nUSRkeep, reload=reload, 
			trim5=trim5, trim3=trim3)
	nPolyNdrops <- ans$nPolyN
	USRfile <- ans$USR_File

	# now look for poly nucleotides
	outTxt[nOut+1] <- " "
	nOut <- nOut + 1

	ans <- noHit.polyBase( sampleID, annotationFile, optionsFile, nPolyN=nPolyNdrops,
				nRawReads=nRawReads)
	outTxt <- base::append( outTxt, ans$text)
	nOut <- length( outTxt)
	nHitPoly <- ans$nHits

	# now look for adapter tails
	ans <- noHit.adapterTails( sampleID, annotationFile, optionsFile,
				nRawReads=nRawReads)
	outTxt <- base::append( outTxt, ans$text)
	nOut <- length( outTxt)
	nHitAdapters <- ans$nHits
	adapterCounts <- ans$adapterCounts

	# stash adapter tail info with the USR, upate and then we're done with USR data
	USR_AdapterHits <<- adapterCounts
	saveUSRcontext( USRfile)
	USR_cleanup()

	# other Plasmodium variants
	nHitVariants <- 0
	#if ( any( regexpr( "Pf", getCurrentTarget()) > 0)) {
#
#		ans <- noHit.PfVariants( sample, annotationFile, optionsFile,
#					nRawReads=nRawReads)
#		outTxt <- base::append( outTxt, ans$text)
#		nOut <- length( outTxt)
#		nHitVariants <- ans$nHits
#	}


	# all done, report it...
	ttlHits <- sum( nHitContaminants, na.rm=T) + nHitPoly + nHitAdapters + sum( nHitVariants, na.rm=T)
	outTxt[nOut+1] <- " "
	outTxt[nOut+2] <- paste( "Total 'NoHit' Hits:\t", ttlHits, 
				"\t", as.percent( ttlHits, big.value=nNoHits),
				"\t", as.percent( ttlHits, big.value=nRawReads))
	nOut <- nOut + 2
	cat( "\n\n", outTxt[nOut])

	summaryFile <- file.path( resultsPath, "summary", paste( sampleID, "noHits", "Summary.txt", sep="."))
	writeLines( outTxt, con=summaryFile)

	if ( verbose) cat( "\n...Done.\n")

	cat( "\n\nFinished 'No Hits Investigation' for Sample:     ", sampleID, "\n\n")

	return()
}


`noHit.contaminants` <- function( sampleID, annotationFile="Annotation.txt",
		optionsFile="Options.txt", nRawReads=1000000) {

	# get the file of noHit reads
	resultsPath <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=FALSE)
	infile <- paste( sampleID, "noHits", "fastq", sep=".") 
	infile <- file.path( resultsPath, "fastq", infile)
	nNoHits <- getFileLineCount( infile, sampleID) / 4

	outTxt <- vector()
	outTxt[1] <- " "
	outTxt[2] <- "Contaminants:"
	cat( "\n\n", outTxt[2])
	nOut <- 2
	
	# see how many hit various contaminants
	targets <- c( "Eco.genomic_idx", "mycoplasma", "Sty.genomic_idx")
	names(targets) <- c( "E.coli      ", "Mycoplasma  ", "S.typhimurium")
	doGeneIDs <- c( TRUE, FALSE, TRUE)
	targetSpeciesID <- c( "Ecoli", "", "Styphi")
	nHitTargets <- rep( 0, time=length( targets))
	outfile <- paste( sampleID, "tempAlign", "bam", sep=".")

	for( i in 1:length( targets)) {
		ans <- fastqToBAM( infile, outfile, m=NULL, k=1, sampleID=sampleID, 
				alignIndex=targets[i], alignPolicy=" ", verbose=FALSE)
		nHitTargets[i] <- ans$ReadsAligned
		outTxt[nOut+i] <- paste( "    ", names(targets)[i], "\t", nHitTargets[i], 
					"\t", as.percent( nHitTargets[i], big.value=nNoHits),
					"\t", as.percent( nHitTargets[i], big.value=nRawReads))
		cat( "\n", outTxt[nOut+i])

		if ( doGeneIDs[i]) {
			savefile <- paste( sampleID,"noHits", targets[i], sep=".")
			savefile <- sub( "_idx", ".align", savefile)
			savefile <- file.path( resultsPath, "align", savefile)
			saveSpecies <- getCurrentSpecies()
			setCurrentSpecies( targetSpeciesID[i])
			alignToAlignID( outfile, savefile, sampleID=sampleID, verbose=F)
			tmp <- readAlignmentFile( savefile, sampleID=sampleID, verbose=F)
			topGenes <- base::table( tmp$GENE_ID)
			if ( length( topGenes) > 0) {
				gprod <- gene2Product( names(topGenes))
				sml <- data.frame( "GENE_ID"=names(topGenes), "HITS"=as.vector(topGenes), 
						"PRODUCT"=gprod)
				ord <- base::order( sml$HITS, decreasing=T)
				sml <- sml[ ord, ]
				rownames(sml) <- 1:nrow(sml)
				N <- min( nrow(sml), 10)
				cat( "\nTop hits:     ", names(targets)[i], "\n")
				print( sml[ 1:N, ])
			}
			setCurrentSpecies( saveSpecies)
		}

	}
	nOut <- nOut + length( targets)
	return( list( "nHits"=nHitTargets, "text"=outTxt))
}


`noHit.polyBase` <- function( sampleID, annotationFile="Annotation.txt",
		optionsFile="Options.txt", nRawReads=1000000, nPolyN=0) {

	# get the file of noHit reads
	resultsPath <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=FALSE)
	infile <- paste( sampleID, "noHits", "fastq", sep=".") 
	infile <- file.path( resultsPath, "fastq", infile)
	nNoHits <- getFileLineCount( infile, sampleID) / 4

	outTxt <- vector()
	outTxt[1] <- " "
	outTxt[2] <- "Poly Base Search:"
	cat( "\n\n", outTxt[2])
	nOut <- 2
	
	ans <- countPolyBases()
	nAs <- ans$nAs
	nCs <- ans$nCs
	nGs <- ans$nGs
	nTs <- ans$nTs
	nNs <- ans$nNs + nPolyN

	outTxt[nOut+1] <- paste( "    ", "Poly 'N'", "\t", nNs, "\t", as.percent( nNs, big.value=nNoHits),
						"\t", as.percent( nNs, big.value=nRawReads))
	outTxt[nOut+2] <- paste( "    ", "Poly 'A'", "\t", nAs, "\t", as.percent( nAs, big.value=nNoHits),
						"\t", as.percent( nAs, big.value=nRawReads))
	outTxt[nOut+3] <- paste( "    ", "Poly 'C'", "\t", nCs, "\t", as.percent( nCs, big.value=nNoHits),
						"\t", as.percent( nCs, big.value=nRawReads))
	outTxt[nOut+4] <- paste( "    ", "Poly 'G'", "\t", nGs, "\t", as.percent( nGs, big.value=nNoHits),
						"\t", as.percent( nGs, big.value=nRawReads))
	outTxt[nOut+5] <- paste( "    ", "Poly 'T'", "\t", nTs, "\t", as.percent( nTs, big.value=nNoHits),
						"\t", as.percent( nTs, big.value=nRawReads))
	cat( "\n", outTxt[nOut+1])
	cat( "\n", outTxt[nOut+2])
	cat( "\n", outTxt[nOut+3])
	cat( "\n", outTxt[nOut+4])
	cat( "\n", outTxt[nOut+5])
	nHitPoly <- nNs + nAs + nCs + nGs + nTs
	nOut <- nOut + 5

	return( list( "nHits"=nHitPoly, "text"=outTxt))
}


`noHit.adapterTails` <- function( sampleID, annotationFile="Annotation.txt",
		optionsFile="Options.txt", logFreq=100000, Nquit=1000000, 
		nRawReads=1000000, verbose=TRUE) {

	# get the file of noHit reads
	optT <- readOptionsTable( optionsFile)
	#resultsPath <- getOptionValue( optT, "results.path", notfound=".", verbose=FALSE)
	#infile <- paste( sampleID, "noHits", "fastq", sep=".") 
	#infile <- file.path( resultsPath, "fastq", infile)
	#nNoHits <- getFileLineCount( infile, sampleID) / 4

	# get the adapters
	adapter1 <- getOptionValue( optT, "ForwardAdapter", 
			notfound="AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT")
	adapter2 <- getOptionValue( optT, "ReverseAdapter", 
			notfound="CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT")
	minOverlap <- as.numeric( getOptionValue( optT, "Adapter.minimumOverlap", notfound="20"))
	adapters <- list( "rc(FwdAdapt)"=myReverseComplement(adapter1),
			"rc(RevAdapt)"=myReverseComplement(adapter2))

	outTxt <- vector()
	outTxt[1] <- " "
	outTxt[2] <- "Adapter Tail Search:"
	nOut <- 2
	if (verbose) cat( "\n\n", outTxt[nOut])
	ans <- countAdapterHits( adapters, minOverlap=minOverlap, logFreq=logFreq, Nquit=Nquit)
	nhits <- ans$nTails
	nTested <- ans$nTested

	for ( j in 1:2) {
		outTxt[nOut+j] <- paste( "    ", names(adapters)[j], "\t", nhits[j], 
					"\t", as.percent( nhits[j], big.value=nTested),
					"\t", as.percent( nhits[j], big.value=nRawReads))
		if (verbose) cat( "\n", outTxt[nOut+j])
	}
	return( list( "adapterCounts"=nhits, "nHits"=sum(nhits), "text"=outTxt))
}


`noHit.PfVariants` <- function( sampleID, annotationFile="Annotation.txt",
		optionsFile="Options.txt", nRawReads=1000000) {

	# get the file of noHit reads
	resultsPath <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=FALSE)
	infile <- paste( sampleID, "noHits", "fastq", sep=".") 
	infile <- file.path( resultsPath, "fastq", infile)
	nNoHits <- getFileLineCount( infile, sampleID) / 4

	outTxt <- vector()
	outTxt[1] <- " "
	outTxt[2] <- "Non-3D7 Strains and Isolates: "
	nOut <- 2
	cat( "\n\n", outTxt[nOut])
	species <- c( "mastercard_idx")
	names(species) <- c( "Other_Strains")
	nHitVariants <- rep( 0, time=length( species))
	outfile <- paste( sampleID, "tempAlign", "tmp", sep=".")

	for( i in 1:length( species)) {
		ans <- fastqToAlign( infile, outfile, m=NULL, k=1, sampleID=sampleID, 
				alignIndex=species[i], 
				alignPolicy=" -n 2 -l 26 -e 75 ", verbose=FALSE)
		nHitVariants[i] <- ans$ReadsAligned
		outTxt[nOut+i] <- paste( "    ", names(species)[i], "\t", nHitVariants[i], 
					"\t", as.percent( nHitVariants[i], big.value=nNoHits),
					"\t", as.percent( nHitVariants[i], big.value=nRawReads))
		cat( "\n", outTxt[nOut+i])

	}
	nOut <- nOut + length( species)
	file.delete( outfile)
	return( list( "nHits"=nHitVariants, "text"=outTxt))
}

