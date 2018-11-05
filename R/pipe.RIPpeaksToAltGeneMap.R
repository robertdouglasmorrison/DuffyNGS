# pipe.RIPpeaksToAltGeneMap.R

# combine multiple samples RIP Peak results, to create one Alternate Gene Map of RIP features
# for use in Differential Expression workflows

`pipe.RIPpeaksToAltGeneMap` <- function( sampleIDset, annotationFile="Annotation.txt", 
				optionsFile="Options.txt", speciesID=NULL, results.path=NULL,
				max.pvalue=0.05) {

	annT <- readAnnotationTable( annotationFile)
	optT <- readOptionsTable( optionsFile)
	if ( is.null(speciesID)) speciesID <- getCurrentSpecies()
	setCurrentSpecies( speciesID)
	speciesPrefix <- getCurrentSpeciesFilePrefix()

	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	rip.path <- file.path( results.path, "RIPpeaks")

	allPeaks <- data.frame()

	for ( sampleID in sampleIDset) {

		fileInPeaks <- paste( sampleID, speciesPrefix, "RIPpeaks.txt", sep=".")
		fileInPeaks <- file.path( rip.path, sampleID, fileInPeaks)
		if ( ! file.exists( fileInPeaks)) {
			cat( "\nExisting RIP Peaks file not found: ", fileInPeaks)
			cat( "\nSkipping...")
			next
		}
		peaks <- read.delim( fileInPeaks, as.is=T)
		cat( "\n  ", basename(fileInPeaks), "   N=", nrow(peaks))
		peaks <- subset( peaks, P_Value <= max.pvalue)
		cat( "   N_Significant=", nrow(peaks))
		if ( ! nrow( peaks)) next
		allPeaks <- rbind( allPeaks, data.frame( "SampleID"=sampleID, peaks, stringsAsFactors=F))
		cat( "  Total=", nrow( allPeaks))
	}

	# with all the peaks together, we will merge to make one common set
	ord <- order( allPeaks$Height, decreasing=T)
	peaksIn <- allPeaks[ ord, ]
	NP <- nrow( peaksIn)
	isAvailable <- rep.int( TRUE, NP)
	peaksOut <- data.frame()
	cat( "\nVisiting", NP, "peaks from", length(sampleIDset), "samples to build 'AltGeneMap' of RIP peaks..\n")

	# visit each peak, and find aall that overlap it (within reason...)
	for ( i in 1:NP) {
		if ( ! isAvailable[i]) next
		peakNow <- peaksIn[ i, ]

		# only consider a small window around this peak
		myCenter <- peakNow$Center[1]
		myLeftEdge <- peakNow$Start[1]
		myRightEdge <- peakNow$Stop[1]
		toVisit <- which( isAvailable == TRUE & peaksIn$Stop > myLeftEdge & peaksIn$Start < myRightEdge)
		if ( length(toVisit)) {
			deltaCenter <- abs( myCenter - peaksIn$Center[toVisit])
			smlOrder <- order( deltaCenter)
			toVisit <- toVisit[ smlOrder]
			# see if we need to extend our footprint to contain other peaks
			for ( k in toVisit) {
				if ( peaksIn$Center[k] > myLeftEdge & peaksIn$Center[k] < myRightEdge) {
					myLeftEdge <- min( myLeftEdge, peaksIn$Start[k])
					myRightEdge <- max( myRightEdge, peaksIn$Stop[k])
					isAvailable[k] <- FALSE
				}
			}
		}

		# this peak may have expanded some, save it, mark it, and continue
		peakNow$POSITION <- myLeftEdge
		peakNow$END <- myRightEdge
		peaksOut <- rbind( peaksOut, peakNow)
		if ( nrow(peaksOut) %% 100 == 0) cat( "\r" , i, peakNow$GENE_ID[1], nrow(peaksOut), sum( isAvailable))
	}

	# we have the set
	# re-format to be a gene map
	gmap <- getCurrentGeneMap()
	gid <- paste( peaksOut$GENE_ID, peaksOut$Center, sep="::")
	wh <- match( peaksOut$GENE_ID, gmap$GENE_ID)
	sid <- gmap$SEQ_ID[ wh]
	prods <- gmap$PRODUCT[ wh]
	strand <- ifelse( peaksOut$Strand == "Plus", "+", "-")
	nExons <- peaksOut$END - peaksOut$POSITION + 1
	who <- peaksOut$SampleID

	out <- data.frame( "GENE_ID"=gid, "POSITION"=peaksOut$POSITION, "END"=peaksOut$END, "SEQ_ID"=sid,
				"STRAND"=strand, "NAME"=peaksOut$GENE_ID, "PRODUCT"=prods, "N_EXON_BASES"=nExons,
				"SampleName"=who, stringsAsFactors=F)

	# let's carry forward some RIP peak features people may want
	otherColumns <- c( "VPM", "P_Value", "Score", "Height", "Center", "Width", 
				"Type", "Start", "Stop", "Model_Error")
	out2 <- peaksOut[ , otherColumns]
	out <- cbind( out, out2, stringsAsFactors=F)

	ord <- order( out$POSITION)
	out <- out[ ord, ]
	if ( nrow(out)) rownames(out) <- 1:nrow(out)
	return( out)
}

