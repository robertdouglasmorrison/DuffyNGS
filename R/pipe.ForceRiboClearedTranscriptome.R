# pipe.ForceRiboClearedTranscriptome.R

# force a transcriptome to have zeroed expression of RiboCleared genes.
# since the ribo clearing is an alignment phase task, and each alignment program has limitations, allow
# an explicit zeroing of genes slated for removal by ribo clearing.  Then force the recalculation of all
# expression units metrics to adjust for the zeroing of ribo cleared genes.

`pipe.ForceRiboClearedTranscriptome` <- function( sampleID, annotationFile="Annotation.txt", optionsFile="Options.txt", 
				speciesID=getCurrentSpecies(), results.path=NULL, verbose=TRUE) {

	if ( speciesID != getCurrentSpecies()) setCurrentSpecies( speciesID)
	speciesPrefix <- getCurrentSpeciesFilePrefix()

	if ( is.null( results.path)) {
		optT <- readOptionsTable( optionsFile)
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}

	# get existing transcriptome
	trans.path <- file.path( results.path, "transcript")
	file <- file.path( trans.path, paste( sampleIDs, speciesPrefix, "Transcript.txt", sep="."))

	if ( ! file.exists( file)) {
		cat( "\nError: required transciptome file not found:  ", file)
		return( NULL)
	}

	rrnaMap <- subset( getCurrentRrnaMap(), CLEAR == TRUE)
	
	# read the original transcriptome 
	tbl <- read.delim( file, as.is=T)
	
	# find the genes that need zeroing out
	rowsToZero <- match( rrnaMap$GENE_ID, tbl$GENE_ID, nomatch=0)
	rowsToZero <- rowsToZero[ rowsToZero > 0]
	nToClear <- length(rowsToZero)
	if ( ! nToClear) {
		cat( "\nWarn: no RiboClear genes found in transcriptome")
		return( NULL)
	}
	
	if (verbose) cat( "\nZeroing ", nToClear, "genes flagged for RiboClearing")
	# we will clip the uncertainies at their lower bound
	minSigma <- min( tbl$SIGMA_M, na.rm=T)
	
	# now visit every expression units column, to zero that expression value
	ExpressColumns <- c("RPKM_M", "RPKM_U", "READS_M", "READS_U", "TPM_M", "TPM_U", "SIGMA_M", "SIGMA_U")
	
	newtbl <- tbl
	for ( j in ExpressColumns) {
		if ( ! (j %in% colnames(tbl))) next
		v <- tbl[[j]]
		vnew <- v
		vnew[ rowsToZer0] <- 0
		newtbl[[j]] <- vnew
	}
	newtbl$SIGMA_M <- ifelse( newtbl$SIGMA_M < minSigma, minSigma, newtbl$SIGMA_M)
	newtbl$SIGMA_U <- ifelse( newtbl$SIGMA_U < minSigma, minSigma, newtbl$SIGMA_U)

	# now re-force all the expression units to recalculate, to adjust for the zeroing of the ribo cleared genes
	nReadsM <- sum( newtbl$READS_M, na.rm=T)
	nReadsU <- sum( newtbl$READS_U, na.rm=T)
	newtbl$RPKM_M <- calculateRPKM( newtbl$READS_M, newtbl$N_EXON_BASES, totalReads=nReadsM)
	newtbl$RPKM_U <- calculateRPKM( newtbl$READS_U, newtbl$N_EXON_BASES, totalReads=nReadsU)
	newtbl$TPM_M <- tpm( newtbl$READS_M, newtbl$N_EXON_BASES)
	newtbl$TPM_U <- tpm( newtbl$READS_U, newtbl$N_EXON_BASES)
	if ( all( c("TPM_M","RANK_M") %in% colnames(newtbl))) {
        newtbl$RANK_M <- rank( newtbl$TPM_M, ties.method="min")
        newtbl$RANK_U <- rank( newtbl$TPM_U, ties.method="min")
	}
		
	# put the new transcripome back into expression order
	exOrd <- order( newtbl$RPKM_M, decreasing=T)
	newtbl <- newtbl[ exOrd, ]
	#W overwrite the existing transcriptome file
	write.table( newtbl, file, sep="\t", quote=F, row.names=F)
	if (verbose) cat( "\nDone.\n")

	# pass back the new transcriptome, invisibly
	out <- newtbl
	return( invisible(out))
}
