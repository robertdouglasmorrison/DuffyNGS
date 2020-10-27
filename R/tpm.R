# tpm.R

# implement the TPM method of gene expression, as from Wagner, Theory of Bioscience, 2012


tpm <- function( readCount, geneLen, readLen=100, minReadsPerSpecies=10000) {

	# the sum of all genes is the T term in the denominator

	# there is a chance of having very small number of aligned reads, which breaks TPM
	readCountForT <- readCount
	totalReads <- sum( readCount)
	if (totalReads < minReadsPerSpecies) {
		if ( all( readCountForT == 0)) readCountForT <- rep.int(1,length(readCount))
		readCountForT <- readCountForT * minReadsPerSpecies / sum( readCountForT)
	}
	# first normalize the read counts by gene length
	perK <- readLen / geneLen
	RperG <- readCountForT * perK
	T <- sum(RperG)


	# make the TPM units for each gene
	tpm <- (readCount * readLen * 1000000) / (geneLen * T)

	# trim to something sensible?...
	return( tpm)
}

