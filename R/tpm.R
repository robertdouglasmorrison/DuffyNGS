# tpm.R

# implement the TPM method of gene expression, as from Wagner, Theory of Bioscience, 2012


tpm <- function( readCount, geneLen, readLen=100, minReadsPerSpecies=10000) {

	# the sum of all genes is the T term in the denominator

	# there is a chance of having very small number of aligned reads, which breaks TPM math
	readCountForT <- readCount
	totalReadsForT <- totalReads <- sum( readCount)
	if (totalReads < minReadsPerSpecies) {
		# scale everyone to look like we got the minimum number of reads
		# special catch if zero reads to all genes
		if ( all( readCountForT == 0)) {
			readCountForT <- rep.int(1,length(readCount))
			totalReadsForT <- sum( readCountForT)
			readCountForT <- readCountForT * minReadsPerSpecies / totalReadsForT
		} else {
			readCountForT <- readCountForT * minReadsPerSpecies / totalReadsForT
			readCount <- readCount * minReadsPerSpecies / totalReads
		}
	}

	# first normalize the read counts by gene length
	perK <- readLen / geneLen
	RperG <- readCountForT * perK
	T <- sum(RperG)

	# make the TPM units for each gene
	tpm <- (readCount * readLen * 1000000) / (geneLen * T)

	return( tpm)
}

