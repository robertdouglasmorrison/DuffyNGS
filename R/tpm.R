# tpm.R

# implement the TPM method of gene expression, as from Wagner, Theory of Bioscience, 2012


tpm <- function( readCount, geneLen, readLen=100) {

	# first normalize the read counts by gene length
	perK <- readLen / geneLen
	RperG <- readCount * perK

	# the sum of all genes is the T term in the denominator
	T <- sum(RperG)
	# tiny chance of having no aligned reads at all...
	if ( T < 1) T <- 1

	# make the TPM units for each gene
	tpm <- (readCount * readLen * 1000000) / (geneLen * T)

	# trim to something sensible?...
	return( tpm)
}

