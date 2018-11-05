`rosettaPvalue` <- function( readsA, sigmaA=rosettaSigma(readsA), normFacA=1, 
				readsB, sigmaB=rosettaSigma(readsB), normFacB=1, 
				nBasesInGene=1000) {

	# based on the Rosetta error model by John Castle
	# turn the raw values into normalized, by length of gene and total reads
	geneFactor <- 1000 / nBasesInGene
	facA <- geneFactor * normFacA
	normCntA <- readsA * facA
	facB <- geneFactor * normFacB
	normCntB <- readsB * facB

	difInten <- normCntA - normCntB 
	sumSigma <- sqrt( ((sigmaA*facA)^2) + ((sigmaB*facB)^2))
	xdev <- difInten / sumSigma
	Pvalue <- erfc( abs( xdev) / sqrt(2))
	return( Pvalue)
}

