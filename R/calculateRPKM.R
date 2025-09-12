`calculateRPKM` <-
function( readCount, exonBases=1000, totalReads=1000000) {

	# turn a gene's total read count into a "Reads Per Thousand_Exon_Bases per Million_Reads"
	
	# allow the tool to do either a single gene or a vector of genes
	if ( length(readCount) < 2) {
		thousandBasesFac <-  max( exonBases, 1) / 1000
	} else {
		thousandBasesFac <-  pmax( exonBases, 1) / 1000
	}
	millionReadsFac <- max( totalReads, 1) / 1000000

	rpkm <- readCount / thousandBasesFac / millionReadsFac
	return( rpkm)
}

