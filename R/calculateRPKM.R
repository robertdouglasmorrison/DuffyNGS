`calculateRPKM` <-
function( readCount, exonBases=1000, totalReads=1000000) {

	# turn a gene's total read count into a "Reads Per Thousand_Exon_Bases per Million_Reads"
	thousandBasesFac <-  exonBases / 1000
	millionReadsFac <- totalReads / 1000000

	rpkm <- readCount / thousandBasesFac / millionReadsFac
	return( rpkm)
}

