# newWB.R

# create a new empty WiggleBin data structure

`newWB` <- function ( speciesID=getCurrentSpecies(), dataType="RNA-seq") {

	# make an empty Wiggle data structure to hold RNA sequencing reads...
	if ( speciesID != getCurrentSpecies()) setCurrentSpecies( speciesID)
	WB_setBinSizeBySpecies( speciesID)

	seqMap <- getCurrentSeqMap()
	nSeq <- nrow( seqMap)

	# each chromosome has a set of bins that span its length
	binList <- vector( mode="list", length=nSeq)
	names( binList) <- seqMap$SEQ_ID

	# each bin holds "how many reads touched this region" for both + and - strands
	nBases <- seqMap$LENGTH
	binSize <- WB_getBinSize()
	nbins <- WB_base2bin( nBases)
	for ( i in 1:nSeq) binList[[i]] <- matrix( 0, nrow=nbins[i], ncol=WB_getBinStrands())

	# also make a list of 'growable' information
	infoList <- list( "Species"=speciesID, "BinSize"=binSize, "FileName"=vector(), 
			"RawUniqueReads"=0, "RawMultiReads"=0, "DataType"=dataType)

	return( list( "BinData"=binList, "Info"=infoList))
}

