# WB_Tools.R

# assorted WiggleBin functions and constants

WB_PLUS_UNIQUE <- 1
WB_MINUS_UNIQUE <- 2
WB_PLUS_MULTI <- 3
WB_MINUS_MULTI <- 4
WB_UNIQUE_CNT <- 1
WB_UNIQUE_PCT <- 2



# set up default constants
`WB.defaults` <- function() {
	assign( "BinSize", 100, envir=WB_Env)
	assign( "BinStrands", 4, envir=WB_Env)		# + and - for unique and Multi
	assign( "SelfDetectability", NULL, envir=WB_Env)
	assign( "OtherDetectability", NULL, envir=WB_Env)
}


# get or set various WiggleBin items...
`WB_getBinSize` <- function() { return( get( "BinSize", envir=WB_Env)) }
`WB_setBinSize` <- function(x) { return( assign( "BinSize", value=x, envir=WB_Env)) }
`WB_setBinSizeBySpecies` <- function(x) { 
	
	mysize <- 100
	if ( x %in% MAMMAL_SPECIES) mysize <- 250
	assign( "BinSize", value=mysize, envir=WB_Env)
	return( mysize) 
}
`WB_getBinStrands` <- function() { return( get( "BinStrands", envir=WB_Env)) }


# get or set various Detectability objects...
`WB_getCurrentSelfDetectability` <- function() { return( get( "SelfDetectability", envir=WB_Env)) }
`WB_getCurrentOtherDetectability` <- function() { return( get( "OtherDetectability", envir=WB_Env)) }


# convert between nucleotide base position and wiggle bin number
`WB_bin2base` <- function( x, edge="left") { 

	binsize <- WB_getBinSize()

	switch( edge,
		"left" = return( (x - 1) * binsize + 1),
		"right" = return( x * binsize),
		"center" = return( x * binsize - round( binsize/2))
	)
	return( NA)
}

`WB_base2bin` <- function( x) { return( ceiling( x / WB_getBinSize())) }


# convert between gene/exon boundaries and bin range, using sensible edge clipping to try to
# minimize the bin induced gene edge problems
`WB_baseRange2binRange` <- function( firstBase, lastBase) {

	# get the first bin that is at least half inside the gene/exon
	bin1 <- rawbin1 <- WB_base2bin( firstBase)
	midBase <- WB_bin2base( bin1, edge="center")
	if ( midBase < firstBase) bin1 <- bin1 + 1

	# get the last bin that is at least half inside the gene/exon
	bin2 <- rawbin2 <- WB_base2bin( lastBase)
	midBase <- WB_bin2base( bin2, edge="center")
	if ( midBase > lastBase) bin2 <- bin2 - 1
	if ( bin1 > bin2) {
		bin1 <- rawbin1
		bin2 <- rawbin2
	}
	return( bin1 : bin2)
}


`WB_multiBaseRange2binRange` <- function( firstBases, lastBases) {

	out <- vector()
	for( i in 1:length( firstBases)) out <- base::append( out, WB_baseRange2binRange( firstBases[i], lastBases[i]))
	return( base::sort( unique.default( out)))
}


# get the stored read counts in the entire data structure
`WB_getTotalRawReads` <- function( WB) {

	out <- list( "Unique"=WB$Info$RawUniqueReads, "Multi"=WB$Info$RawMultiReads)
	return( out)
}


# add up all the reads in the entire data structure
`WB_sumAllBins` <- function( WB) {

	# accumulate each strand separately
	strandTotals <- rep( 0, times=WB_getBinStrands())
	for ( i in 1:length( WB$BinData)) strandTotals <- strandTotals + apply( WB$BinData[[i]], MARGIN=2, FUN=sum)

	out <- list( "Unique"=(strandTotals[WB_PLUS_UNIQUE] + strandTotals[WB_MINUS_UNIQUE]), 
		 "Multi"=(strandTotals[WB_PLUS_MULTI] + strandTotals[WB_MINUS_MULTI])) 
	return( out)
}


# add up all the bin counts in a specific bin set range
`WB_sumBinSetRange` <- function( binSet, bin1, bin2) {

	# accumulate each strand separately
	strandTotals <- rep( 0, times=WB_getBinStrands())

	# this range may be a single row, or a subset matrix
	if ( bin1 == bin2) {
		strandTotals <- binSet[ bin1, ]
	} else if ( bin2 > bin1) {
		strandTotals <- apply( binSet[ bin1:bin2, ], MARGIN=2, FUN=sum)
	}

	out <- list( "PlusUnique"=strandTotals[WB_PLUS_UNIQUE], 
		 "MinusUnique"=strandTotals[WB_MINUS_UNIQUE],
		 "PlusMulti"=strandTotals[WB_PLUS_MULTI],
		 "MinusMulti"=strandTotals[WB_MINUS_MULTI])
	return( out)
}


# find a BinSet by its SeqID name
`WB_getBinSetPtrFromSeqID` <- function( WB, seqid) { return( base::match( seqid, names( WB$BinData), nomatch=0)) }



# grab the detectability for a species and its 'other' if present
`WB_setSpeciesDetectability` <- function( ) {

	# make sure we have the right annotation objects in place
	curDetectability <- WB_getCurrentSelfDetectability()
	curOther <- WB_getCurrentOtherDetectability()

	curSpecies <- selfID <- getCurrentSpecies()
	otherID <- NULL
	otherSpecies <- setdiff( getCurrentTargetSpecies(), selfID)
	if ( length( otherSpecies) == 1) otherID <- otherSpecies[1]

	# allow the default 'other species' to be human...
	if ( is.null(otherID)) {
		# try to be a bit smarter...
		if ( selfID %in% c( "PbANKA", "PyAABL", "Py17X", "PyYM", "PchAS")) {
			otherID <- "Mmu_grc"
		} else if ( selfID %in% c( "MT_H37", "BCG", "Msmeg_mc2_155", "Mabsc", "Mchel", "MT_HN878")) {
			otherID <- "Mmu_grc"
		} else if (selfID %in% c( "PkH", "PCO", "PcoAH", "Pcy")) {
			otherID <- "MacMu"
		} else if (selfID %in% c( "Pf3D7", "PvSal1", "PvP01")) {
			otherID <- "Hs_grc"
		} else if (selfID %in% c( "Mmu_grc")) {
			otherID <- "Py17X"
		} else if (selfID %in% c( "MacMu", "MacFas")) {
			otherID <- "PcoAH"
		} else if (selfID %in% MAMMAL_SPECIES) {
			otherID <- "Pf3D7"
		} else if (selfID %in% PARASITE_SPECIES) {
			otherID <- "Hs_grc"
		} else {
			otherID <- "Hs_grc"
		}
	}

	needLoad <- ( is.null( curDetectability) || ( curDetectability$Species != selfID))
	 
	# also see if the 'other' species forces a need to re-load
	if ( ! is.null( otherID)) {
		if ( is.null( curOther) || (curOther$Species != otherID))  needLoad <- TRUE
	}

	if ( needLoad ) {
		cat( "..loading detectability" , selfID, otherID,"...")
		ans <- WB_loadSpeciesDetectability( selfID, otherID)
		cat( "   Done.\n")
	}

	if ( curSpecies != getCurrentSpecies()) setCurrentSpecies( curSpecies)

	return()
}


`WB_loadSpeciesDetectability` <- function( speciesID, otherID=NULL) {

	if ( getCurrentSpecies() != speciesID) setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()
	file <- paste( prefix, "selfUnique", "WB", sep=".")

	# load that 'selfUnique' detectability object,  this object is called 'WB'
	WB <- NA
	try( data( list=list(file), package="DuffyNGS", envir=environment()), silent=TRUE)
	WBself <- WB
	selfAnswer <- list( "Species"=speciesID, "Unique"=WBself)

	# load that 'otherDetectable' detectability object,  this object is called 'WB'
	otherAnswer <- NULL
	WB <- NA
	if ( ! is.null( otherID)) { 
		setCurrentSpecies( otherID)
		otherPrefix <- getCurrentSpeciesFilePrefix()
		otherFile <- paste( prefix, "by", otherPrefix,"detectable", "WB", sep=".")
		try( data( list=list(otherFile), package="DuffyNGS", envir=environment()), silent=TRUE)
	}
	WBother <- WB
	# if we got a good WB data set, then pass it back
	if ( ! any(is.na( WBother))) otherAnswer <- list( "Species"=otherID, "Detectable"=WBother)

	#out <- list( "Self"=selfAnswer, "Other"=otherAnswer)
	
	# set this as the current detectability
	assign( "SelfDetectability", value=selfAnswer, envir=WB_Env)
	assign( "OtherDetectability", value=otherAnswer, envir=WB_Env)

	#return( out)
	return( )
}
