# DuffyNGS_Envir.R

# set up and manage the 'global environment' for the DuffyNGS package

# one shared environment for all DuffyNGS objects
DuffyNGSEnv <- new.env( hash=TRUE, parent=emptyenv())

# global constants
SPLICEMAP_COLUMNS <- c( "GENE_ID", "SPLICE_ID", "SEQ_ID", "SIZES", "POSITIONS", "ENDS")
RIBOCLEAR_COLUMNS <- c( "GENE_ID", "GROUP", "SEQ_ID", "START", "SEQ_END", "FASTA_POSITION", "FASTA_END")

verboseOutputDivider <- "\n=====================================\n"


# initial Bowtie parameters environment
BowtieEnv <- new.env( hash=TRUE, parent=emptyenv())
Bowtie2Env <- new.env( hash=TRUE, parent=emptyenv())
BOWTIE_COLUMNS <- c( "READ_ID", "STRAND", "SEQ_ID", "START", "READ_SEQ", "SCORE", "N_MATCH", "MISMATCHES")
SAM_COLUMNS <- c( "READ_ID", "FLAG", "SEQ_ID", "START", "MAPQ", "CIGAR", "MATE_SEQ_ID", "MATE_START", "TLEN", "READ_SEQ", "SCORE")

BAMTAG_ALIGNSCORE <- "AS"
BAMTAG_READWEIGHT <- "WT"
BAMTAG_GENEID <- "GI"
BAMTAG_SPLICEID <- "SP"
BAMTAG_MISMATCHSTRING <- "MD"


# WiggleBin environment
WB_Env <- new.env( hash=TRUE, parent=emptyenv())
WIG_Env <- new.env( hash=TRUE, parent=emptyenv())


# .onLoad() function is called when the package get loaded at runtime
# any needed setup goes here...
`.onLoad` <- function( libname, pkgname) {
}

`.onUnload` <- function( libpath) {
}


# .onAttach() function is called when the package gets attached,
# i.e. at the time the user first has access to the package
`.onAttach`  <- function( libname, pkgname) {

	# wake-up message
	cat( "\nPackage: \t\t", pkgname, "\n")

	# save the library and package name...
	assign( "LibraryName", value=libname, envir=DuffyNGSEnv)
	assign( "PackageName", value=pkgname, envir=DuffyNGSEnv)

	# initialize...
	DuffyNGS.defaults()
}


`DuffyNGS.defaults` <- function() {

	# bowtie2Par.defaults()
	WB.defaults()
	WIG.defaults()
}

