# readVCFfile.R -- read up a generic VCF file, as created by BCFTOOLS Call


readVCF <- function( file, verbose=TRUE) {

	# try to generically read a VCF file
	# It contains any number of comment records, starting with "##", 
	# and then a true header row starting with "#", then rows of variant calls

	# so read the same content 2 ways, one for the data, and once again to extract the column names
	tbl <- read.delim( file, header=F, comment.char="#", sep="\t")
	if (verbose) cat( "\nNumber of VCF rows: ", nrow(tbl))

	txt <- readLines( file)
	commentRows <- grep( "^##", txt)
	if ( length(commentRows)) txt <- txt[ -commentRows]
	headerTxt <- sub( "^#", "", txt[1])
	headerTerms <- strsplit( headerTxt, split="\t")[[1]]

	# clean and relabel as needed
	colnames(tbl) <- headerTerms

	out <- tbl
	return(out)
}

