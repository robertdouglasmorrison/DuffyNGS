# fastCS2DNA.R

fastCS2DNA <- function( txt) {

	nReads <- length(txt)
	lengths <- base::nchar( txt)

	# make a copy to be the result storage
	txtOut <- txt

	ans <- .C( "fastColorSpace2DNA", 
			nReads=as.integer(nReads), 
			lengths=as.integer(lengths), 
			txtIn=txt,
			txtOut=txtOut)
	return( ans$txtOut)
}

