# cutadapt.R -- wrapper to call Marcel Martin's CUTADAPT base trimmer

`cutadapt` <- function( file1, file2=NULL, 
			adapt1=myReverseComplement("CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT"), 
			adapt2=myReverseComplement("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"), 
			path=".", min.length=NULL, cutadaptArgs="", cutadaptProgram="~/bin/cutadapt",
			cores=4) {

	# at least one file is required, use it to deduce suffix/compression,
	# and generate out output filename(s)
	if ( isFASTQ <- (length( grep( "\\.fastq", file1[1])) == 1)) {
		outfile1 <- sub( "\\.fastq", ".trimmed.fastq", file1[1])
		if ( ! is.null(file2)) outfile2 <- sub( "\\.fastq", ".trimmed.fastq", file2[1])
	} else if ( isFASTA <- (length( grep( "\\.fasta", file1[1])) == 1)) {
		outfile1 <- sub( "\\.fasta", ".trimmed.fasta", file1[1])
		if ( ! is.null(file2)) outfile2 <- sub( "\\.fasta", ".trimmed.fasta", file2[1])
	} else if ( isFQ <- (length( grep( "\\.fq", file1[1])) == 1)) {
		outfile1 <- sub( "\\.fq", ".trimmed.fq", file1[1])
		if ( ! is.null(file2)) outfile2 <- sub( "\\.fq", ".trimmed.fq", file2[1])
	} else if ( isFA <- (length( grep( "\\.fa", file1[1])) == 1)) {
		outfile1 <- sub( "\\.fa", ".trimmed.fa", file1[1])
		if ( ! is.null(file2)) outfile2 <- sub( "\\.fa", ".trimmed.fa", file2[1])
	} else {
		stop( paste( "Could not resolve FASTQ/FASTA file extension type: ", file1))
	}
	isPaired <- exists( "outfile2")

	# prepend the path if one was given
	if ( !is.null(path) && !is.na(path) && (path != "")) {
		file1 <- file.path( path, file1)
		outfile1 <- file.path( path, outfile1)
		if ( isPaired) {
			file2 <- file.path( path, file2)
			outfile2 <- file.path( path, outfile2)
		}
	}
	if ( ! file.exists( file1)) stop( paste( "Input FASTQ/FASTA file not found: ", file1))
	if (isPaired) {
		if ( ! file.exists( file2)) stop( paste( "Input mate2 FASTQ/FASTA file not found: ", file2))
	}

	# discard 'too small' reads after trimming?
	minLenTrimText <- ""
	if ( ! is.null( min.length)) {
		minLenTrimText <- paste( " --minimum-length=", as.integer(min.length), " ", sep="")
	}

	# allow multi-threading
	coresText <- if (is.null(cores)) "" else paste( "--cores", as.character(cores))

	# make that command line
	if ( isPaired) {
		cmdLine <- paste( cutadaptProgram, "-a", adapt1, "-A", adapt2, coreText,
				"-o", outfile1, "-p", outfile2, minLenTrimText, cutadaptArgs, file1, file2, sep=" ")
	} else {
		cmdLine <- paste( cutadaptProgram, "-a", adapt1, coreText, "-o", outfile1, minLenTrimText, cutadaptArgs, file1, sep=" ")
	}
	catch.system( cmdLine)

	out <- outfile1
	if (isPaired) out <- c( out, outfile2)
	return(out)
}
