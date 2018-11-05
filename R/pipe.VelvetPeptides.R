# pipe.VelvetPeptides.R


`pipe.VelvetPeptides` <- function( sampleID, annotationFile="Annotation.txt",
		optionsFile="Options.txt", velvetPath="~/Velvet/velvet_1.2.10",
		kmerSize=25, minLength=200, minCoverage=3, maxCoverage=3000,
		useNoHits=TRUE, fastq.keywords=NULL, buildHash=TRUE, buildContigs=TRUE,
		verbose=TRUE) {

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\n\nStarting 'VelvetPeptides' for Sample:     ", sampleID, "\n\n")
	}

	optT <- readOptionsTable( optionsFile)
	resultsPath <- getOptionValue( optT, "results.path", notfound=".")

	outpath <- file.path( resultsPath, "VelvetPeptides", sampleID)
	if ( ! file.exists( outpath)) dir.create( outpath, recursive=T)

	# gather all the fastq files
	keywords <- vector()
	filesArgs <- ""
	if ( ! is.null( fastq.keywords)) keywords <- c( keywords, fastq.keywords)
	if ( useNoHits) keywords <- unique( c( keywords, "noHits"))

	for ( key in keywords) {
		f <- paste( sampleID, key, "fastq", sep=".")
		f <- file.path( resultsPath, "fastq", f)
		f <- allowCompressedFileName(f)
		fmt <- if ( regexpr( "gz$", f) > 0) "-fastq.gz" else "-fastq"
		filesArgs <- paste( filesArgs, fmt, f, sep=" ")
	}
	if ( filesArgs == "") {
		cat( "\nNo Fastq files and/or keywords given...  Quitting.")
		return()
	}

	# step 1:  VelvelH to build the Velvet hash table
	if ( buildHash) {
		cmdline <- paste( file.path( velvetPath, "velveth"), outpath, kmerSize, filesArgs, sep=" ")
		cat( "\n\nCalling VelvetH to build kmer hash table...\n")
		system( cmdline)
	}

	# step 2:  VelvetG to build the contigs
	if (buildContigs) {
		cmdline <- paste( file.path( velvetPath, "velvetg"), outpath, " -exp_cov auto")
		cat( "\n\nCalling VelvetG to build contigs...\n")
		system( cmdline)
	}

	# final call is to extract just the 'good' contigs
	cmdline <- paste( file.path( velvetPath, "velvetg"), outpath, "-min_contig_lgth", minLength, 
			"-exp_cov auto", "-cov_cutoff", minCoverage) 
	cat( "\n\nCalling VelvetG to extract 'good' contigs...\n")
	system( cmdline)

	# grap the final set
	contigs <- loadFasta( file.path( outpath, "contigs.fa"))
	cat( "\nN_Contigs: ", length( contigs$desc))

	# convert to peptides
	peps <- DNAtoBestPeptide( contigs$seq, clipAtStop=FALSE)
	faOut <- as.Fasta( paste( sampleID, contigs$desc, sep="_"), peps)

	outfile <- paste( sampleID, keywords[1], "velvetPeptides.fasta", sep=".")
	outfile <- file.path( outpath, outfile)
	writeFasta( faOut, file=outfile, line.width=60) 
	cat( "\nWrote Peptides .FASTA file: ", outfile)
	cat( "\n\nFinished 'VelvetPeptides' for Sample:     ", sampleID, "\n\n")

	return()
}


