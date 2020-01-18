# compressResults.R -- compress the entire 'Results' folder tree into a smaller more managable format
#			for copying, saving, transferring, etc.

`compressResults` <- function( path=NULL, overwrite=FALSE) {

	if ( is.null( path)) path <- getOptionValue( "Options.txt", "results.path", verbose=F, notfound="results")

	# intended to be run from the main Experiment level folder, with the intent of combining
	# subfolders into one 'CompressedResults.tar.gz' file for each results subfolder
	if ( ! file.exists( path)) stop( paste( "Results folder not found: ", path))

	# Step 1:  remove any clean-able files first
	cat( "\nBAM file cleanup:\n")
	cleanupBAMfiles( path)

	cat( "\nVelvet file cleanup:\n")
	cleanupVelvetFiles( path=file.path( path, "VelvetContigs"))

	cat( "\nHLA Typing file cleanup:\n")
	cleanupHLAtypingFiles( path=file.path( path, "HLA.typing"))

	# Step 2:  get the list of directories in this results folder
	dirSet <- dir( path, include.dirs=T, full.names=T)
	isDir <- file.info( dirSet)$isdir
	dirSet <- dirSet[ isDir]
	if ( !length( dirSet)) return(NULL)

	# Step 3:  visit each one, and act various ways depending on the contents
	out <- data.frame()
	for ( d in dirSet) {
		baseDir <- basename( d)

		# some are just giant compressed data that can't be further compressed
		if ( baseDir %in% c( "align", "fastq", "Viral.BAMS")) {
			cat( "\nIgnoring already compressed folder: ", baseDir)
			next
		}
		# some are folders of folders, where each subfolder should get compressed
		if ( baseDir %in% c( "DESeq", "EdgeR", "MetaResults", "RankProduct", "RoundRobin", "SAM", 
					"ConsensusProteins", "SieveAnalysis", "VariantCalls", "VelvetContigs")) {
			cat( "\nCompressing folders down inside: ", baseDir)
			ans <- compressFolderSubfolders( d, overwrite=overwrite)
			if ( ! is.null(ans)) out <- rbind( out, ans)
			next
		}
		# there are a few other 'folder of folders' that often have versions of names
		if ( grepl( "ChIPpeak|RIPpeak", baseDir )) {
			cat( "\nCompressing folders down inside: ", baseDir)
			ans <- compressFolderSubfolders( d, overwrite=overwrite)
			if ( ! is.null(ans)) out <- rbind( out, ans)
			next
		}
		# default is to compress the full folder
		cat( "\nCompressing folder: ", baseDir)
		ans <- compressFolder(d, overwrite=overwrite)
		if ( ! is.null(ans)) out <- rbind( out, ans)
	}
	
	# summarice the answer
	out <- apply( as.matrix( out), MARGIN=2, sum)
	out <- round( out)
	cat( "\nDone Compressing Folders: \n")
	out
}


`compressedFileName` <- function( path) {

	# given a full path to a folder, turn it into a file name that will hold the compressed version thereof...
	path <- gsub( "/+$", "", path)
	mydir <- dirname( path)
	myfile <- basename(path)
	outfile <- paste( myfile, "CompressedResults.tar.gz", sep=".")
	out.path <- file.path( mydir, outfile)
	out.path
}


`compressFolder` <- function( path, overwrite=FALSE, verbose=FALSE) {

	# given the full path name of a folder, compress the entire thing into a single file
	# then delete the original

	# first make sure its the folder spec itself
	path <- gsub( "/+$", "", path)
	compFile <- compressedFileName( path)

	# see if we can/must overwrite
	if ( file.exists( compFile)) {
		if ( ! overwrite) {
			cat( "\nCompressed File already exists: ", compFile, "  Can't overwrite..")
			return( NULL)
		}
		file.delete( compFile)
	}

	# make sure there is something to compress
	fset <- dir( path, recursive=T, full.names=T)
	if ( ! length( fset)) return( NULL)

	# create the compress command
	# paths 'could' have embadded blanks...
	if (verbose) {
		cmdline <- paste( "tar -czv -f '", compFile, "'  '", path, "'", sep="")
	} else {
		cmdline <- paste( "tar -cz -f '", compFile, "'  '", path, "'", sep="")
	}

	# do it
	system( cmdline, wait=TRUE)
	Sys.sleep( 0.1)

	# verify we see the new file
	if ( ! file.exists( compFile)) {
		cat( "\nCompression error. Result file not detected: ", compFile)
		return(NULL)
	}

	# let's count/measure before and after
	mbBefore <- sum( file.info( fset)$size, na.rm=T) / 1000000
	mbAfter <- file.info( compFile)$size / 1000000
	savings <- round( mbBefore - mbAfter, digits=3)

	# OK, clear to delete that folder
	unlink( path, recursive=TRUE)
	
	return( data.frame( "N_Folders"=1, "N_Files"=length(fset), "MB_Savings"=savings))
}


`compressFolderSubfolders` <- function( path, overwrite=FALSE, verbose=TRUE) {

	# given the full path name of a folder, only compress the folders under this folder

	# first make sure its the folder spec itself
	path <- gsub( "/+$", "", path)

	# get the set of subfolders to compress
	dirSet <- dir( path, recursive=F, full.names=T, include.dirs=T)
	isDir <- file.info( dirSet)$isdir
	dirSet <- dirSet[ isDir]
	if ( ! length( dirSet)) return(NULL)

	out <- data.frame()
	for ( d in dirSet) {
		baseDir <- basename( d)
		if (verbose) cat( "\nCompressing folder: ", baseDir)
		ans <- compressFolder(d, overwrite=overwrite, verbose=F)
		if ( ! is.null(ans)) out <- rbind( out, ans)
	}
	return( out)
}

