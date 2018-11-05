# getBowtie2AlignMetrics.R  -- read and parse the readd counts after the Bowtie2 run

`getBowtie2AlignMetrics` <- function( metricsFile, pairedEnd=FALSE, bowtieTiming=NULL, verbose=TRUE) {

	die <- FALSE
	if ( ! file.exists( metricsFile)) {
		die <- TRUE
		ans <- vector()
	} else {
		ans <- readLines( metricsFile)
		if ( length( ans) < 1) die <- TRUE
	}
	if ( die) {
		cat( "\nBowtie2 Metrics file not found or empty:  ", metricsFile)
		out <- list( "RawReads"=0, "UniqueReads"=0, "MultiReads"=0, 
			"NoHitReads"=0, "Time"=NULL)
		return(out)
	}

	if (verbose) {
		cat( "\nBowtie2 done:    \t", date(), "\n\n")
		print( ans);
		cat( "\n")
	}

	nReadsIn <- nUnique <- nMulti <- nReadsNoHit <- 0

	who <- grep( "reads; of these:", ans, fixed=T)
	if ( length( who) > 0) {
		nReadsIn <- tryCatch( as.integer( sub( " reads.+", "", ans[who[1]])), 
					warning = function(e) { return( 0)} )
	}
	if ( pairedEnd) {
		# when paired end, the raw read count is pairs
		nReadsIn <- nReadsIn * 2

		nUnique1 <- nUnique2 <- nUnique3 <- 0
		who <- grep( "aligned concordantly exactly 1 time", ans, fixed=T)
		if ( length( who) > 0) {
			nUnique1 <- tryCatch( as.integer( sub( " \\(.+", "", ans[who[1]])), 
					warning = function(e) { return( 0)} )
		}
		who <- grep( "aligned discordantly 1 time", ans, fixed=T)
		if ( length( who) > 0) {
			nUnique2 <- tryCatch( as.integer( sub( " \\(.+", "", ans[who[1]])), 
					warning = function(e) { return( 0)} )
		}
		who <- grep( "aligned exactly 1 time", ans, fixed=T)
		if ( length( who) > 0) {
			nUnique3 <- tryCatch( as.integer( sub( " \\(.+", "", ans[who[1]])), 
					warning = function(e) { return( 0)} )
		}
		# case 1 and 2 are pairs, case 3 is singletons
		nUnique <- nUnique1*2 + nUnique2*2 + nUnique3

		who <- grep( "%) aligned 0 times", ans, fixed=T)
		if ( length( who) > 0) {
			nReadsNoHit <- tryCatch( as.integer( sub( " \\(.+", "", ans[who[1]])), 
					warning = function(e) { return( 0)} )
		}
		 
		nMulti1 <- nMulti2 <- 0
		who <- grep( "aligned concordantly >1 times", ans, fixed=T)
		if ( length( who) > 0) {
			nMulti1 <- tryCatch( as.integer( sub( " \\(.+", "", ans[who[1]])), 
					warning = function(e) { return( 0)} )
		}
		who <- grep( "aligned >1 times", ans, fixed=T)
		if ( length( who) > 0) {
			nMulti2 <- tryCatch( as.integer( sub( " \\(.+", "", ans[who[1]])), 
					warning = function(e) { return( 0)} )
		}
		# case 1 is pairs, case 2 is singletons
		nMulti <- nMulti1*2 + nMulti2
	} else {
		who <- grep( "aligned exactly 1 time", ans, fixed=T)
		if ( length( who) > 0) {
			nUnique <- tryCatch( as.integer( sub( " \\(.+", "", ans[who[1]])), 
					warning = function(e) { return( 0)} )
		}
		who <- grep( "aligned 0 times", ans, fixed=T)
		if ( length( who) > 0) {
			nReadsNoHit <- tryCatch( as.integer( sub( " \\(.+", "", ans[who[1]])), 
					warning = function(e) { return( 0)} )
		}
		who <- grep( "aligned >1 times", ans, fixed=T)
		if ( length( who) > 0) {
			nMulti <- tryCatch( as.integer( sub( " \\(.+", "", ans[who[1]])), 
					warning = function(e) { return( 0)} )
		}
	}

	timeAns <- NULL
	if ( ! is.null( bowtieTiming)) {
		timeAns <- elapsedProcTime( bowtieTiming$timeIn, bowtieTiming$timeOut, N=nReadsIn)
	}

	out <- list( "RawReads"=nReadsIn, "UniqueReads"=nUnique, "MultiReads"=nMulti, 
			"NoHitReads"=nReadsNoHit, "Time"=timeAns)

	# if we got mostly positive numbers, we can delete the file
	if ( sum( c( nReadsIn, nUnique, nMulti, nReadsNoHit) > 0) >= 2) file.delete( metricsFile)

	return( out);
}
