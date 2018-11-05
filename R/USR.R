# USR.R
#
#	USR	- Unique Short Reads
#
#	The universe of all unique short reads, to be searched for Consensus Reads.
#



`USR_setup` <- function( file, sampleID, resultsPath=".", trim5=0, trim3=0, Nkeep=NULL, 
			reload=FALSE, dropPolyN=TRUE) {

	filePath <- file.path( resultsPath, "USR")
	if ( ! file.exists( filePath)) dir.create( filePath, showWarnings=FALSE)
	USRfile <- file.path( filePath, paste( sampleID, "USR.rda", sep="."))

	if ( (!file.exists( USRfile)) || reload) {
		# create the universe of Unique Short Reads
		cat( "\n\nCreating USR dataset for:  ", sampleID)
	
		newUSR( file, trim5=trim5, trim3=trim3, Nkeep=Nkeep, dropPolyN=dropPolyN)

		# when we don't read up the whole file, some counts are off
		if ( ! is.null( Nkeep)) {
			trueTotalReads <- getFileLineCount( file, sampleID) / 4
			totalSeen <- USR_GrandTotal + USR_nDropPolyN
			USR_pctCoverage <<- totalSeen / trueTotalReads
		}

		saveUSRcontext( USRfile)
	} else {
		cat( "\n\nLoading existing USR dataset for:  ", sampleID)
		load( USRfile, envir=.GlobalEnv)
		cat( "   N_USR: ", USR_nTotal)
	}

	return( list( "nUSR"=USR_nTotal, "nReads"=USR_GrandTotal, "nPolyN"=USR_nDropPolyN, "USR_File"=USRfile))
}


`USR_cleanup` <- function() {

	# clean up
	if ( ! exists( "USR_seq")) return()
	rm( USR_seq, USR_nTotal, USR_counts, USR_GrandTotal, envir=.GlobalEnv)
	rm( list=c("USR_nDropPolyN", "USR_pctCoverage", "USR_AdapterHits"),  envir=.GlobalEnv)
	if ( exists("USR_colors")) rm( USR_colors, envir=.GlobalEnv)
}


newUSR <- function( file, trim5=0, trim3=0, Nkeep=NULL, dropPolyN=TRUE) {

	ans <- buildUSRfromFile( file, trim5, trim3, Nkeep=Nkeep, dropPolyN=dropPolyN)
	if ( length(ans$seq) < 10) stop( "newUSR:  not enough short reads to proceed.")
	
	# for every unique read, get its repeat count and factoring of its bases
	seq <- toupper( ans$seq)
	cnts <- ans$counts

	# let's not do the colors till we need them...
	#baseColors <- lapply( seq, FUN=myBaseFactorsFunc)
	# order from 'most frequent down'
	ord <- base::order( cnts, decreasing=TRUE)

	# save to global storage for speed
	USR_nTotal <<- ans$N
	USR_seq <<- seq[ord]
	USR_counts <<- cnts[ord]
	USR_GrandTotal <<- sum( USR_counts)
	USR_nDropPolyN <<- ans$nPolyN
	USR_pctCoverage <<- 1.0
	USR_AdapterHits <<- NULL
	return()
}


getUSRentry <- function( id) {
	if ( id < 1 | id > USR_nTotal) stop( paste( "getUSRentry:  invalid USR id numer: ", id))
	return( list( "seq"=USR_seq[id], "count"=USR_counts[id], "colors"=USR_BaseFactorsFunc(USR_seq[id])))
}


writeUSRasFasta <- function( fileout, N=NULL) {

	if ( is.null(N)) N <- length( USR_counts)
	cat( "\nWriting ", N, "Unique Short Reads to file: ", fileout)
	nams <- base::paste( "USR_", 1:N, "::", USR_counts[1:N], sep="")
	txt <- base::paste( ">", nams, "\n", USR_seq[1:N], sep="")
	writeLines( txt, con=fileout)
	cat( "  Done.\n")
	return()
}


buildUSRfromFile <- function( file, trim5=0, trim3=0, Nkeep=NULL, dropPolyN=TRUE) {
	
	#  force evaluation of arguments now...
	trim5 <- trim5
	trim3 <- trim3

	file <- allowCompressedFileName(file)
	cat("\nReading file: ", file, "\n")
	conIn <- openCompressedFile( file, open="r")
	chunkSize <- 1000000
	more <- TRUE
	allTables <- vector( mode="list")
	nTables <- 0
	nPolyN <- 0

	cat("\n")

	while (more) {
		raw <- readLines( conIn, n=chunkSize)
		NR <- length(raw)
		if ( NR < chunkSize) more <- FALSE
		if ( NR < 2) break

		# get just the sequence line
		lines <- vector()
		isFastq <- ( base::substr( raw[1], 1,1) == "@")
		isFasta <- ( base::substr( raw[1], 1,1) == ">")
		if ( isFastq) lines <- seq( 2, NR, by=4)
		if ( isFasta) lines <- seq( 2, NR, by=2)
		if ( length( lines) < 1) stop( "extractUSRfromFile: unknown file type.")
	
		thisSet <- raw[ lines]
	
		if ( trim5 > 0 || trim3 > 0) {
			cat( "  trim5=", trim5, " trim3=", trim3)
			from <- 1
			to <- max( base::nchar( thisSet))
			from <- from + trim5
			to <- to - trim3
			thisSet <- base::substr( thisSet, from, to)
		}

		allSeqs <- base::table( thisSet)
		if ( dropPolyN) {
			nNs <- sapply( names(allSeqs), FUN=function(x) {
				 	chr <- strsplit( x, split="")[[1]]
					return( ( sum( chr == "N") + sum( chr == ".")))
				})
			chrLen <- base::nchar( names(allSeqs)[1])
			drops <- which( nNs > (chrLen * 0.2))
			if ( length(drops) > 0) {
				allSeqs <- allSeqs[ -drops]
				class( allSeqs) <- "table"
				cat( "  PolyN: ", length(drops), "  ")
			}
			nPolyN <- nPolyN + length(drops)
		}

		if ( length(allSeqs) < 1) {
			cat( "  No reads left in chunk...\n")
			next
		}

		nTables <- nTables + 1
		allTables[[ nTables]] <- allSeqs

		thisTotal <- sum( sapply( allTables, length))
		curTotal <- sum( sapply( allTables, sum))

		cat( "  Chunk: ",nTables, "  Unique:", formatC(thisTotal,format="d",big.mark=","),
			"  Total: ", formatC(curTotal,format="d",big.mark=","), "\n")

		if ( !is.null(Nkeep)) {
			if ( curTotal >= (Nkeep * 5)) {
				cat("\nLimiting out...\n")
				break
			}
		}
	} # while...

	close( conIn)

	# now merge subsets...
	if ( length( allTables) < 2) {
		newAll <- vector( mode="list")
		newAll[[1]] <- allTables[[1]]
	}
	while ( length( allTables) > 1) {
		cat( "\nMerging ", length(allTables), "subsets...")
		newAll <- vector( mode="list")
		ones <- seq( 1, length(allTables), by=2)
		for ( i in 1:length(ones)) {
			thisOne <- ones[i]
			if ( thisOne < length(allTables)) {
				newOne <- mergeTables( allTables[[thisOne]], allTables[[thisOne+1]] )
			} else {
				newOne <- allTables[[thisOne]]
			}
			cat(".")

			# allow a trimming to keep the total size down..
			if ( ! is.null( Nkeep)) {
				trimSize <- round( Nkeep)
				if ( length( newOne) > trimSize) {
					cat( "trim")
					ord <- base::order( newOne, decreasing=TRUE)
					keep <- base::sort( ord[ 1:trimSize])
					newOne <- newOne[ keep]
				}
			}
			newAll[[i]] <- newOne
		}
		allTables <- newAll
		cat( "  N_Unique: ", formatC(sum(sapply(allTables,length)),format="d",big.mark=","),
			"  N_Total: ", formatC(sum(sapply(allTables,sum)),format="d",big.mark=","))
	}
	cat("\nDone.\n")
	
	newOne <- newAll[[ 1]]

	if ( ! is.null( Nkeep)) {
		trimSize <- round( Nkeep)
		if ( length( newOne) > trimSize) {
			cat( "\nFinal trimming..")
			ord <- base::order( newOne, decreasing=TRUE)
			keep <- base::sort( ord[ 1:trimSize])
			newOne <- newOne[ keep]
		}
	}

	return( list( "N"=length(newOne), "seq"=names( newOne), "counts"=newOne, "nPolyN"=nPolyN))
}


saveUSRcontext <- function( fileout) {
	# cat( "\nSaving USR context...\n")
	save( USR_nTotal, USR_GrandTotal, USR_seq, USR_counts, USR_nDropPolyN, USR_pctCoverage, 
			USR_AdapterHits, file=fileout)
	cat( "\nWrote USR file: ", fileout, "\n")
	return()
}


calc2USRalignment <- function( id1, id2) {

	GAP_OPEN_COST <- (-30)
	MISMATCH_COST <- (-12)
	return( pairwiseAlignment( USR_seq[id1], USR_seq[id2], type="overlap", gapOpening=GAP_OPEN_COST))
}


`USR_BaseFactorsFunc` <- function(bases, needSplit=TRUE) {

	# turn a character string of DNA into a vector of base IDs
	if ( needSplit) bases <- strsplit(bases, split="")[[1]]

	# this may be faster...
	bases[ bases == "."] <- "N"
	return( base::match( bases, c( "A", "C", "G", "T", "N", "-"), nomatch=0))
}

