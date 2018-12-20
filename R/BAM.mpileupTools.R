# BAM.mpileupTools.R -- little routines to manipulate some of the SAMTOOLS mpileup data


GEN_CALL <- ","     # the SAMtools character for 'matches the genome'
INDEL_CHARS <- c( "+", "-")
DIGITS <- as.character( 0:9)
GSUB <- base::gsub
LAPPLY <- base::lapply
MATCH <- base::match
NCHAR <- base::nchar
ORDER <- base::order
PASTE <- base::paste
SAPPLY <- base::sapply
SORT <- base::sort
SUBSTR <- base::substr
TOUPPER <- base::toupper
UNION <- base::union
WHICH <- base::which
WHICH.MAX <- base::which.max


# translate the pileup string of base calls into one called base at each genomic location


# turn the cryptic MPILEUP base call text into the best base and a table of all calls

MPU.callBases <- function( baseCalls, referenceBase,  minDepth=3, debug=F,
				mode=c("counts","percents")) {
	
	refs <- referenceBase
	strs <- baseCalls
	N <- length(strs)
	mode <- match.arg( mode)

	# parse the SAMtools pileup string to call the actual bases

	# the quality call of the "begin read" flag might be a "$", so find the starts first!
	# drop the 'begin of read segment' flags
	if (debug) cat( "  ^beg..")
	strs <- GSUB( "\\^.", "", strs)

	# drop the 'end of read segment' flags
	if (debug) cat( "  $end..")
	strs <- GSUB( "$", "", strs, fixed=T)

	# turn all genome hits to forward strand
	if (debug) cat( "  strand..")
	strs <- GSUB( ".", GEN_CALL, strs, fixed=T)

	# turn all mismatches to forward strand
	strs <- TOUPPER( strs)

	# check for indels, that will be much slower/harder
	if (debug) cat( "  indels..")
	hasDelete <- WHICH(regexpr( "-", strs, fixed=T) > 0)
	hasInsert <- WHICH(regexpr( "+", strs, fixed=T) > 0)
	hasIndel <- SORT( UNION( hasDelete, hasInsert))

	# make the simple call based on what's most common
	if (debug) cat( "  split..")
	allbases <- strsplit( strs, split="")

	if (debug) cat( "  count..")
	myTABLE <- oneDtable

	baseCntsTable <- LAPPLY( allbases, FUN=myTABLE)
	if ( mode == "percents") {
		baseCntsTable <- LAPPLY( baseCntsTable, FUN=function(x) round( x * 100 / sum(x)))
	}
	bestbase <- SAPPLY( baseCntsTable, function(tbl) {
					if ( sum(tbl) < 1) return(GEN_CALL)
					names(tbl)[ WHICH.MAX(tbl)]
				})
	TotalBases <- length( bestbase)

	# the insertions are much harder to parse
	if ( length( hasIndel) > 0) {

		# first they get tabled differently
		if (debug) cat( "\nIndel Evaluations:", length( hasIndel), "\n")
		indelCntsTable <- LAPPLY( allbases[hasIndel], FUN=MPU.tablifyIndelCalls)
		bestIndelBase <- SAPPLY( indelCntsTable, function(tbl) {
						if ( sum(tbl) < minDepth) return(GEN_CALL)
						names(tbl)[ WHICH.MAX(tbl)]
					})
		baseCntsTable[ hasIndel] <- indelCntsTable
		if ( mode == "percents") {
			baseCntsTable[ hasIndel] <- LAPPLY( baseCntsTable[ hasIndel], FUN=function(x) round( x * 100 / sum(x)))
		}
		bestbase[hasIndel] <- bestIndelBase

		# next we need to step along and see how they impact their neighboring calls
		foundInserts <- vector()
		mustSkips <- vector()
		for (k in hasIndel) {
			# if already covered by a previous indel, ignore
			if ( k %in% mustSkips) next
			# if the top choice has no indel, ignore
			mybest <- bestbase[k]
			if ( regexpr( "-|\\+", mybest) < 1) next
			
			# which sign won?
			myFirstBase <- SUBSTR(mybest,1,1)
			mysign <- SUBSTR( mybest, 2,2)
			myNbase <- NCHAR(mybest) - 2
			indelLetters <- SUBSTR( mybest, 3, NCHAR(mybest))
			if ( myNbase < 1) next

			# an 'insert' means extra letters after me...
			if ( mysign == "+") {
				# since we're building the string, put in the'real base' now
				if ( myFirstBase == GEN_CALL) myFirstBase <- refs[k]
				myAns <- PASTE( myFirstBase, indelLetters, sep="")
				bestbase[k] <- myAns
				foundInserts <- c( foundInserts, mybest)
			} else {
				# deletion is much tougher, we need to make sure that 'most' of
				# the deleted bases are also called 'deleted' too...
				toCheck <- (k+1) : (k+myNbase)
				# don't check past the end of the data
				toCheck <- toCheck[ toCheck <= TotalBases]
				if (length( toCheck) > 0) {
					nDel <- sum( (bestbase[ toCheck] == "*"), na.rm=T)
					if (debug) cat( "\nChecking downstream bases: ", k, "|", myNbase,
							"|", toCheck, "|", nDel)
					if ( nDel > (myNbase/2)) {
						# YES, really delete...'after me'
						bestbase[k] <- myFirstBase
						bestbase[toCheck] <- "*"
						foundInserts <- c( foundInserts, mybest)
					} else {
						# no, ignore... not overwhelming enough
						if (debug) cat( "  not enough downsteam agree: ", nDel, myNbase)
						bestbase[k] <- myFirstBase
						# for now, just revert to genomic call... (not ideal?)
						bestbase[toCheck] <- refs[ toCheck]
					}
				}
				# these downstream ones have been done, so skip em
				mustSkips <- toCheck
			}
		}
		if (debug) {
			cat( "\nTable of found indels:\n")
			print( table( foundInserts))
		}
	}

	if (debug) cat( "  final calls..")
	# the genome is the default
	out <- refs

	# SNPs are just a substitution
	isSNP <- WHICH( bestbase %in% c( "A","C","G","T"))
	if ( length( isSNP) > 0) out[ isSNP] <- bestbase[ isSNP]

	# deletions
	isDEL <- WHICH( bestbase == "*")
	if ( length( isDEL) > 0) out[ isDEL] <- ""

	# insertions
	isIN <- WHICH( NCHAR( bestbase) > 1)
	if ( length( isIN) > 0) out[ isIN] <- bestbase[ isIN]

	return( list( "call"=out, "depth.table"=baseCntsTable))
}


MPU.strandDepth <- function( baseCalls, pos=NULL, debug=F) {
	
	strs <- baseCalls
	N <- length(strs)

	# parse the SAMtools pileup string to call the actual bases

	# drop the 'end of read segment' flags
	if (debug) cat( "  $end..")
	strs <- GSUB( "$", "", strs, fixed=T)

	# drop the 'begin of read segment' flags
	if (debug) cat( "  ^beg..")
	strs <- GSUB( "\\^.", "", strs)

	# check for indels, that will be much slower/harder
	if (debug) cat( "  indels..")
	hasDelete <- WHICH(regexpr( "-", strs, fixed=T) > 0)
	hasInsert <- WHICH(regexpr( "+", strs, fixed=T) > 0)
	hasIndel <- SORT( UNION( hasDelete, hasInsert))

	# make the simple call based on what's most common
	if (debug) cat( "  split..")
	allbases <- strsplit( strs, split="")

	if (debug) cat( "  count..")
	myTABLE <- oneDtable

	# comma and dot are not valid names, so alter
	FWD_MARKS <- c( ".", "A","C","G", "T", "N")
	REV_MARKS <- c( ",", "a","c","g", "t", "n")

	baseCntsTable <- LAPPLY( allbases, FUN=myTABLE)

	# the ones with indels need a more precise tool
	if ( length( hasIndel) > 0) {
		if (debug) cat( "  indel counting..")
		indelCntsTable <- LAPPLY( allbases[hasIndel], FUN=MPU.tablifyIndelCalls)
		baseCntsTable[ hasIndel] <- indelCntsTable
	}

	if (debug) cat( "  strand depths..")
	strandCounts <- SAPPLY( baseCntsTable, function(tbl) {
					fwdBins <- MATCH( FWD_MARKS, names(tbl), nomatch=0)
					fwdCnt <- sum( tbl[ fwdBins], na.rm=T)
					revBins <- MATCH( REV_MARKS, names(tbl), nomatch=0)
					revCnt <- sum( tbl[ revBins], na.rm=T)
					return( c( fwdCnt, revCnt))
				})

	out <- t( strandCounts)
	colnames( out) <- c( "FWD_DEPTH", "REV_DEPTH")
	if ( ! is.null( pos)) rownames(out) <- pos

	return( out)
}


# make a table from a single pileup that has indels -- needs special parsing...

MPU.tablifyIndelCalls <- function( callChars) {

	# given a vector of the 1-letter pileup calls that have indel strings, parse it into counts of each
	# accumulate tokens as we go
	calls <- vector()
	nout <- 0
	N <- length( callChars)

	# we must start at an actual call, never at an indel
	atCall <- TRUE
	checkIndel <- atNumber <- FALSE

	i <- 1
	while ( i <= N) {
		if (atNumber) {
			if (callChars[i] %in% DIGITS) {
				size <- (size * 10) + as.numeric( callChars[i])
				i <- i + 1
				next
			} else {
				atNumber <- FALSE
				# we know know how many characters the indel text is
				indelText <- callChars[ (i) : (i+size-1)]
				# overwrite the original call with this composite
				calls[ nout] <- PASTE( c( calls[nout], indelSign, indelText), collapse="")
				i <- i + size
				atChar <- TRUE
				next
			}
		}
		if (checkIndel) {
			if ( callChars[i] %in% INDEL_CHARS) {
				atNumber <- TRUE
				checkIndel <- FALSE
				indelSign <- callChars[i]
				size <- 0
				i <- i + 1
				next
			} else {
				checkIndel <- FALSE
				atCall <- TRUE
			}
		}
		if (atCall) {
			nout <- nout + 1
			calls[ nout] <- callChars[i]
			checkIndel <- TRUE
			i <- i + 1
			next
		}
		stop( "never should get here...")
	}
	return( oneDtable( calls))
}


# convert a table of base call counts into a character string

MPU.callTableToString <- function( callTables) {

	if ( is.list( callTables)) {
		#cat( "  toStrings..")
		ans <- SAPPLY( callTables, function(x) {
				vals <- as.character(x)
				nams <- as.character( names(x))
				ord <- ORDER( x, decreasing=T)
				PASTE( nams[ord], vals[ord], sep="=", collapse=";")
			})
	} else {
		x <- callTables
		vals <- as.character(x)
		nams <- as.character( names(x))
		ord <- ORDER( x, decreasing=T)
		ans <- PASTE( nams[ord], vals[ord], sep="=", collapse=";")
	}
	ans
}


# un-convert a table of base call counts from it's character string format

MPU.callStringToTable <- function( callStrings) {

	terms <- strsplit( callStrings, split=";")
	callTables <- lapply( terms, function(x) {
					if ( length(x) == 0) return( NULL)
					subterms <- unlist( strsplit( x, split="="))
					if ( length(subterms) == 0) return( NULL)
					namptrs <- seq.int( 1, length(subterms), by=2)
					nams <- subterms[ namptrs]
					vals <- as.numeric( subterms[ namptrs + 1])
					dim(vals) <- length(vals)
					dimnames(vals) <- list( nams)
					class(vals) <- "table"
					vals
				})
	callTables
}


MPU.callStringsToMatrix <- function( callStrings) {

	N <- length(callStrings)
	CALLS <- c( ",", "A", "C", "G", "T", "N", "Indel")
	out <- matrix( 0, nrow=N, ncol=length( CALLS))
	colnames(out) <- CALLS
	if ( N < 1) return( out)
	
	listOfTables <- MPU.callStringToTable( callStrings)

	lapply( 1:N, function(x) {
			thisTable <- listOfTables[[x]]
			if ( is.null( thisTable)) return()
			where <- MATCH( names( thisTable), CALLS, nomatch=0)
			out[ x, where] <<- thisTable[ where > 0]
			# any that don't match CALLS are indels
			if ( any( where == 0)) {
				nIndels <- sum( thisTable[ where == 0], na.rm=T)
				out[ x, length(CALLS)] <<- nIndels
			}
			return()
		})
	return(out)
}	
	

`MPU.callMatrixToBaseCounts` <- function( m, referenceBase, indelsToo=FALSE, normalize=FALSE) {

	# given a matrix of base counts with column names ", A C G T N Indel", turn it to just A,C,G,T

	# step 1: extract the main part for the answer
	if (indelsToo) {
		out <- m[ , c(2:5,7), drop=FALSE]
		ignore <- m[ ,6, drop=FALSE]
		Nout <- 5
	} else {
		out <- m[ , 2:5, drop=FALSE]
		ignore <- m[ ,6:7, drop=FALSE]
		Nout <- 4
	}

	# step 2: assign the 'genomic' to the right column
	isA <- WHICH( referenceBase == "A")
	if ( length(isA)) out[ isA, 1] <- m[ isA, 1]
	isC <- WHICH( referenceBase == "C")
	if ( length(isC)) out[ isC, 2] <- m[ isC, 1]
	isG <- WHICH( referenceBase == "G")
	if ( length(isG)) out[ isG, 3] <- m[ isG, 1]
	isT <- WHICH( referenceBase == "T")
	if ( length(isT)) out[ isT, 4] <- m[ isT, 1]

	# step 3:  normalize?
	if ( normalize) {
		rowSums <- apply( out, MARGIN=1, sum, na.rm=T)
		others <- apply( ignore, MARGIN=1, sum, na.rm=T)
		ncalls <- rowSums + others
		for (j in 1:Nout) out[ , j] <- out[ ,j] * 100 / ncalls
		out <- round( out, digits=2)
	}

	out
}


# calculate a P-value from the base call counts, given the total counts from the entire 
# population (chromosome)

callTable.Pvalue <- function( callTable, callSumsTable) {

	# given a table of call frequencies for one base position,
	# see how unlikely given the population frequencies

	# to make this fair to the ChiSq test, let's scale the counts down to a
	# sane maximum
	nCnts <- sum( callTable)
	SANE_MAX <- 200
	if ( nCnts > SANE_MAX) {
		f <- SANE_MAX / nCnts
		callTable <- round( callTable * f)
		callTable <- callTable[ callTable > 0]
	}
	if ( length( callTable) < 2) return(1)

	# first make sure all of my terms are in the big table
	where <- MATCH( names( callTable), names( callSumsTable), nomatch=0)
	if ( any( where == 0)) {
		extras <- names(callTable)[where == 0]
		smlPop <- rep( 0, times=length(extras))
		smlPop <- callTable[ where == 0]
		names(smlPop) <- extras
		callSumsTable <- mergeTables( callSumsTable, as.table(smlPop))
	}

	# make sure the genome call is in our table
	if ( ! (GEN_CALL %in% names(callTable))) { 
		callTable <- mergeTables( callTable, as.table( c( ","=0)))
	}

	# now trim the big table to the intersection of me and the big table
	where <- MATCH( names( callSumsTable), names( callTable), nomatch=0)
	callSumsTable <- callSumsTable[ where]

	popProbs <- callSumsTable / sum( callSumsTable, na.rm=T)
	ans <- suppressWarnings( chisq.test( callTable, p=popProbs, simulate=F))
	return( sqrt(ans$p.value))
}


# integrate all the counts for all base calls into one total counts table for the population (chromosome)

callTable.totalSum <- function( callTables, verbose=FALSE) {

	if ( ! is.list( callTables) || length(callTables) < 2) {
		cat( "\nExpected a list of length > 1")
		return(NULL)
	}
	N <- length(callTables)
	
	# do via merge by 2's
	if (verbose) cat( "  summary..")
	while ( N > 1) {
		who <- seq.int( 1, N, by=2)
		newN <- length(who)
		if ( who[ newN] == N) {
			odd <- TRUE
			length(who) <- newN <- newN - 1
		} else { 
			odd <- FALSE
		}
		newTables <- lapply( who, function(x) {
					return( mergeTables( callTables[[x]], callTables[[x+1]]))
					})
		if (odd) newTables[[ newN + 1]] <- callTables[[N]]

		callTables <- newTables
		N <- length(callTables)
	}
	return( callTables[[1]])
}


`MPU.getIndelDetails` <- function( flips, curMPU) {

	indelText <- rep.int("", nrow(flips))
	whoIndel <- vector()
	bases <- type <- vector()
	ans <- list( "nIndels"=0, "who"=whoIndel, "bases"=vector(), "type"=vector(), "indelText"=indelText)

	# the 'flips' matrix has colnames of {',',A,C,G,T,N,Indel}
	whoMax <- apply( flips, MARGIN=1, WHICH.MAX)
	INDEL_COLUMN <- WHICH( colnames(flips) == "Indel")
	nIndels <- sum( whoMax == INDEL_COLUMN)
	if ( nIndels > 0) {

		# extract the details for these
		whoIndel <- WHICH( whoMax == INDEL_COLUMN)
		locs <- as.integer( rownames( flips)[whoIndel])
		where <- MATCH( locs, curMPU$POSITION)
		ref <- curMPU$REF_BASE[where]
		call <- curMPU$CALL_BASE[where]
	
		# the context determines insertion from deletion
		bases <- rep.int( "", length(whoIndel))
		type <- ifelse( NCHAR(ref) > NCHAR(call), "delete", "insert")
		isIns <- WHICH( type == "insert")
		isDel <- WHICH( type == "delete")
		# for insertions, the 'call' has both what was there and the extra bases
		bases[isIns] <- call[isIns]
		bases[isDel] <- call[isDel]
	}

	# also get the top 3 indel calls from any base that had any indels at all
	who2 <- WHICH( flips[ ,INDEL_COLUMN] > 0)
	if ( length( who2)) {
		top3indels <- SAPPLY( curMPU$BASE_TABLE[who2], function( txt) {
					terms <- strsplit( txt, split=";", fixed=T)[[1]]
					# reference and SNPs have a easy format
					notIndel <- grep( "^[,ACGT]=", terms)
					isIndel <- setdiff( 1:length(terms), notIndel)
					if ( length(isIndel)) {
						nshow <- min( length(isIndel), 3)
						return( PASTE( terms[isIndel[1:nshow]], collapse=";"))
					}
					return( "")
				})
		indelText[ who2] <- top3indels
	}

	ans <- list( "nIndels"=nIndels, "who"=whoIndel, "bases"=bases, "type"=type, "indelText"=indelText)
	return( ans)
}

