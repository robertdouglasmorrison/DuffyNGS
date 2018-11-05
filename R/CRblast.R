# CRblast.R

# Blast the CR constructs to see what they are


CRblaster <- function( sampleID, crIDs=NULL, nBest=10, xmlOutFile=NULL, doBlast=TRUE,
			optionsFile="Options.txt", blastIndex="nt",
			evalue=1.0, wordsize=8, nHSP=3, verbose=TRUE) {

	if ( doBlast && verbose) cat( "\n\nBlasting Consensus Reads:  ", crIDs, "\n")

	ansOut <- vector( mode="list")

	# we will blast to XML output file, then extract what we want
	if ( is.null( xmlOutFile)) xmlOutFile <- paste( sampleID, "blastOutput.xml", sep=".")
	if ( doBlast) file.delete( xmlOutFile)

	# get the Blast program and index name
	optT <- readOptionsTable( optionsFile)
	blastProgram <- getOptionValue( optT, "blastProgram", notfound="./blastall")
	blastIndexPath <- getOptionValue( optT, "blastIndex.path", notfound=".")
	blastIndex <- file.path( blastIndexPath, blastIndex)

	# verify...
	if (doBlast) {
		if ( ! file.exists( blastProgram)) {
			cat( "\n\nBlast Program not found.  Tried: ", blastProgram, "\n")
			return( NULL)
		}
		if ( ! file.exists( blastIndexPath)) {
			cat( "\n\nBlast Index folder not found.  Tried: ", blastIndexPath, "\n")
			return( NULL)
		}
	}

	# decide which CR constructs to Blast
	if ( is.null( crIDs)) {
		whoLen <- CRplotter( seconds=NULL, plotOrder="length", N=nBest, doPlot=FALSE)
		whoCnt <- CRplotter( seconds=NULL, plotOrder="count", N=nBest, doPlot=FALSE)
		crIDs <- base::sort( base::union( whoLen, whoCnt))
	}

	if (doBlast) {
		# turn the given CR sequences into a .fasta on standard input
		stdin <- vector()
		for ( i in crIDs) {
			stdin <- c( stdin, base::paste( ">cr_",i, sep=""))
			stdin <- c( stdin, base::paste( CRT_List[[i]]$bases, collapse=""))
		}
	
		# allow for a retry if no hits
		useEvalue <- evalue
		useWordsize <- wordsize
		for ( k in 1:3) {
			# OK, build the command line
			blastFlags <- paste( " -evalue ", useEvalue, " -word_size ", useWordsize,
					" -outfmt 5 -dust no ", " -num_alignments ", nHSP, " -num_threads 4 ")
			cmdLine <- paste( blastProgram, " -db ", blastIndex, " -out ", xmlOutFile, 
					blastFlags, sep=" ")
			if ( verbose) cat( "\nTry: ",k, "\tBLAST command line:\n", cmdLine)
			
			system( cmdLine,  input=stdin)
		
			# extract the answer
			ansOut <- extractBlastXMLdetails( xmlOutFile, crIDs=crIDs)
			if ( length( ansOut) < 1) {
				cat( "\n\nNo Hits Found!!   used Evalue cutoff of: ", useEvalue, "\n")
				# remove any old results, and run it again!
				file.delete( xmlOutFile)
			} else {
				break
			}
			useEvalue <- useEvalue * 100
			useWordsize <- if ( useWordsize > 7) useWordsize - 2 else useWordsize
		}
	} else {
		cat( "\nUsing previously calculated Blast answer..")
		ansOut <- extractBlastXMLdetails( xmlOutFile, crIDs=crIDs)
	}

	return(ansOut)
}


`printCRblastOutput` <- function( ans) {

	if ( length( ans) < 1) return()

	for ( i in 1:length( ans)) {
		cat( "\n\nCR id = ", ans[[i]]$id, "\tLength =", CRT_List[[ N[i]]]$len, "bases\n")
		cat( "\nN_Distinct_Reads =", length( CRT_List[[ N[i]]]$USRid), 
			"\tN_Total_Reads =", CRT_List[[ N[i]]]$totReads, "\tPct_All_Reads =",
			formatC( (CRT_List[[ N[i]]]$totReads * 100 / USR_GrandTotal), digits=3, 
			format="f"),"\n")

		if ( ans[[i]]$n < 1) {
			cat( "\nNo Hits Found!!   used Evalue cutoff of: ", evalue, "\n\n")
			next
		}

		for( j in 1:ans[[i]]$n) {
			cat( "\nScore:", round(ans[[i]]$score[j]), "\tEval:", formatC(ans[[i]]$evalue[j], 
				digits=2,format="e"), "\t\tAccession:", 
				ans[[i]]$accession[j], "\n", ans[[i]]$definition[j], "\n")
			cat( ans[[i]]$match[j], "\n")
		}
	}
}


extractBlastXMLdetails <- function( filein, crIDs) {

	# we want to extract the details that show the 'best hits'...
	out <- vector( mode="list")
	txt <- readLines( con=filein)
	iterDefLines <- grep( "Iteration_query-def", txt, fixed=TRUE)
	if ( length( iterDefLines) < 1) return( out)

	for( i in 1:length(iterDefLines)) {
		if ( i < length( iterDefLines)) {
			lastLine <- iterDefLines[i+1]
		} else {
			lastLine <- length( txt)
		}
		out[[i]] <- extractOneXMLiterationSet( txt[ iterDefLines[i] : lastLine] )
	}

	# some CRs may not return anything, but try to make the answer more complete
	idsOut <- sapply( out, function(x) x$id)
	finalOut <- vector( mode="list")
	for ( j in 1:length(crIDs)) {
		i <- crIDs[j]
		thisID <- paste( "cr_",i,sep="")
		where <- match( thisID, idsOut, nomatch=0)
		if ( where > 0) {
			finalOut[[j]] <- out[[where]]
		} else {
			finalOut[[j]] <- list( "n"=0, "id"=thisID, "accession"=NA, "definition"=NA, "evalue"=NA, 
						"score"=NA, "match"=NA, "hitFrom"=NA, "hitTo"=NA)
		}
	}

	return( finalOut)
}


extractOneXMLiterationSet <- function( txt) {

	# given a subset of text lines that contain one XML blast output alignment...
	id <- sub( "( *<Iteration_query-def>)(.+)(</Iteration_query-def>)", "\\2", txt[1])
	hitLines <- grep( "<Hit>", txt, fixed=TRUE)
	nHits <- length( hitLines)
	if ( nHits < 1) return( list( "n"=0, "id"=id))

	hitACClines <- grep( "<Hit_accession>", txt, fixed=TRUE)
	hitACCs <- sub( "( *<Hit_accession>)(.+)(</Hit_accession>)", "\\2", txt[ hitACClines])
	hitDEFlines <- grep( "<Hit_def>", txt, fixed=TRUE)
	hitDEFs <- sub( "( *<Hit_def>)(.+)(</Hit_def>)", "\\2", txt[ hitDEFlines])
	
	# use just the first HSP for each hit, static offsets from hit line...
	hitEvalue <- hitScore <- hitMatch <- hitFrom <- hitTo <- vector()
	for ( i in 1:length( hitLines)) {
		line <- hitLines[i]
		hitScore[i] <- as.numeric( sub( "( +<Hsp_bit-score>)(.+)(</Hsp_bit-score>)", "\\2", txt[line+9]))
		hitEvalue[i] <- as.numeric( sub( "( +<Hsp_evalue>)(.+)(</Hsp_evalue>)", "\\2", txt[line+11]))
		qfrom <- as.integer( sub( "( +<Hsp_query-from>)(.+)(</Hsp_query-from>)", "\\2", txt[line+12]))
		qto <- as.integer( sub( "( +<Hsp_query-to>)(.+)(</Hsp_query-to>)", "\\2", txt[line+13]))
		hfrom <- as.integer( sub( "( +<Hsp_hit-from>)(.+)(</Hsp_hit-from>)", "\\2", txt[line+14]))
		hto <- as.integer( sub( "( +<Hsp_hit-to>)(.+)(</Hsp_hit-to>)", "\\2", txt[line+15]))
		hitFrom[i] <- hfrom
		hitTo[i] <- hto

		# there may be a 'gaps' line... there may not be...
		hasGaps <- ( regexpr( "Hsp_gaps", txt[ line+20], fixed=TRUE) > 0)
		if ( hasGaps) line <- line + 1

		seq <- sub( "( +<Hsp_qseq>)(.+)(</Hsp_qseq>)", "\\2", txt[line+21])
		qtxt <- base::paste( formatC( qfrom, digits=12, format="d"), seq, formatC( qto, digits=12, 
				format="d"), sep=" ")
		seq <- sub( "( +<Hsp_hseq>)(.+)(</Hsp_hseq>)", "\\2", txt[line+22])
		htxt <- base::paste( formatC( hfrom, digits=12, format="d"), seq, formatC( hto, digits=12, 
				format="d"), sep=" ")
		seq <- sub( "( +<Hsp_midline>)(.+)(</Hsp_midline>)", "\\2", txt[line+23])
		mtxt <- base::paste( "             ", seq, "             ", sep=" ")
		hitMatch[i] <- base::paste( qtxt, mtxt, htxt, sep="\n")
	}

	out <- list( "n"=nHits, "id"=id, "accession"=hitACCs, "definition"=hitDEFs, "evalue"=hitEvalue, 
			"score"=hitScore, "match"=hitMatch, "hitFrom"=hitFrom, "hitTo"=hitTo)

	return( out)
}


`summarizeCRblastOutput` <- function( ans, fileout="summaryBlast.txt") {

	# given the results of a CR blast, turn it into a table text file...
	N_CR <- length( ans)
	file.delete( fileout)
	if ( N_CR < 1) return()

	# build one DF of all the CR blast results...
	idSet <- lenSet <- totReadsSet <- pctReadsSet <- scoreSet <- evalueSet <- hitSet <- hitTextSet <- vector()
	chrSet <- posSet <- endSet <- seqSet <- vector()

	nout <- 0
	for ( j in 1:N_CR) {
		thisAns <- ans[[j]]
		if ( is.null( thisAns)) next
		nout <- nout + 1
		idSet[nout] <- thisAns$id
		# force every thing to get an entry
		lenSet[nout] <- totReadsSet[nout] <- scoreSet[nout] <- evalueSet[nout] <- NA
		seqSet[nout] <- chrSet[nout] <- hitSet[nout] <- ""
		hitTextSet[nout] <- "No Blast Hits Found"
		posSet[nout] <- endSet[nout] <- NA
		# OK, fill in what we found
		thisIDnum <- as.integer( sub( "cr_", "", thisAns$id))
		thisCRT <- CRT_List[[ thisIDnum]]
		lenSet[nout] <- thisCRT$len
		totReadsSet[nout] <- thisCRT$totReads
		seqSet[nout] <- base::paste( thisCRT$bases, collapse="")
		if ( thisAns$n < 1) next
		scoreSet[nout] <- thisScore <- thisAns$score[1]
		evalueSet[nout] <- thisAns$evalue[1]

		# some parsing only is doable for a true hit...
		if ( ! is.na( thisScore)) {
			hitTxt <- massageHitDefinition( thisAns$definition)
			hitSet[nout] <- hitTxt$brief
			hitTextSet[nout] <- hitTxt$full
	
			# see if we got a chromosomal hit we can use
			if ( hitTxt$chromosome > 0) {
				chrSet[ nout] <- hitTxt$chromosome
				posSet[ nout] <- thisAns$hitFrom[ hitTxt$textLineNumber]
				endSet[ nout] <- thisAns$hitTo[ hitTxt$textLineNumber]
			}
		}
	}

	pctReadsSet <- round( totReadsSet * 100 / USR_GrandTotal, digits=4)
	out <- data.frame( idSet, lenSet, totReadsSet, pctReadsSet, hitSet, scoreSet, evalueSet, chrSet, 
				posSet, endSet, seqSet, hitTextSet, stringsAsFactors=FALSE)
	colnames( out) <- c( "CR_ID", "N_BASES", "N_READS", "PCT_READS", "BLAST_HIT", "SCORE", "E_VALUE", 
				"CHROMO", "START", "STOP", "CONSENSUS_SEQUENCE", 
				"TEXT_DESCRIPTION_OF_BEST_HITS")
	ord <- base::order( out$N_READS, decreasing=TRUE)
	out <- out[ ord, ]
	rownames( out) <- 1:nrow(out)
	write.table( out, file=fileout, sep="\t", quote=FALSE, row.names=FALSE)
	cat( "\nWrote CR blast summary: ", fileout)

	return( out)
}


massageHitDefinition <- function( txt) {

	# given a set of blast hit definitions, try to get one key nickname...
	full <- base::paste( 1:length(txt), txt, sep=")  ", collapse=" |   ")

	hasHsap <- any( regexpr( "Homo sapiens", txt) > 0)
	hasPfal <- any( regexpr( "Plasmodium falciparum", txt) > 0)
	hasMycoplasma <- any( regexpr( "Mycoplasma", txt) > 0)
	#hasRRNA <- any( (regexpr( "ribosomal RNA", txt) > 0) | (regexpr( "18S|28S|5\\.8S", txt, extended=TRUE) > 0))
	#hasGlobin <- any( regexpr( "globin", txt) > 0)
	#hasMitochon <- any( regexpr( "mitochondrion", txt) > 0)

	nickname <- ""
	if ( hasHsap) nickname <- "H.sap"
	if ( hasPfal) nickname <- "P.fal"
	if ( hasMycoplasma) nickname <- "Mycoplasma"
	#if ( hasRRNA) nickname <- paste( nickname, "(ribosomal RNA)")
	#if ( hasGlobin) nickname <- paste( nickname, "(globin)")
	#if ( hasMitochon) nickname <- paste( nickname, "(mitochondrion)")

	# try to draw out the Pfal geneID if possible
	if ( hasPfal) {
		isPfal <- ( regexpr( "Plasmodium falciparum", txt) > 0)
		tmps <- sub( "(^.+)(\\(PF)(.{4,9}\\))(.*$)", "\\2\\3", txt)
		best <- which( ( base::substr( tmps, 1,1) == "(") & isPfal)
		if ( length( best) > 0) {
			nickname <- base::paste( nickname, tmps[ best[1]])
		} else {
			tmps <- sub( "(^.+)(\\(MAL)(.{4,9}\\))(.*$)", "\\2\\3", txt)
			best <- which( ( base::substr( tmps, 1,1) == "(") & isPfal)
			if ( length( best) > 0) {
				nickname <- base::paste( nickname, tmps[ best[1]])
			}
		}
	}
	if ( !hasPfal && hasHsap) {
		isHsap <- (regexpr( "Homo sapiens", txt) > 0)
		tmps <- sub( "(^.+)(\\([A-Z]+.*\\))(.*$)", "\\2", txt)
		best <- which( ( base::substr( tmps, 1,1) == "(") & isHsap)
		if ( length( best) > 0) {
			nickname <- base::paste( nickname, tmps[ best[1]])
		}
	}

	# see if we have a 3D7 chromosomal hit that we could map
	chr <- chrTextLine <- 0
	who3D7 <- which( (regexpr( "3D7", txt, fixed=TRUE) > 0) & (regexpr( "chromosome", txt, fixed=TRUE) > 0))
	if ( length( who3D7) > 0) {
		line <- who3D7[1]
		chrText <- sub( "(^.+)(chromosome )([0-9]{1,2})(.*)", "\\3", txt[line])
		try( chrText <- as.integer( chrText))
		if ( ! is.na( chrText)) {
			chr <- chrText
			chrTextLine <- line
		}
	}

	# if we still have 'nothing', just get the first two words from the first line
	if ( nickname == "") {
		terms <- strsplit( txt[1], split=" ", fixed=TRUE)[[1]]
		if ( terms[1] == "PREDICTED:") terms <- terms[ 2: length(terms)]
		nickname <- base::paste( terms, collapse=" ")
		nickname <- sub( ",.+", "", nickname)
	}

	return( list( "brief"=nickname, "full"=full, "chromosome"=chr, "textLineNumber"=chrTextLine))
}
