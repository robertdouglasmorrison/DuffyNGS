# extractPipelineSummaryDetails.R -- pull the key values out of the summary text file

`extractPipelineSummaryDetails` <- function( sampleID, optionsFile="Options.txt", 
			results.path=NULL, verbose=FALSE) {

	nRaw <- nNoHit <- nRibo <- nGenomic <- nSplice <- vector()
	pRaw <- pNoHit <- pRibo <- pGenomic <- pSplice <- vector()
	nout <- 0

	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=FALSE)
	}
	target <- getAndSetTarget( optionsFile)
	speciesIDs <- getCurrentTargetSpecies()
	prefixes <- getCurrentTargetFilePrefix()
	doWIGs <- ( (NSP <- length(speciesIDs)) > 1)
	if (doWIGs) {
		speciesCnts <- matrix( 0, nrow=length(sampleID), ncol=NSP)
		speciesPcts <- matrix( 0, nrow=length(sampleID), ncol=NSP)
		colnames(speciesCnts) <- paste( "N", speciesIDs, sep="_")
		colnames(speciesPcts) <- paste( "P", speciesIDs, sep="_")
	}


	extractReadCountFromTextLine <- function( txt) {

		if ( is.null( txt) || length( txt) < 1) return( 0)

		txt <- txt[1]

		n <- sub( "(^.+Reads:)([ \t]+)([0-9][0-9,]*)(.*)", "\\3", txt)
		n <- gsub( ",", "", n)

		return( as.numeric(n))
	}


	for ( s in sampleID) {

		file <- paste( s, "pipeline.Summary.txt", sep=".")
		file <- file.path( results.path, "summary", file)
		nraw <- nnohit <- nribo <- ngenomic <- nsplice <- NA
		praw <- pnohit <- pribo <- pgenomic <- psplice <- NA

		if ( file.exists( file)) {
			txt <- readLines( file)
			rawLine <- grep( "N_Raw Reads", txt)
			nohitLine <- grep( "N_NoHit Reads", txt)
			riboLine <- grep( "All Ribo Reads", txt)
			genomicLine <- grep( "All Genomic Reads", txt)
			spliceLine <- grep( "All Splice Reads", txt)

			nraw <- extractReadCountFromTextLine( txt[rawLine])
			nnohit <- extractReadCountFromTextLine( txt[nohitLine])
			nribo <- extractReadCountFromTextLine( txt[riboLine])
			ngenomic <- extractReadCountFromTextLine( txt[genomicLine])
			nsplice <- extractReadCountFromTextLine( txt[spliceLine])

			# if any failed to extract # from text summary, look in files
			if ( is.na( nnohit)) {
				f <- file.path( results.path, "fastq", paste( s, "noHits.fastq", sep="."))
				nnohit <- getFileLineCount( f, sampleID=s, verbose=F, what="lineCount") / 4
			}
		} else {
			# this may be a paired end strand specific sample with 2 sets of summary
			nraw <- nnohit <- nribo <- ngenomic <- nsplice <- 0
			for (spair in paste( s, 1:2, sep="_")) {
				file <- paste( spair, "pipeline.Summary.txt", sep=".")
				file <- file.path( results.path, "summary", file)
				if ( file.exists( file)) {
					txt <- readLines( file)
					rawLine <- grep( "N_Raw Reads", txt)
					nohitLine <- grep( "N_NoHit Reads", txt)
					riboLine <- grep( "All Ribo Reads", txt)
					genomicLine <- grep( "All Genomic Reads", txt)
					spliceLine <- grep( "All Splice Reads", txt)

					nraw <- nraw + extractReadCountFromTextLine( txt[rawLine])
					nnohit <- nnohit + extractReadCountFromTextLine( txt[nohitLine])
					nribo <- nribo + extractReadCountFromTextLine( txt[riboLine])
					ngenomic <- ngenomic + extractReadCountFromTextLine( txt[genomicLine])
					nsplice <- nsplice + extractReadCountFromTextLine( txt[spliceLine])
				}
			}
		}
		if (verbose) cat( "\r", s, nraw, nnohit, nribo, ngenomic, nsplice)

		nout <- nout + 1
		nRaw[nout] <- nraw
		nNoHit[nout] <- nnohit
		nRibo[nout] <- nribo
		nGenomic[nout] <- ngenomic
		nSplice[nout] <- nsplice
		pNoHit[nout] <- as.percent( nnohit, big.value=nraw, percentSign=F)
		pRibo[nout] <- as.percent( nribo, big.value=nraw, percentSign=F)
		pGenomic[nout] <- as.percent( ngenomic, big.value=nraw, percentSign=F)
		pSplice[nout] <- as.percent( nsplice, big.value=nraw, percentSign=F)

		if (doWIGs) {
			cnts <- rep( NA, times=NSP)
			for( j in 1:NSP) {
				file <- paste( s, prefixes[j], "WIG.rda", sep=".")
				file <- file.path( results.path, "wig", file)
				if ( ! file.exists(file)) next
				load( file)
				cnts[j] <- wiggles$Info$TotalReads
			}
			pcts <- ( cnts * 100 / sum(cnts, na.rm=T))
			speciesCnts[ nout, ] <- round( cnts)
			speciesPcts[ nout, ] <- round( pcts, digits=2)
		}
	}

	out <- data.frame( "SampleID"=sampleID, "RawReads"=nRaw, 
			 "N_RiboClear"=nRibo, "P_RiboClear"=pRibo,
			 "N_Genomic"=nGenomic, "P_Genomic"=pGenomic,
			 "N_Splice"=nSplice, "P_Splice"=pSplice, 
			 "N_NoHit"=nNoHit, "P_NoHit"=pNoHit,
			 stringsAsFactors=FALSE)
	rownames(out) <- 1:nrow(out)

	if ( doWIGs) out <- cbind( out, speciesCnts, speciesPcts, stringsAsFactors=FALSE)

	return( out)
}
