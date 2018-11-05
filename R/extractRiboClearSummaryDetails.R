# extractRiboClearSummaryDetails.R -- pull the key values out of the summary text file

`extractRiboClearSummaryDetails` <- function( sampleID, optionsFile="Options.txt", 
			results.path=NULL, speciesID="Hs_grc", cleared=c("18S/28S", "Globin", "Mito_RNA")) {

	nRaw <- nRibo <- vector()
	pRaw <- pRibo <- vector()
	nCLEAR <- matrix( 0, nrow=length(sampleID), ncol=length(cleared))
	pCLEAR <- matrix( "", nrow=length(sampleID), ncol=length(cleared))
	colnames(nCLEAR) <- paste( "N", cleared, sep="_")
	colnames(pCLEAR) <- paste( "P", cleared, sep="_")
	nout <- 0

	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=FALSE)
	}


	extractReadCountFromTextLine <- function( txt) {

		if ( is.null( txt) || length( txt) < 1) return( 0)

		txt <- txt[1]

		n <- sub( "(^.+Reads:)([ \t]+)([1-9][0-9,]*)(.*)", "\\3", txt)
		n <- gsub( ",", "", n)

		return( as.numeric(n))
	}


	extractReadCountFromClearLines <- function( txt, speciesID, cleared="18S/28S") {

		out <- rep.int( 0, length(cleared))
		names(out) <- cleared

		if ( is.null( txt) || length( txt) < 1) return( out)

		# only for this SpeciesID
		keep <- grep( speciesID, txt)
		txt <- txt[ keep]
		if ( length( txt) < 1) return( out)

		for ( i in 1:length(txt)) {
		    onetxt <- txt[i]
		    for ( j in 1:length(cleared)) {
		    	idx <- regexpr( cleared[j], onetxt, fixed=T)
			if ( idx < 1) next
			nskip <- idx + nchar(cleared[j])
			subtxt <- substr( onetxt, nskip, nchar(onetxt))
			n <- sub( "(^[ \t]+)([1-9][0-9,]*)(.+)", "\\2", subtxt)
			n <- gsub( ",", "", n)
			out[ j] <- as.numeric(n)
		    }
		}

		return( out)
	}



	for ( s in sampleID) {

		file <- paste( s, "ribo.Summary.txt", sep=".")
		file <- file.path( results.path, "summary", file)
		nraw <- nribo <- NA
		praw <- pribo <- NA
		nclears <- rep.int( 0, length( cleared))

		if ( file.exists( file)) {
			txt <- readLines( file)
			rawLine <- grep( "^Raw Reads", txt)
			riboLine <- grep( "^Total Cleared Reads", txt)
			clearLines <- grep( "^Cleared", txt)

			nraw <- extractReadCountFromTextLine( txt[rawLine])
			nribo <- extractReadCountFromTextLine( txt[riboLine])
			nclears <- extractReadCountFromClearLines( txt[clearLines], speciesID=speciesID, 
						cleared=cleared)
		} else {
			# this may be a paired end strand specific sample with 2 sets of summary
			nraw <- nribo <- 0
			for (spair in paste( s, 1:2, sep="_")) {
				file <- paste( spair, "ribo.Summary.txt", sep=".")
				file <- file.path( results.path, "summary", file)
				if ( file.exists( file)) {
					txt <- readLines( file)
					rawLine <- grep( "^Raw Reads", txt)
					riboLine <- grep( "^Total Cleared Reads", txt)
					clearLines <- grep( "^Cleared", txt)

					nraw <- nraw + extractReadCountFromTextLine( txt[rawLine])
					nribo <- nribo + extractReadCountFromTextLine( txt[riboLine])
					nclears <- nclears + extractReadCountFromClearLines( txt[clearLines], 
								speciesID=speciesID, cleared=cleared)
				}
			}
		}

		nout <- nout + 1
		nRaw[nout] <- nraw
		nRibo[nout] <- nribo
		pRibo[nout] <- as.percent( nribo, big.value=nraw, percentSign=F)
		nCLEAR[ nout, ] <- nclears
	}
	nCLEAR <- nCLEAR[ 1:nout, ]
	for (j in 1:ncol(nCLEAR)) {
		pCLEAR[ , j] <- as.percent( nCLEAR[ ,j], big.value=nraw, percentSign=F)
	}

	out <- data.frame( "SampleID"=sampleID, "RawReads"=nRaw, 
			 "N_RiboClear"=nRibo, "P_RiboClear"=pRibo,
			 nCLEAR, pCLEAR, stringsAsFactors=FALSE)
	rownames(out) <- 1:nrow(out)

	return( out)
}
