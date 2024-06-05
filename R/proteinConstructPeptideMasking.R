# protienConstructPeptideMasking.R -- functions to allow masking off of some pileup framents, typically from Introns


`IntronMaskSetup` <- function( intronMaskFasta, refAA, geneName, path=".", verbose=TRUE) {

	# we were given some intron masking info, in the form of a FASTA file or object
	
	# step 1:  let the FASTA be a filename, or an object with 'desc','seq' elements.
	# and if a filename, let a local copy in the results folder take precedent 
	intronFA <- NULL
	if ( is.character(intronMaskFasta)) {
		localFastaFile <- file.path( path, basename( intronMaskFasta))
		if ( file.exists( localFastaFile)) {
			cat( "\nUsing sample-specifc Intron Mask FASTA file: ", localFastaFile)
			intronFA <- loadFasta( localFastaFile, verbose=F)
		} else {
			cat( "\nUsing global Intron Mask FASTA file: ", intronMaskFasta)
			intronFA <- loadFasta( intronMaskFasta, verbose=F)
		}
	}
	if ( is.list(intronMaskFasta)) {
		if ( all( c("desc","seq") %in% names(intronMaskFasta))) {
			cat( "\nUsing explicit Intron Mask FASTA object ")
			intronFA <- intronMaskFasta
		}
	}
	if ( is.null( intronFA)) {
		cat( "\nError:  failed to turn 'intronMaskFasta' argument into masking data object")
		return( NULL)
	}
	
	# step 2:  verify that all the masks meet our requirements.
	# A)  length of each name (the reference sequence) matches each mask sequence (after stripping blank characters)
	maskDesc <- gsub( " ", "", toupper(intronFA$desc))
	maskSeq <- gsub( " ", "", toupper(intronFA$seq))
	maskLen <- nchar( maskDesc)
	nMask <- length( maskSeq)
	if ( ! all( maskLen == nchar(maskSeq))) {
		cat( "\nError:  length of each reference motif must exactly match length of mask sequence\n")
		probs <- which( maskLen != nchar(maskSeq))
		badDF <- data.frame( "Reference.Motif"=maskDesc[probs], "Mask.Seq"=maskSeq[probs])
		print(badDF)
		return( NULL)
	}
	# B)  at each amino acid, the mask AA should never be an exact match to the expected reference AA
	badCh <- FALSE
	for ( i in 1:nMask) {
		myDesc <- maskDesc[i]
		mySeq <- maskSeq[i]
		for ( j in 1:maskLen[i]) {
			refCh <- substr( myDesc, j, j)
			maskCh <- substr( mySeq, j, j)
			if ( refCh == maskCh) {
				cat( "\nWarning:  mask AA exactly matches reference AA at location", j, refCh, maskCH)
				cat( "\n  Double check mask for motif: ", myDesc, " <-> ", mySeq)
				badCh <- TRUE
			}
		}
	}
	if ( badCh) {
		cat( "\nError:  correct mask sequences to not match expected reference AAs.")
		return( NULL)
	}
	
	# step 3:  find the one best location for each mask
	maskStart <- maskScore <- rep.int(NA,nMask)
	for ( i in 1:nMask) {
		pa <- pairwiseAlignment( maskDesc[i], refAA, type='global-local', scoreOnly=FALSE)
		maskScore[i] <- scorePerAA <- round( score(pa) / nchar(maskDesc[i]), digits=3)
		if ( scorePerAA < 1) {
			cat( "\nWarning:  mask name ", maskDesc[i], "has no good scoring match to Reference AA")
			next
		} 
		maskStart[i] <- start( subject(pa))
	}
	maskInfo <- data.frame( "Reference.Motif"=maskDesc, "Mask.Seq"=maskSeq, "Mask.Start"=maskStart, 
				"Mask.Len"=maskLen, "ScorePerAA"=maskScore, stringsAsFactors=F)
	
	if (verbose) {
		cat( "\nParsed ", nrow(maskInfo), "Intron Masking fragments for gene: ", geneName, "\n")
		print( maskInfo)
	}
	return( maskInfo)
}


`IntronMaskDisplay` <- function( maskInfo, cex=1) {

	# render the intron masking key on the pileups image
	if ( is.null(maskInfo)) return(NULL)
	nMask <- nrow( maskInfo)
	if ( ! nMask) return(NULL)
	text( 0, -1.75, "Masks:", cex=cex, col=2, pos=2)
	for ( i in 1:nMask) {
		myStart <- maskInfo$Mask.Start[i]
		if (is.na(myStart)) next
		myLen <- maskInfo$Mask.Len[i]
		myStop <- myStart + myLen - 1
		myChars <- strsplit( maskInfo$Mask.Seq[i], split="")[[1]]
		text( myStart:myStop, rep.int(-1.75,myLen), myChars, cex=cex, col=2)
	}
}


`IntronMaskAdjustPileups` <- function( obj, maskInfo, verbose=TRUE) {

	# given the full pileup result, extract the call details
	aaCalls <- obj$AA_Calls
	aaWeights <- obj$AA_Weights
	
	if ( is.null(maskInfo)) return(NULL)
	nMask <- nrow( maskInfo)
	if ( ! nMask) return(NULL)
	nMaskedBits <- 0
	if (verbose) cat( "\nDoing Intron Masking with", nMask, "masks.")
	for ( i in 1:nMask) {
		myStart <- maskInfo$Mask.Start[i]
		if (is.na(myStart)) next
		myLen <- maskInfo$Mask.Len[i]
		myStop <- myStart + myLen - 1
		myChars <- strsplit( maskInfo$Mask.Seq[i], split="")[[1]]
		myProteinLocs <- myStart:myStop
		for ( j in 1:myLen) {
			myMaskChar <- myChars[j]
			myColumn <- myProteinLocs[j]
			if (myColumn > ncol(aaCalls)) next
			myCalls <- aaCalls[ , myColumn]
			myWeights <- aaWeights[ , myColumn]
			toMask <- which( myCalls == myMaskChar)
			if ( length(toMask)) {
				myCalls[toMask] <- ""
				myWeights[toMask] <- 0
				aaCalls[ , myColumn] <- myCalls
				aaWeights[ , myColumn] <- myWeights
				nMaskedBits <- nMaskedBits + length(toMask)
			}
		}
	}
	
	# all done
	if (nMaskedBits) {
		totalBits <- nrow(aaCalls) * (ncol(aaCalls)-30)		# the extra tail at C-term end
		pctMask <- round( nMaskedBits * 100 / totalBits, digits=2)
		if (verbose) cat( "\nMasking results:  removed ", nMaskedBits, " intron AA calls from ", totalBits, 
			" pileup locations (", pctMask, "%)\n", sep="")
		# overwrite if we modified anything
		obj$AA_Calls <- aaCalls
		obj$AA_Weights <- aaWeights
	}
	return(obj)
}
