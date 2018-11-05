# pipe.RemoveReadsFromNoHits.R


`pipe.RemoveReadsFromNoHits` <- function( sampleID=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt",
				results.path=NULL, indexSet=c("IlluminaAdapters_idx", "Ecoli_O157H7_idx"), 
				indexMode=c( "local", "genomic"), saveBAMs=TRUE, verbose=TRUE) {

	if ( is.null( results.path)) results.path <- getOptionValue( optionsFile, "results.path", 
				notfound=".", verbose=F)

	# make sure we have the existing no hits data
	noHits.path <- file.path( results.path, "fastq")
	noHitsFile <- file.path( noHits.path, paste( sampleID, "noHits.fastq.gz", sep="."))
	if ( ! file.exists( noHitsFile)) {
		cat( "\nNoHits file not found: ", noHitsFile)
		return(NULL)
	}

	# can be given a set of more than one index to align against
	if ( length( indexSet) != length( indexMode)) {
		cat( "\nLengths of Index Set and Index Mode mismatch..")
		stop()
	}
	N_Index <- length(indexSet)
	out <- vector( mode='list', length=N_Index)

	for ( i in 1:N_Index) {
		myIndex <- indexSet[i]
		myIndexName <- sub( "_idx$", "", myIndex)
		myMode <- indexMode[i]
		if ( ! (myMode %in% c( "genomic", "local"))) {
			cat( "\nIndex Mode must be one of 'genomic', 'local'")
			stop()
		}

		# align against this 'unwanted' genome
		bamOutfile <- paste( sampleID, "ReadHits", myIndexName, "bam", sep=".")
		newNoHitsFile <- paste( sampleID, "NewNoHits.fastq.gz", sep=".")

		if ( myMode == "genomic") {
			myPolicy <- getOptionValue( optionsFile, "GenomicAlignmentPolicy", verbose=F)
		} else {
			myPolicy <- " --very-sensitive-local  --score-min L,50,0 "
		}
		ans <- fastqToBAM( noHitsFile, bamOutfile, sampleID=sampleID, optionsFile=optionsFile, 
				noHitsFile=newNoHitsFile, alignIndex=myIndex, alignPolicy=myPolicy, verbose=TRUE)
		nNoHits <- ans$NoHitReads
		nHits <- ans$UniqueReads + ans$MultiReads

		# put this new smaller no hits results and line count back where it goes...
		file.rename( newNoHitsFile, noHitsFile)
		quickFileLineCountRecord( noHitsFile, sampleID, lineCount=nNoHits*4, readCount=nNoHits)

		# either store or delete the new BAM file of hits
		if (saveBAMs) {
			file.rename( bamOutfile, file.path( results.path, "align", bamOutfile))
		} else {
			file.delete( bamOutfile)
		}

		cat( "\n\nReads that aligned to: ", myIndexName, " \t", nHits)
		cat( "\nReads still 'NoHits':    ", " \t", nNoHits)

		# done with this index
		out[[i]] <- ans
	}
	return( out)
}
