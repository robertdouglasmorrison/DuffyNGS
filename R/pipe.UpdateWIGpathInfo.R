# pipe.UpdateWIGpathInfo.R -- change the embedded file path info in a WIG data structure

`pipe.UpdateWIGpathInfo` <- function( sampleID,  old.results.path, 
				annotationFile="Annotation.txt", optionsFile="Options.txt",
				verbose=TRUE) {

	# make folders for all results...
	results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	target <- getAndSetTarget( optionsFile, sampleID=sampleID, annotationFile=annotationFile)
	prefixes <- getCurrentTargetFilePrefix()
	
	thisFolder <- file.path( results.path, "wig")

	# make sure the final path separator is not on either 'old' or 'new' to keep it correct during substitution
	results.path <- sub( "/$", "", results.path)
	old.results.path <- sub( "/$", "", old.results.path)

	for ( prefix in prefixes) {
		wigfile <- file.path( thisFolder, paste( sampleID, prefix, "WIG.rda", sep="."))

		who <- load( wigfile)
		if( who != "wiggles") next
		myInfo <- wiggles$Info
		if (is.null( myInfo)) next
		mySubFolder <- wiggles$SubWigFolder
		if (is.null( mySubFolder)) next
		cat( "\nRenaming file paths in WIG file: ", basename(wigfile))

		# first verify that the substitution will work
		if ( regexpr( old.results.path, myInfo$FileName[1], fixed=T) < 1) {
			cat( "\nError:  'old.results.path' not seen in WIG path info")
			cat( "\nold.results.path expected:  ", old.results.path)
			cat( "\nfile path currently in WIG: ", dirname( myInfo$FileName[1]))
			return()
		}

		myInfo$FileName <- sub( old.results.path, results.path, myInfo$FileName, fixed=T)
		mySubFolder <- sub( old.results.path, results.path, mySubFolder, fixed=T)
		wiggles$Info <- myInfo
		wiggles$SubWigFolder <- mySubFolder

		save( wiggles, file=wigfile)
	}
	return()
}
