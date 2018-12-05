# pipe.ChIPseq.R

# do the A to Z alignment pipeline on a ChIPseq sample

`pipe.ChIPseq` <- function( sampleID=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt",
				...) {

	checkSampleName( sampleID)
	target <- getAndSetTarget( optionsFile, verbose=T)

	# note the time we start
	startTime <- proc.time()

	# setup...
	ans0 <- pipe.PreAlignTasks( sampleID, annotationFile, optionsFile, dataType="ChIP-seq")

	# do the alignment
	ans1 <- pipe.ChIPalignment( sampleID, annotationFile, optionsFile)

	# make pileups but not a transcriptome
	ans2 <- pipe.AlignToWig( sampleID, annotationFile, optionsFile, dataType="ChIP-seq")

	# do the peak picking...
	ans3 <- pipe.ChIPpeaks( sampleID, annotationFile, optionsFile, ...)

	# make any extra wanted display files
	ans4 <- pipe.PostAlignTasks( sampleID, annotationFile, optionsFile, dataType="ChIP-seq")

	# done...   report total time used
	cat( "\n", verboseOutputDivider, "\n\nChIP-seq PIPELINE DONE:   ", sampleID, "\n\n")
	print( elapsedProcTime( startTime, proc.time(), N=ans1$nReadsIn))

	return()
}

