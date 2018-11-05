# pipe.RIPseq.R

# do the A to Z alignment pipeline on a RIPseq sample

`pipe.RIPseq` <- function( sampleID=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt",
				...) {

	target <- getAndSetTarget( optionsFile, verbose=T)

	# note the time we start
	startTime <- proc.time()

	# setup...
	ans0 <- pipe.PreAlignTasks( sampleID, annotationFile, optionsFile, dataType="RIP-seq")

	# do the alignment
	ans1 <- pipe.RIPalignment( sampleID, annotationFile, optionsFile)

	# make pileups but not a transcriptome
	ans2 <- pipe.AlignToWig( sampleID, annotationFile, optionsFile, dataType="RIP-seq")

	# make transcriptomes for all species...
	ans3 <- pipe.Transcriptome( sampleID, annotationFile, optionsFile, dataType="RIP-seq")

	# do the peak picking...
	ans3 <- pipe.RIPpeaks( sampleID, annotationFile, optionsFile, ...)

	# make any extra wanted display files
	ans4 <- pipe.PostAlignTasks( sampleID, annotationFile, optionsFile, dataType="RIP-seq")

	# done...   report total time used
	cat( "\n", verboseOutputDivider, "\n\nRIP-seq PIPELINE DONE:   ", sampleID, "\n\n")
	print( elapsedProcTime( startTime, proc.time(), N=ans1$nReadsIn))

	return()
}

