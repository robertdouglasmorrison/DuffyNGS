# pipe.RNAseq.R

# do the A to Z alignment pipeline on a RNAseq sample

`pipe.RNAseq` <- function( sampleID=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt") {

	checkSampleNames( sampleID)
	target <- getAndSetTarget( optionsFile, sampleID=sampleID, annotationFile=annotationFile)

	# note the time we start
	startTime <- proc.time()

	# setup...
	ans0 <- pipe.PreAlignTasks( sampleID, annotationFile, optionsFile, dataType="RNA-seq")

	# do the alignment
	ans1 <- pipe.RNAalignment( sampleID, annotationFile, optionsFile)

	# build the wiggle track data
	ans2 <- pipe.AlignToWig( sampleID, annotationFile, optionsFile, dataType="RNA-seq")

	# make transcriptomes for all species...
	ans3 <- pipe.Transcriptome( sampleID, annotationFile, optionsFile, dataType="RNA-seq")

	# make any extra wanted display files
	ans4 <- pipe.PostAlignTasks( sampleID, annotationFile, optionsFile, dataType="RNA-seq")

	# done...   report total time used
	cat( "\n", verboseOutputDivider, "\n\nRNA-seq PIPELINE DONE:   ", sampleID, "\n\n")
	print( elapsedProcTime( startTime, proc.time(), N=ans1$nReadsIn))

	return()
}

