# pipe.DNAseq.R

# do the A to Z alignment pipeline on a DNAseq sample

`pipe.DNAseq` <- function( sampleID=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt") {

	target <- getAndSetTarget( optionsFile, sampleID=sampleID, annotationFile=annotationFile)

	# note the time we start
	startTime <- proc.time()

	# setup...
	ans0 <- pipe.PreAlignTasks( sampleID, annotationFile, optionsFile, dataType="DNA-seq")

	# do the alignment
	ans1 <- pipe.DNAalignment( sampleID, annotationFile, optionsFile)

	# build the wiggle track data
	ans2 <- pipe.AlignToWig( sampleID, annotationFile, optionsFile, dataType="DNA-seq")

	# make transcriptomes for all species...
	ans3 <- pipe.Transcriptome( sampleID, annotationFile, optionsFile, dataType="DNA-seq")

	# make any extra wanted display files
	ans3 <- pipe.PostAlignTasks( sampleID, annotationFile, optionsFile, dataType="DNA-seq")

	# done...   report total time used
	cat( "\n", verboseOutputDivider, "\n\nDNA-seq PIPELINE DONE:   ", sampleID, "\n\n")
	print( elapsedProcTime( startTime, proc.time(), N=ans1$nReadsIn))

	return()
}

