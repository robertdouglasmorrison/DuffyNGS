# cleanupFiles.R -- delete unneed files from various processing steps


`cleanupFiles` <- function( path="results") {

	cat( "\nBAM file cleanup:\n")
	cleanupBAMfiles( path)

	cat( "\nVelvet file cleanup:\n")
	cleanupVelvetFiles( path=file.path( path, "VelvetContigs"))

	cat( "\nHLA Typing file cleanup:\n")
	cleanupHLAtypingFiles( path=file.path( path, "HLA.typing"))

}
