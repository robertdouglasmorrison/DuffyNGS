# checkDuffyExecutables.R -- verify needed third party software

`checkDuffyPrograms` <- function() {

	# explicit code to detect presence and version of the various tools we rely on
	Nfound <- Nfail <- 0

	cat( "\nChecking for third party software programs used by DuffyNGS pipeline tools:\n")

	# main workhorse first
	ans <- checkBowtie2Version()
	if ( is.null( ans)) {
		Nfail <- Nfail + 1
	} else {
		cat( "\nBowtie2:   version =", ans, "   path =", as.character(Sys.which("bowtie2")))
		Nfound <- Nfound + 1
	}

	# samtools
	exe <- Sys.which( "samtools")
	if ( exe == "") {
		cat( "\n'samtools' is not in the current search path..")
		Nfail <- Nfail + 1
	} else {
		ans <- system( "samtools --version", intern=T)
		if ( length(ans) == 3) {
			cat( "\nsamtools:  version =", sub("samtools ","",ans[1]), "       path =", exe)
			Nfound <- Nfound + 1
		} else {
			cat( "\nsamtools:  unexpected reply;\n", ans)
			Nfail <- Nfail + 1
		}
	}

	# bcftools
	exe <- Sys.which( "bcftools")
	if ( exe == "") {
		cat( "\n'bcftools' is not in the current search path..")
		Nfail <- Nfail + 1
	} else {
		ans <- system( "bcftools --version", intern=T)
		if ( length(ans) == 6) {
			cat( "\nbcftools:  version =", sub("bcftools ","",ans[1]), "       path =", exe)
			Nfound <- Nfound + 1
		} else {
			cat( "\nbcftools:  unexpected reply;\n", ans)
			Nfail <- Nfail + 1
		}
	}

	# python
	exe <- Sys.which( "python2")
	if ( exe == "") {
		cat( "\n'python2' is not in the current search path..")
		Nfail <- Nfail + 1
	} else {
		ans <- system( "python2 --version  2>&1", intern=T)
		#cat( "\nDebug: ", exe, "|", length(ans), "|", ans, "\n")
		if ( length(ans) == 1) {
			cat( "\npython2:   version =", sub("Python ","",ans[1]), "       path =", exe)
			Nfound <- Nfound + 1
		} else {
			cat( "\npython2:  unexpected reply;\n", ans)
			Nfail <- Nfail + 1
		}
	}

	# perl
	exe <- Sys.which( "perl")
	if ( exe == "") {
		cat( "\n'perl' is not in the current search path..")
		Nfail <- Nfail + 1
	} else {
		ans <- system( "perl --version", intern=T)
		if ( length(ans) > 1) {
			ans <- sub( " built.+", "", sub("This is perl ", "", ans[2]))
			cat( "\nperl:      version =", ans, "       path =", exe)
			Nfound <- Nfound + 1
		} else {
			cat( "\nperl:  unexpected reply;\n", ans)
			Nfail <- Nfail + 1
		}
	}

	# mafft
	exe <- Sys.which( "mafft")
	if ( exe == "") {
		cat( "\n'mafft' is not in the current search path..")
		Nfail <- Nfail + 1
	} else {
		ans <- system( "mafft --version  2>&1", intern=T)
		if ( length(ans) == 1) {
			cat( "\nmafft:     version =", ans, "    path =", exe)
			Nfound <- Nfound + 1
		} else {
			cat( "\nmafft:     unexpected reply;\n", ans)
			Nfail <- Nfail + 1
		}
	}

	# velveth
	exe <- Sys.which( "velveth")
	if ( exe == "") {
		cat( "\n'velveth' is not in the current search path..")
		Nfail <- Nfail + 1
	} else {
		ans <- system( "velveth", intern=T) # no explict 'version arg, grab 2nd line
		if ( length(ans) > 1) {
			cat( "\nvelveth:   version =", sub("Version ","",ans[2]), "    path =", exe)
			Nfound <- Nfound + 1
		} else {
			cat( "\nvelveth:   unexpected reply;\n", ans)
			Nfail <- Nfail + 1
		}
	}
	exe <- Sys.which( "velvetg")
	if ( exe == "") {
		cat( "\n'velvetg' is not in the current search path..")
		Nfail <- Nfail + 1
	} else {
		ans <- suppressWarnings( system( "velvetg", intern=T)) # no explict 'version arg, grab 2nd line
		if ( length(ans) > 1) {
			cat( "\nvelvetg:   version =", sub("Version ","",ans[2]), "    path =", exe)
			Nfound <- Nfound + 1
		} else {
			cat( "\nvelvetg:   unexpected reply;\n", ans)
			Nfail <- Nfail + 1
		}
	}

	# cutadapt
	exe <- Sys.which( "cutadapt")
	if ( exe == "") {
		cat( "\n'cutadapt' is not in the current search path..")
		Nfail <- Nfail + 1
	} else {
		ans <- system( "cutadapt --version", intern=T) 
		if ( length(ans) == 1) {
			cat( "\ncutdapt:   version =", ans[1], "      path =", exe)
			Nfound <- Nfound + 1
		} else {
			cat( "\ncutadapt:   unexpected reply;\n", ans)
			Nfail <- Nfail + 1
		}
	}

	# bowtie 1
	ans <- checkBowtieVersion()
	if ( is.null( ans)) {
		Nfail <- Nfail + 1
		cat( "\n\nNote:  Bowtie1 only used for building genome detectability files when adding new genomes...\n")
	} else {
		cat( "\nBowtie1:   version =", ans, "     path =", as.character(Sys.which("bowtie")))
		Nfound <- Nfound + 1
	}

	if ( Nfail == 0) {
		cat( "\n\nAll programs used by DuffyNGS are found and ready.\n")
		return()
	}
	if ( Nfail > 0) {
		cat( "\n\nFailled to find and verify ", Nfail, " programs wanted by DuffyNGS.\n")
		cat( "\nCheck your PATH environment variable etc...\n")
	}
	
	return( list( "N_Found"=Nfound, "N_Fail"=Nfail))
}
